# GPU acceleration plan — ParallelSearch / ManySearchTask

Branch: `ManySearchTask-gpu` (off `nbollis:ManySearchTask`).
Goal: offload the hot path of the many-database search to the NVIDIA Quadro M4000
(CUDA via **ILGPU**), keeping a CPU path for correctness comparison and fallback.

## Hardware on the dev box
- NVIDIA Quadro M4000 — Maxwell, 8 GB, CUDA cap 5.2, ~2.6 TFLOPS FP32. **Primary target.**
- AMD Radeon PRO W6300 — RDNA2, ~4 GB, OpenCL only. Out of scope for v1 (ILGPU/OpenCL
  can add it later with no kernel rewrite).
- 60+ CPU threads / 512 GB. The GPU does **not** dominate the CPU on raw FLOPS, so the win
  comes from (a) restructuring the inner loop to batched struct-of-arrays kernels and
  (b) keeping the experimental spectra resident in VRAM across all 30k databases.

Decisions locked with the user:
- NVIDIA only for v1.
- Build **comparison-first**: CPU and GPU behind one interface + a `preferCpu` toggle +
  a diff harness, so we can measure speedup *and* verify the GPU result against CPU before
  trusting it. Bit-identical vs statistically-equivalent is TBD — the harness decides.
- Profile before committing to the big port (hot spot not yet pinned: search vs statistics).

## Where the time goes (the two compute regimes)
1. **Search** — `EngineLayer/ParallelSearch/TransientClassicSearchEngine.cs`.
   Per database: digest → fragment → for each theoretical fragment,
   `Ms2ScanWithSpecificMass.GetClosestExperimentalIsotopicEnvelope` (binary search +
   tolerance + charge filter) → `score += 1 + intensity/TIC`.
   Runs across ~30k databases × ~10 raw files. Canonical GPU target.
2. **Statistics** — `TaskLayer/ParallelSearch/Statistics/*` (permutation tests, isolation
   forest, multiple-testing correction). Self-contained, embarrassingly parallel, lower risk.

Phase timing (added in this branch) logs the real split at the end of a run.

## The key structural fact that makes the GPU sing
The experimental spectra (`ParallelSearchTask.AllSortedMs2Scans`) are loaded **once** and are
**read-only and identical across all 30k transient searches** for a given file set. Each
`Ms2ScanWithSpecificMass` already exposes a **sorted `double[]` of deconvoluted fragment
masses** (`DeconvolutedMonoisotopicMasses`) plus parallel charge/intensity and per-scan
`TotalIonCurrent` / `PrecursorCharge`.

So: upload the spectra to VRAM **once** (CSR layout — concatenated sorted mass slices +
per-scan offsets), then for each database stream only the theoretical fragments. PCIe
transfer is amortized over 30k databases instead of paid per search.

VRAM budget: one file-set's deconvoluted peak lists are well under 1 GB of float32 — fits
the M4000's 8 GB with room for theoretical-fragment batches.

## Proposed architecture (mirrors mzLib PR #1027)
Zero compile-time dependency on ILGPU in the core assembly; GPU type reflection-loaded;
auto CPU fallback. New code under `EngineLayer/ParallelSearch/Gpu/`:

- `ISpectralScorer` — batched interface. Input: candidate `(scanIndex, fragment-mass slice)`
  pairs for a database; output: per-candidate `score` + `matchedCount` (objects materialized
  on CPU only for surviving hits).
- `CpuSpectralScorer` — current loop, extracted verbatim (the comparison baseline).
- `GpuSpectralScorer` — ILGPU CUDA kernel: resident spectra arrays + streamed fragment
  batches; per-thread does binary-search-in-slice + tolerance/charge test + score reduce.
- `GpuDeviceDetector` / `SpectralScorerFactory` — copied pattern from #1027 (detect, pick
  backend, `preferCpu`, fallback, `DescribeBackend`).
- `ScorerComparison` harness — run both on the same input, diff scores/match-counts per scan,
  report max/mean abs deviation and any PSM-call disagreements.

`TransientClassicSearchEngine` keeps digestion+fragmentation on the CPU (chemistry, branchy)
and delegates the match/score inner loop to the injected `ISpectralScorer`.

## Phased delivery
- **P0 (this commit):** phase timing in `ParallelSearchTask` → reveals search-vs-stats split.
- **P1:** extract `ISpectralScorer` + `CpuSpectralScorer`; wire into the engine; no behavior
  change. Add the comparison harness skeleton.
- **P2:** ILGPU CUDA `GpuSpectralScorer` (resident spectra, batched fragments) + detector +
  factory. Validate against CPU with the harness on real data.
- **P3:** benchmark CPU vs GPU on a representative subset; tune batch size; decide the
  correctness bar from the harness numbers.
- **P4 (optional):** add the AMD card (ILGPU OpenCL) and/or offload the permutation tests.
