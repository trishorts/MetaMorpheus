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

## Phased delivery — status
- **P0 (DONE):** phase timing. Search engine **18 CPU-s** vs post-analysis **0.1 CPU-s** (~200:1).
  Hot spot = the match/score loop. See `E:\CodeReview\metaProteomics\PROFILING_RESULTS.md`.
- **P1 (DONE):** `ISpectralScorer` + `CpuSpectralScorer` + flattened resident `SpectralScoringData`;
  wired into the engine; verified bit-identical.
- **P2 (DONE — correctness ✅, performance ❌):** ILGPU CUDA FP64 `GpuSpectralScorer` (shared-memory
  scan-slice kernel) + detector + factory. **Byte-identical end-to-end on real data.** Microbench
  (one big batch, kernel-only): **8.6x**. But the real run is **SLOWER** (see below).
- **P2.5 (DONE):** GPU gated behind `MM_PARALLELSEARCH_GPU=1`; **default is CPU** so the branch is
  not "GPU = slower" while the redesign is pending.
- **P3 (NEXT — the redesign that makes GPU actually win):** see below.
- **P4 (later):** AMD card (ILGPU OpenCL); compose with mzLib #1036 precomputed fragmentation.

## Why the integrated GPU path is slower (the lesson)
The microbench measured **one large batch, kernel-only, no lock, no per-call overhead** → 8.6x.
The real engine is the opposite shape:
- `TransientClassicSearchEngine` is constructed **per database** and runs `Parallel.ForEach`
  internally. With 12 databases concurrent × ~5 threads = **~60 CPU threads**.
- Each thread flushes batches via a **synchronous, locked** GPU round-trip: scan-sort → upload →
  launch → `Synchronize` → `GetAsArray1D` (allocates) → host compaction. One global GPU lock →
  **full serialization** + per-call overhead + lock-wait.
- First cut also created a CUDA context + re-uploaded the 45 MB spectra **per database** (fixed by
  sharing one scorer — but that *concentrated* all 60 threads on one lock, making it worse: 95 → 239 CPU-s).

Measured (12 db / 1 raw): CPU baseline 18 CPU-s (1.5 s wall) → GPU 239 CPU-s (20 s wall), byte-identical.

## P3 redesign — single GPU consumer + pipeline (required before GPU wins)
Goal: keep the fast kernel fed without 60 threads blocking on it. Decouple producers from the GPU.

1. **One GPU consumer thread (or a small fixed pool), not 60 producers calling the GPU.**
   - CPU search threads do digestion + fragmentation + candidate-scan selection and **enqueue**
     work items into a concurrent queue/ring buffer. They never touch the GPU and never block on it.
   - A single consumer pulls work, packs **large coalesced batches** (tens–hundreds of k work items),
     launches the kernel, and dispatches results back. This removes lock contention entirely and
     amortizes per-call overhead over huge batches.

2. **Pipeline / double-buffer with CUDA streams.** Overlap H2D upload of batch N+1, kernel of
   batch N, and D2H download of batch N-1 (ILGPU streams + 2–3 rotating buffer sets). The GPU never
   idles waiting on transfers; the bench's "kernel-only" number becomes the realistic number.

3. **Result routing.** Each work item carries its `(scanIndex, notch, peptide, products)` identity;
   results return asynchronously, so the consumer (or a completion callback) must apply
   `AddPeptideCandidateToPsm` for the right PSM. Preserve per-scan candidate ordering for
   byte-identity — easiest is to keep, per producing thread, the original work order and apply
   updates in that order once a batch completes (single-threaded apply, or per-scan locks as today).

4. **Cheap wins to fold in:** reuse a host output buffer via `CopyToCPU` (stop allocating with
   `GetAsArray1D` every batch); drop or cheapen the per-batch scan-sort (real data is already
   roughly scan-clustered; or sort once per large batch); raise/auto-tune the flush threshold so
   batches are big; size `MaxSlice` from the real per-scan fragment-count distribution (fall back to
   global search above it).

5. **Threading-model interaction with the engine.** The current `Parallel.ForEach`-per-database +
   `Parallel.ForEach`-over-proteins nesting is the source of the 60-way contention. The consumer
   model wants: many CPU producers (fine to keep the nested parallelism for digestion) but a single
   GPU drain. Consider building the producer/consumer at the `ParallelSearchTask` level (across
   databases) so one consumer serves all databases, rather than per-engine.

6. **Re-validate** byte-identity (run both baseline and GPU **single-threaded**, `MaxThreadsToUsePerFile=1`,
   for a guaranteed-deterministic diff) and re-measure wall-clock. Target: approach the ~6–9x
   kernel/upload microbench numbers on the dominant transient-search cost.

### Amdahl + future factor
- `Initialize` (~37 s one-time: base human search + decon + DB load) is **not** GPU-accelerated;
  overall run speedup is bounded by it on small runs but it amortizes at 30k-DB scale.
- **mzLib #1036 (binary-indexed spectral library)** would precompute theoretical fragmentation and
  populate `ScoringBatch.FragmentNeutralMasses` directly (no per-peptide `Fragment()`), shrinking the
  CPU producer cost and shifting the producer/consumer balance — design P3 with this in mind.

## How to run the GPU path (for P3 work)
`set MM_PARALLELSEARCH_GPU=1` (PowerShell: `$env:MM_PARALLELSEARCH_GPU=1`) before the CMD run.
Standalone kernel harness for fast iteration: `E:\CodeReview\metaProteomics\GpuMicrobench`.
