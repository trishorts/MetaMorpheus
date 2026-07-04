# ISD search reproduction kit

Everything needed to reproduce the ISD-vs-DDA proteoform-ID searches from this branch: the exact search/GPTMD
configs (`tomls/`), parameterized scripts (`scripts/`), and how the pieces fit together. Results these produce
are written up in the manuscript repo `collab_zs/analysis/RESULTS_isd_vs_dda_gptmd.md` and
`RESULTS_consensus_xic_ids.md`.

## DEFAULT DDA strategy = the consensus-paper way
**We search DDA by default with precursor feature tracing + FromFile assembly** (consensus-trace the MS1 →
`_ms1.feature` → `ms2.GetIsolatedMassesAndCharges(ms1, FromFileDeconvolutionParameters)` → real fragments →
MGF → search). Script `05_dda_fromfile_search.sh`; in-code `SearchDdaWithFromFileConsensusFeatures_Complete_FromEnv`.
This makes DDA and ISD share the **same** search config (`td_pseudoMS2_search.toml`, provided precursors) — the
only difference between them is how the fragments are obtained (consensus-traced for ISD vs the real isolated MS2
for DDA). **Normal live-deconvolution DDA (`02_dda_search.sh`, `dda_topdown_search.toml`) is kept only as a
labeled baseline.**

## The method, end to end
Two stages: **(A) generate pseudo-MS2 spectra** (code in this branch), **(B) search them** (MetaMorpheus CMD).

**A. Pseudo-MS2 generation (C#, in the Test project — env-driven harnesses in
`MetaMorpheus/Test/DIATests/IsdMsAlignExportTest.cs`):**
- `ConsensusIdSweep_FromEnv` — ISD path. Reads an ISD raw/mzML, splits scans by source voltage (`sid=`),
  builds XICs with **consensus feature tracing** (`ConsensusMassXicConstructor` = mzLib MassTraceBuilder →
  TraceCorrector → MassFeatureBuilder), groups fragments to precursors by **apex-RT + Pearson correlation**
  (`XicGroupingEngine`), builds pseudo-MS2 (`ISDEngine`/`GetPseudoMs2ScanFromPfGroup`), and writes one MGF per
  correlation cutoff. Env knobs: `ISD_MZML`, `ISD_OUTDIR`, `ISD_CORRS` (e.g. `0,0.4,0.6,0.7`),
  `ISD_PREC_MINMASS`/`ISD_PREC_MINCHARGE` (precursor intact filter), `ISD_FRAG_AGG` (`false` = fragment mass
  tracing only, no charge aggregation).
- `SearchDdaWithConsensusPrecursors_FromEnv` — DDA-with-consensus-precursors path (hand-rolled nearest-feature
  match). Consensus-traces the DDA MS1, pairs each real DDA MS2 to its consensus precursor (isolation m/z + RT),
  keeps the real fragments, writes an MGF. Env knobs: `DDA_MZML`, `DDA_OUT`.
- `SearchDdaWithFromFileConsensusFeatures_FromEnv` — **the consensus-paper DDA method.** Consensus-traces the
  DDA MS1, writes an external `*_ms1.feature` file, and assembles each MS2's precursor via mzLib's real FromFile
  join `ms2.GetIsolatedMassesAndCharges(ms1, FromFileDeconvolutionParameters)` (matching features to the
  isolation window + precursor RT — the same mechanism MetaMorpheus uses), then attaches real fragments -> MGF.
  Env knobs: `DDAFF_MZML`, `DDAFF_OUT`. (Script `05_dda_fromfile_search.sh`.)
- `SearchDdaWithFromFileConsensusFeatures_Complete_FromEnv` — **fully self-contained**: does the above AND runs a
  real top-down `SearchTask` (config built in C#, no external toml/CMD) against a supplied protein DB, then
  reports proteoform + PSM counts at 1% FDR and writes the normal MetaMorpheus output (`AllProteoforms.psmtsv`,
  `AllPSMs.psmtsv`, `results.txt`). Env knobs: `DDAFF_MZML`, `DDAFF_DB`, `DDAFF_OUT`. Run directly with:
  `DDAFF_MZML=... DDAFF_DB=... DDAFF_OUT=... dotnet test Test.csproj --filter FullyQualifiedName~SearchDdaWithFromFileConsensusFeatures_Complete_FromEnv`

The pseudo-MS2 is written as MGF by `IsdMsAlignExporter.WriteMgf` (fragments as neutral-mass+proton). It can
also be written as TopPIC `.msalign` via `IsdMsAlignExporter.WriteMs2Align`/`WriteMs1Align`.

**ISD and DDA assemble the MS2 scan differently — on purpose.** The PRECURSOR is built the same way in both
(a consensus-traced MS1 feature → m/z at its charge). The FRAGMENTS differ by necessity: ISD fragments are
consensus-traced across the all-MS1 voltage cycle (they are spread over multiple scans), while DDA has one real
isolated MS2 scan per precursor and so deconvolutes that scan directly. Do not unify these.

**B. Search (MetaMorpheus CMD):** a normal top-down `SearchTask` on the MGF (ISD paths) or the real `.mzML`
(normal DDA), using the tomls in `tomls/`. Both ISD and DDA use **Classic** deconvolution so the comparison is
fair. IDs are counted from `AllProteoforms.psmtsv` at `QValue <= 0.01`.

## Configs (`tomls/`)
| file | used for | key settings |
|---|---|---|
| `td_pseudoMS2_search.toml` | searching ISD / DDA-consensus **MGFs** | Classic, top-down, `UseProvidedPrecursorInfo=true`, product charge 1 |
| `dda_topdown_search.toml` | searching **normal DDA** `.mzML` | Classic, top-down, `DoPrecursorDeconvolution=true` |
| `gptmd_isd.toml` | GPTMD on ISD MGFs | GPTMD mod list + the ISD CommonParameters |
| `gptmd_dda.toml` | GPTMD on DDA mzML | GPTMD mod list + the DDA CommonParameters |

## Build once
```
dotnet build MetaMorpheus/MetaMorpheus/CMD/CMD.csproj   -c Release
dotnet build MetaMorpheus/MetaMorpheus/Test/Test.csproj -c Debug
```

## Run (edit paths in `scripts/config.sh` first)
```bash
cd isd-search-reproduction/scripts

# --- ISD (consensus tracing) ---
./01_isd_consensus_search.sh  "$RAWDIR/09-18-25_YC_81min_ISD60-80-100_preFilter700-900-1100_rep1.raw"  YC  "0,0.4,0.6,0.7"  true

# --- DDA: DEFAULT = consensus-paper FromFile method ---
./05_dda_fromfile_search.sh   "$RAWDIR/09-18-25_YC_81min_DDA_1-5iso_mscan4_rep1.raw"  YC
# GPTMD on the default DDA (FromFile MGF) uses the ISD/MGF search config:
./04_gptmd.sh  isd  "$OUT/YC/dda_fromfile.mgf"             YC_dda

# --- DDA baselines (for comparison only) ---
./02_dda_search.sh            "$RAWDIR/09-18-25_YC_81min_DDA_1-5iso_mscan4_rep1.raw"  YC   # live deconvolution
./03_dda_consensus_search.sh  "$RAWDIR/09-18-25_YC_81min_DDA_1-5iso_mscan4_rep1.raw"  YC   # hand-rolled match

# --- ISD + GPTMD ---
./04_gptmd.sh  isd  "$OUT/YC/consensus_corr0.mgf"          YC
```
Each script prints `RESULT ... (total ptm)= <N> <M>` = proteoforms and PTM-bearing proteoforms at 1% FDR.

## Fractions used (Test_runs/09-17-25)
YC/YD/YE, rep1 — each has a matched DDA + ISD (preFilter) run. Same yeast UniProt DB.

## Notes
- ISD reads `.raw` directly in the test harness; the CMD DDA search needs `.mzML` (the `.raw` path triggers an
  interactive Thermo-license prompt that crashes headless), so `02`/`03` convert first with ThermoRawFileParser.
- Absolute counts are Classic-deconvolution-limited (below the FLASHDeconv pipeline). The ISD-vs-DDA and
  ±GPTMD *comparisons* are the point.
