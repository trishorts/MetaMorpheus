using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test.DIATests
{
    [TestFixture]
    public static class IsdMsAlignExportTest
    {
        /// <summary>
        /// End-to-end: take a real ISD mzML, build pseudo-MS2 scans using CONSENSUS FEATURE TRACING
        /// (<see cref="ConsensusMassXicConstructor"/> = mzLib MassTraceBuilder→TraceCorrector→MassFeatureBuilder)
        /// as the XIC source, combined with apex-RT/Pearson precursor-fragment grouping and pseudo-MS2 generation,
        /// then export them to ms1/ms2.msalign and confirm the ms2.msalign round-trips through the MsAlign reader
        /// (i.e. it is a valid file TopPIC could consume).
        /// </summary>
        [Test]
        public static void ExportIsdPseudoScansToMsAlign()
        {
            string mzmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "DIA",
                "08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_RT8.1-9.23.mzML");
            Assert.That(File.Exists(mzmlPath), $"missing test fixture: {mzmlPath}");

            var dataFile = MsDataFileReader.GetDataFile(mzmlPath);
            dataFile.LoadAllStaticData();

            // Consensus feature tracing as the XIC front-end for both the intact (precursor) and fragment channels.
            var diaParams = new DIAparameters(
                DIAanalysisType.ISD,
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)),
                new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)),
                new XicGroupingEngine(0.2f, 0.5, 0.5, 1),
                PseudoMs2ConstructionType.Mass,
                combineFragments: true);
            var commonParams = new CommonParameters { DIAparameters = diaParams };

            var pseudoScans = MetaMorpheusTask.GetMs2Scans(dataFile, mzmlPath, commonParams).ToArray();
            Assert.That(pseudoScans.Length, Is.GreaterThan(0), "no pseudo scans produced from the ISD fixture");

            string ms2Path = Path.Combine(TestContext.CurrentContext.TestDirectory, "isd_pseudo_ms2.msalign");
            string ms1Path = Path.Combine(TestContext.CurrentContext.TestDirectory, "isd_pseudo_ms1.msalign");
            IsdMsAlignExporter.WriteMs2Align(pseudoScans, ms2Path);
            IsdMsAlignExporter.WriteMs1Align(pseudoScans, ms1Path);
            Assert.That(File.Exists(ms2Path));
            Assert.That(File.Exists(ms1Path));

            // number of BEGIN IONS blocks == number of pseudo scans
            int blocks = File.ReadLines(ms2Path).Count(l => l.StartsWith("BEGIN IONS"));
            Assert.That(blocks, Is.EqualTo(pseudoScans.Length));

            // round-trip through the MsAlign reader: it must parse back the same number of scans
            var reread = new Ms2Align(ms2Path);
            reread.LoadAllStaticData();
            Assert.That(reread.Scans.Count(s => s != null), Is.EqualTo(pseudoScans.Length));
            // and every reread scan is MS2 with a precursor mass carried through
            Assert.That(reread.Scans.Where(s => s != null).All(s => s.MsnOrder == 2));

            File.Delete(ms2Path);
            File.Delete(ms1Path);
        }

        /// <summary>
        /// Environment-driven harness: generate consensus-path pseudo-MS2 MGFs at several precursor-fragment
        /// correlation thresholds, reusing one consensus-XIC build. Set ISD_MZML (input), ISD_OUTDIR (output),
        /// and optionally ISD_CORRS (comma list, default "0,0.4,0.6,0.7"). Writes consensus_corr&lt;t&gt;.mgf per
        /// threshold; each MGF is then searched externally with the top-down td_pseudoMS2 config to count IDs.
        /// Ignored unless ISD_MZML is set, so it never runs in normal CI.
        /// </summary>
        [Test]
        public static void ConsensusIdSweep_FromEnv()
        {
            string mzml = Environment.GetEnvironmentVariable("ISD_MZML");
            if (string.IsNullOrEmpty(mzml)) { Assert.Ignore("set ISD_MZML to run this harness"); return; }
            string outDir = Environment.GetEnvironmentVariable("ISD_OUTDIR") ?? TestContext.CurrentContext.TestDirectory;
            Directory.CreateDirectory(outDir);
            var corrs = (Environment.GetEnvironmentVariable("ISD_CORRS") ?? "0,0.4,0.6,0.7")
                .Split(',').Select(s => double.Parse(s, System.Globalization.CultureInfo.InvariantCulture)).ToArray();
            float apex = 0.2f;
            double overlap = 0.5;
            var decon = new ClassicDeconvolutionParameters(1, 30, 4, 3);

            var dataFile = MsDataFileReader.GetDataFile(mzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var isdMap = ISDEngine.ConstructIsdGroups(allScans, out var ms1Scans);
            ISDEngine.ReLabelIsdScans(isdMap, allScans);

            // consensus feature tracing -> XICs, built ONCE and reused across thresholds.
            // PRECURSOR channel: keep only intact proteoforms (>=3 kDa, >=3 charge states) so small noise
            // features do not become precursors. FRAGMENT channel: unfiltered (fragments are low mass/charge).
            double precMinMass = double.Parse(Environment.GetEnvironmentVariable("ISD_PREC_MINMASS") ?? "3000",
                System.Globalization.CultureInfo.InvariantCulture);
            int precMinCharge = int.Parse(Environment.GetEnvironmentVariable("ISD_PREC_MINCHARGE") ?? "3");
            var ms1Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon,
                minMass: precMinMass, minChargeCount: precMinCharge);
            var ms2Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon);
            // ISD_FRAG_AGG=false -> fragments use mass tracing only (one XIC per charge-locked trace, no
            // cross-charge aggregation); precursors always use full feature tracing.
            bool fragAgg = (Environment.GetEnvironmentVariable("ISD_FRAG_AGG") ?? "true").ToLowerInvariant() != "false";
            ms2Ctor.AggregateCharges = fragAgg;
            var ms1Xics = ms1Ctor.GetAllXics(ms1Scans);
            var ms2Xics = isdMap.Values.SelectMany(v => ms2Ctor.GetAllXics(v.ToArray())).ToList();
            TestContext.WriteLine($"consensus XICs: {ms1Xics.Count} precursor (>= {precMinMass} Da, >= {precMinCharge} charges), {ms2Xics.Count} fragment (chargeAggregation={fragAgg})");

            foreach (var corr in corrs)
            {
                var grouper = new XicGroupingEngine(apex, overlap, corr, Environment.ProcessorCount);
                var groups = grouper.PrecursorFragmentGrouping(ms1Xics, ms2Xics).ToList();
                int idx = 0;
                var scans = new List<Ms2ScanWithSpecificMass>();
                foreach (var g in groups)
                {
                    idx++;
                    g.PFgroupIndex = idx;
                    scans.Add(g.GetPseudoMs2ScanFromPfGroup(PseudoMs2ConstructionType.Mass, new CommonParameters(), mzml));
                }
                string outMgf = Path.Combine(outDir, $"consensus_corr{corr.ToString(System.Globalization.CultureInfo.InvariantCulture)}.mgf");
                IsdMsAlignExporter.WriteMgf(scans, outMgf);
                TestContext.WriteLine($"corr={corr}: {scans.Count} pseudo scans -> {outMgf}");
            }
        }

        /// <summary>
        /// Search DDA runs using CONSENSUS MASS TRACING for the precursors: trace the DDA MS1 scans with
        /// <see cref="ConsensusMassXicConstructor"/> to get clean intact precursor features, pair each real DDA
        /// MS2 scan to the consensus feature matching its isolation m/z + retention time, and build a pseudo scan
        /// from the consensus precursor + the scan's REAL (deconvoluted) fragments. Writes an MGF that is then
        /// searched with the same td_pseudoMS2 config as the ISD path (apples-to-apples ISD-consensus vs
        /// DDA-consensus). Env-driven: DDA_MZML (input), DDA_OUT (output MGF). Ignored unless DDA_MZML is set.
        ///
        /// NOTE (intentional design — do not "unify"): ISD and DDA assemble the Ms2ScanWithSpecificMass
        /// PRECURSOR identically (a consensus-traced MS1 feature -> m/z at its charge), but assemble the
        /// FRAGMENTS differently, by necessity. ISD fragments are synthesized from consensus-traced
        /// fragmentation-channel features (they are spread across the all-MS1 voltage cycle, so they must be
        /// traced across scans). DDA has a single real isolated MS2 scan per precursor -> nothing to trace, so
        /// its fragments come from deconvoluting that real scan (GetNeutralExperimentalFragments). The two paths
        /// SHOULD differ on the fragment side.
        /// </summary>
        [Test]
        public static void SearchDdaWithConsensusPrecursors_FromEnv()
        {
            string mzml = Environment.GetEnvironmentVariable("DDA_MZML");
            if (string.IsNullOrEmpty(mzml)) { Assert.Ignore("set DDA_MZML to run this harness"); return; }
            string outMgf = Environment.GetEnvironmentVariable("DDA_OUT")
                ?? Path.Combine(TestContext.CurrentContext.TestDirectory, "dda_consensus.mgf");
            const double proton = 1.007276466879;
            var decon = new ClassicDeconvolutionParameters(1, 60, 4, 3);
            // real DDA MS2 fragments are multiply charged -> deconvolute them to neutral masses
            var fragParams = new CommonParameters(productDeconParams: new ClassicDeconvolutionParameters(1, 20, 4, 3));

            var dataFile = MsDataFileReader.GetDataFile(mzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var ms1Scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
            var ms2Scans = allScans.Where(s => s.MsnOrder == 2).ToArray();

            // consensus mass tracing on the MS1 scans -> intact precursor features (>=3 kDa, >=3 charges)
            var ms1Ctor = new ConsensusMassXicConstructor(new PpmTolerance(20), 2, 1, 3, decon,
                minMass: 3000, minChargeCount: 3);
            var precursorXics = ms1Ctor.GetAllXics(ms1Scans);

            // pair each real DDA MS2 scan with the consensus feature matching its isolation m/z + RT
            var pseudoScans = new List<Ms2ScanWithSpecificMass>();
            int assigned = 0;
            foreach (var ms2 in ms2Scans)
            {
                if (ms2.IsolationMz == null) continue;
                double isoMz = ms2.IsolationMz.Value;
                double rt = ms2.RetentionTime;
                ExtractedIonChromatogram best = null; double bestErr = double.MaxValue; int bestZ = 1;
                foreach (var xic in precursorXics)
                {
                    if (rt < xic.StartRT - 0.5 || rt > xic.EndRT + 0.5) continue;
                    double mass = xic.ApexPeak.M;
                    int z = (int)Math.Round(mass / (isoMz - proton)); // charge that places this mass at the iso m/z
                    if (z < 1 || z > 60) continue;
                    double err = Math.Abs(mass.ToMz(z) - isoMz);
                    if (err < bestErr && err < 2.5) { bestErr = err; best = xic; bestZ = z; }
                }
                if (best == null) continue;
                assigned++;
                double precMz = ((double)best.ApexPeak.M).ToMz(bestZ);
                var frags = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2, fragParams);
                pseudoScans.Add(new Ms2ScanWithSpecificMass(ms2, precMz, bestZ, mzml, fragParams, frags));
            }
            IsdMsAlignExporter.WriteMgf(pseudoScans, outMgf);
            TestContext.WriteLine($"DDA consensus: {precursorXics.Count} MS1 features, {ms2Scans.Length} MS2, {assigned} assigned -> {outMgf}");
            Assert.That(assigned, Is.GreaterThan(0), "no MS2 scans matched a consensus precursor feature");
        }

        /// <summary>
        /// Search DDA the CONSENSUS-PAPER way: precursor feature tracing + the built-in "FromFile" precursor
        /// assembly. Consensus-trace the DDA MS1 into features, write them as a `.ms1.feature` file, then use
        /// mzLib's real <c>FromFileDeconvolutionParameters</c> join — `ms2.GetIsolatedMassesAndCharges(ms1, ff)`
        /// — to assemble each MS2's precursor exactly as MetaMorpheus's precursor pipeline does (matching
        /// features to the isolation window + precursor RT), then attach the scan's real fragments. This is the
        /// mechanism from the consensus-deconvolution paper's ExternalMs1Features / FromFile search, as opposed
        /// to the hand-rolled nearest-feature matching in SearchDdaWithConsensusPrecursors_FromEnv.
        /// Env-driven: DDAFF_MZML (input), DDAFF_OUT (output MGF). Ignored unless DDAFF_MZML is set.
        /// </summary>
        [Test]
        public static void SearchDdaWithFromFileConsensusFeatures_FromEnv()
        {
            string mzml = Environment.GetEnvironmentVariable("DDAFF_MZML");
            if (string.IsNullOrEmpty(mzml)) { Assert.Ignore("set DDAFF_MZML to run this harness"); return; }
            string outMgf = Environment.GetEnvironmentVariable("DDAFF_OUT")
                ?? Path.Combine(TestContext.CurrentContext.TestDirectory, "dda_fromfile.mgf");
            string featPath = Path.Combine(Path.GetDirectoryName(outMgf) ?? ".",
                Path.GetFileNameWithoutExtension(outMgf) + "_ms1.feature"); // reader keys on the "_ms1.feature" suffix
            var decon = new ClassicDeconvolutionParameters(1, 60, 4, 3);
            var fragParams = new CommonParameters(productDeconParams: new ClassicDeconvolutionParameters(1, 20, 4, 3));

            var dataFile = MsDataFileReader.GetDataFile(mzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var ms1Scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
            var ms2Scans = allScans.Where(s => s.MsnOrder == 2).ToArray();

            // (A) precursor feature tracing (consensus) -> external _ms1.feature file (intact precursors only)
            var features = ConsensusMassXicConstructor.TraceFeatures(ms1Scans, decon);
            IsdMsAlignExporter.WriteMs1FeatureFile(features, featPath, minMass: 3000, minChargeCount: 3);

            // (B) assemble Ms2ScanWithSpecificMass via mzLib's FromFile join (the consensus-paper mechanism)
            var ff = new Readers.FromFileDeconvolutionParameters(featPath, 1, 60);
            var pseudoScans = new List<Ms2ScanWithSpecificMass>();
            int assigned = 0;
            foreach (var ms2 in ms2Scans)
            {
                if (!ms2.OneBasedPrecursorScanNumber.HasValue) continue;
                var ms1 = dataFile.GetOneBasedScan(ms2.OneBasedPrecursorScanNumber.Value);
                IsotopicEnvelope[] frags = null;
                foreach (var env in ms2.GetIsolatedMassesAndCharges(ms1, ff))
                {
                    frags ??= Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2, fragParams);
                    assigned++;
                    double precMz = env.MonoisotopicMass.ToMz(env.Charge);
                    pseudoScans.Add(new Ms2ScanWithSpecificMass(ms2, precMz, env.Charge, mzml, fragParams, frags,
                        env.TotalIntensity, env.Peaks.Count));
                }
            }
            IsdMsAlignExporter.WriteMgf(pseudoScans, outMgf);
            TestContext.WriteLine($"FromFile consensus: {features.Count} features, {assigned} precursor assignments over {ms2Scans.Length} MS2 -> {outMgf}");
            Assert.That(assigned, Is.GreaterThan(0), "no MS2 scans got a FromFile precursor");
        }

        /// <summary>
        /// COMPLETE, self-contained DDA search via the consensus-paper method — no external toml or CMD needed.
        /// Does the whole pipeline in-process: consensus-trace the DDA MS1 → write a `_ms1.feature` → assemble
        /// each MS2's precursor with mzLib's FromFile join → write an MGF → run a real MetaMorpheus top-down
        /// <see cref="SearchTask"/> (config built in code below) against a protein database → report proteoform
        /// and PSM counts at 1% FDR, exactly like a normal MetaMorpheus search.
        /// Env-driven: DDAFF_MZML (DDA input), DDAFF_DB (protein .xml/.fasta), DDAFF_OUT (output dir, optional).
        /// Ignored unless both DDAFF_MZML and DDAFF_DB are set.
        /// </summary>
        [Test]
        public static void SearchDdaWithFromFileConsensusFeatures_Complete_FromEnv()
        {
            string mzml = Environment.GetEnvironmentVariable("DDAFF_MZML");
            string db = Environment.GetEnvironmentVariable("DDAFF_DB");
            if (string.IsNullOrEmpty(mzml) || string.IsNullOrEmpty(db))
            { Assert.Ignore("set DDAFF_MZML and DDAFF_DB to run the complete in-process search"); return; }
            string outDir = Environment.GetEnvironmentVariable("DDAFF_OUT")
                ?? Path.Combine(TestContext.CurrentContext.TestDirectory, "DdaFromFileSearch");
            Directory.CreateDirectory(outDir);
            string outMgf = Path.Combine(outDir, "dda_fromfile.mgf");
            string featPath = Path.Combine(outDir, "dda_fromfile_ms1.feature");

            // ---- (1) precursor feature tracing + FromFile assembly -> MGF ----
            var decon = new ClassicDeconvolutionParameters(1, 60, 4, 3);
            var fragParams = new CommonParameters(productDeconParams: new ClassicDeconvolutionParameters(1, 20, 4, 3));
            var dataFile = MsDataFileReader.GetDataFile(mzml);
            dataFile.LoadAllStaticData();
            var allScans = dataFile.GetAllScansList().ToArray();
            var ms1Scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
            var ms2Scans = allScans.Where(s => s.MsnOrder == 2).ToArray();
            var features = ConsensusMassXicConstructor.TraceFeatures(ms1Scans, decon);
            IsdMsAlignExporter.WriteMs1FeatureFile(features, featPath, minMass: 3000, minChargeCount: 3);
            var ff = new FromFileDeconvolutionParameters(featPath, 1, 60);
            var pseudoScans = new List<Ms2ScanWithSpecificMass>();
            foreach (var ms2 in ms2Scans)
            {
                if (!ms2.OneBasedPrecursorScanNumber.HasValue) continue;
                var ms1 = dataFile.GetOneBasedScan(ms2.OneBasedPrecursorScanNumber.Value);
                IsotopicEnvelope[] frags = null;
                foreach (var env in ms2.GetIsolatedMassesAndCharges(ms1, ff))
                {
                    frags ??= Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2, fragParams);
                    pseudoScans.Add(new Ms2ScanWithSpecificMass(ms2, env.MonoisotopicMass.ToMz(env.Charge),
                        env.Charge, mzml, fragParams, frags, env.TotalIntensity, env.Peaks.Count));
                }
            }
            IsdMsAlignExporter.WriteMgf(pseudoScans, outMgf);

            // ---- (2) a complete top-down MetaMorpheus search (config in code — the "toml" as C#) ----
            var searchTask = new SearchTask
            {
                CommonParameters = new CommonParameters(
                    digestionParams: new DigestionParams(protease: "top-down"),
                    doPrecursorDeconvolution: false,      // precursor is provided in the MGF (FromFile assembly)
                    useProvidedPrecursorInfo: true,
                    deconvolutionMaxAssumedChargeState: 60,
                    productMassTolerance: new PpmTolerance(20),
                    productDeconParams: new ClassicDeconvolutionParameters(1, 1, 4, 3), // MGF fragments are neutral (z=1)
                    // pseudo-MS2 peaks are all real deconvoluted fragments -> do NOT trim them
                    // (matches td_pseudoMS2.toml; CommonParameters defaults trimMsMsPeaks=true otherwise)
                    trimMs1Peaks: false,
                    trimMsMsPeaks: false),
                SearchParameters = new SearchParameters
                {
                    SearchType = SearchType.Classic,
                    MassDiffAcceptorType = MassDiffAcceptorType.PlusOrMinusThreeMM,
                    DecoyType = DecoyType.Reverse,
                    DoLabelFreeQuantification = false,
                    WriteMzId = false,
                }
            };
            searchTask.RunTask(outDir, new List<DbForTask> { new DbForTask(db, false) },
                new List<string> { outMgf }, "DdaFromFileConsensusSearch");

            // ---- (3) report proteoform + PSM counts at 1% FDR, like a normal MetaMorpheus search ----
            int proteoforms = CountAtFdr(Directory.GetFiles(outDir, "AllProteoforms.psmtsv", SearchOption.AllDirectories).FirstOrDefault());
            int psms = CountAtFdr(Directory.GetFiles(outDir, "AllPSMs.psmtsv", SearchOption.AllDirectories).FirstOrDefault());
            TestContext.WriteLine($"COMPLETE DDA FromFile-consensus search: {pseudoScans.Count} pseudo-MS2 | " +
                                  $"{proteoforms} proteoforms, {psms} PSMs at 1% FDR | features={features.Count}");
            Assert.That(proteoforms, Is.GreaterThanOrEqualTo(0));
        }

        /// <summary>Count rows at QValue &lt;= q that are Target/None (not Decoy/Contaminant) in a MetaMorpheus psmtsv.</summary>
        private static int CountAtFdr(string psmtsv, double q = 0.01)
        {
            if (psmtsv == null || !File.Exists(psmtsv)) return -1;
            var lines = File.ReadAllLines(psmtsv);
            if (lines.Length < 2) return 0;
            var header = lines[0].Split('\t');
            int qi = Array.IndexOf(header, "QValue");
            int ti = Array.IndexOf(header, "Decoy/Contaminant/Target");
            int c = 0;
            for (int i = 1; i < lines.Length; i++)
            {
                var f = lines[i].Split('\t');
                if (qi >= 0 && qi < f.Length && double.TryParse(f[qi], out var qv) && qv <= q &&
                    (ti < 0 || ti >= f.Length || f[ti] == "T" || f[ti] == "N"))
                    c++;
            }
            return c;
        }
    }
}
