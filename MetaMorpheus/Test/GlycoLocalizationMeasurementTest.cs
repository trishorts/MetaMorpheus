using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Measurement harness for the site-level FLR work. NOT a correctness test and NOT part of CI:
    /// every test here is [Explicit] and exists to produce numbers that decide a design question.
    ///
    /// Question being answered (PLAN_GLYCO_FLR Q1): LocalizationGraph.GetAllPaths_CalP enumerates every
    /// legal path through the DAG with no cap. Lifting the Level1/Level2 gate on probability computation
    /// means running it on the ambiguous GSMs, which are exactly the ones with the most routes. Is that a
    /// real problem or a theoretical one? Measure before choosing between a path budget and a per-site
    /// marginalization.
    /// </summary>
    [TestFixture]
    public static class GlycoLocalizationMeasurementTest
    {
        private static GlycanBox[] OGlycanBoxes { get; set; }

        [OneTimeSetUp]
        public static void Setup()
        {
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(
                GlobalVariables.OGlycanDatabasePaths.First(p => p.Contains("OGlycan.gdb")), true, true).ToArray();
            OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(3).OrderBy(p => p.Mass).ToArray();
        }

        /// <summary>
        /// Builds a tryptic peptide carrying exactly <paramref name="serineCount"/> candidate O-glycosites,
        /// spaced by alanines so every site is separated by a backbone cleavage position.
        /// </summary>
        private static PeptideWithSetModifications MakePeptideWithSites(int serineCount)
        {
            var sequence = new StringBuilder("AA");
            for (int i = 0; i < serineCount; i++)
            {
                sequence.Append("SA");
            }
            sequence.Append('K');

            var protein = new Protein(sequence.ToString(), "measure");
            return protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
        }

        /// <summary>
        /// Route count as a function of candidate-site count and glycan-box size. Route enumeration legality
        /// is combinatorial, so the spectrum sets the costs but not how many routes exist; a real spectrum is
        /// still used so the cost path is exercised exactly as it is in a search.
        /// </summary>
        [Test]
        [Explicit("Measurement harness, not a correctness test. Run deliberately.")]
        public static void Measure_RouteCount_ByCandidateSitesAndBoxSize()
        {
            var commonParameters = new CommonParameters(dissociationType: DissociationType.ETD, trimMsMsPeaks: false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"GlycoTestData\181217_Fusion_(LC2)_NewObj_Serum_deSA_Jacalin_HRM_4h_ETD_HCD_DDA_mz(400_1200)_21707.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scan = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).First();

            var rows = new List<string> { "SiteCount,BoxSize,Routes,HighestScorePaths,ElapsedMs" };
            TestContext.WriteLine("SiteCount  BoxSize  Routes  HighestScorePaths  ElapsedMs");

            foreach (int siteCount in new[] { 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 25, 30, 40 })
            {
                var peptide = MakePeptideWithSites(siteCount);
                var products = new List<Product>();
                peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

                var modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
                Assert.That(modPos.Count, Is.EqualTo(siteCount), "Peptide construction did not yield the intended site count.");

                foreach (int boxSize in new[] { 1, 2, 3 })
                {
                    if (boxSize > siteCount)
                    {
                        continue;
                    }

                    var glycanBox = OGlycanBoxes.First(p => p.NumberOfMods == boxSize);
                    var childBoxes = GlycanBox.BuildChildOGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds).ToArray();

                    var localizationGraph = new LocalizationGraph(modPos, glycanBox, childBoxes, -1);
                    LocalizationGraph.LocalizeOGlycan(localizationGraph, scan, commonParameters.ProductMassTolerance, products);

                    var stopwatch = System.Diagnostics.Stopwatch.StartNew();
                    var routes = LocalizationGraph.GetAllPaths_CalP(localizationGraph, 0.1, products.Count);
                    stopwatch.Stop();

                    var highestScorePaths = LocalizationGraph.GetAllHighestScorePaths(localizationGraph.array, localizationGraph.ChildModBoxes);

                    TestContext.WriteLine($"{siteCount,9}  {boxSize,7}  {routes.Count,6}  {highestScorePaths.Count,17}  {stopwatch.ElapsedMilliseconds,9}");
                    rows.Add($"{siteCount},{boxSize},{routes.Count},{highestScorePaths.Count},{stopwatch.ElapsedMilliseconds}");
                }
            }

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "route_count_measurement.csv");
            File.WriteAllLines(outputPath, rows);
            TestContext.WriteLine("Wrote " + outputPath);
        }

        /// <summary>
        /// Phase B, first increment. Two questions on a real peptide/spectrum pair:
        /// (1) do decoy sites actually compete inside the graph, or are they inert;
        /// (2) B2a -- does adding decoy positions shift the probabilities of the REAL sites? GetLocalFragment
        /// brackets ions by ADJACENT ModPos entries, so a decoy inserted between two real sites re-partitions
        /// the fragment evidence. If that perturbation is large the decoy method measures something it is
        /// itself changing, and paired graphs become necessary.
        /// </summary>
        [Test]
        [Explicit("Measurement harness, not a correctness test. Run deliberately.")]
        public static void Measure_DecoySiteCompetition_AndTargetPerturbation()
        {
            var commonParameters = new CommonParameters(dissociationType: DissociationType.ETD, trimMsMsPeaks: false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"GlycoTestData\181217_Fusion_(LC2)_NewObj_Serum_deSA_Jacalin_HRM_4h_ETD_HCD_DDA_mz(400_1200)_21707.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scan = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).First();

            // A peptide/spectrum pair with genuine matching fragments, so costs are real rather than all-zero.
            var protein = new Protein("AATVGSLAGQPLQER", "P16150");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            var products = new List<Product>();
            peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

            var glycanBox = OGlycanBoxes[1];
            var childBoxes = GlycanBox.BuildChildOGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds).ToArray();

            Dictionary<int, double> RunGraph(SortedDictionary<int, string> modPos)
            {
                var graph = new LocalizationGraph(modPos, glycanBox, childBoxes, -1);
                LocalizationGraph.LocalizeOGlycan(graph, scan, commonParameters.ProductMassTolerance, products);
                var routes = LocalizationGraph.GetAllPaths_CalP(graph, 0.1, products.Count);
                var pairs = routes.SelectMany(p => p.ModSitePairs).Distinct().ToList();
                LocalizationGraph.CalProbabilityForModSitePair(routes, pairs);
                return pairs.GroupBy(p => p.SiteIndex).ToDictionary(g => g.Key, g => g.Sum(p => p.Probability));
            }

            // Arm 1: targets only, the current behaviour.
            var targetsOnly = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
            var targetsOnlyProbabilities = RunGraph(targetsOnly);

            // Arm 2: targets plus alanine decoys. The decoy label must match the target motif of the glycans
            // in this box, because MotifCheck admits a glycan onto a position by comparing the glycan's Target
            // to the recorded motif. O-Pair loads each composition twice, once targeting S and once T, so a
            // fixed label would only ever compete against half the glycan instances.
            string boxTargetMotif = GlycanBox.GlobalOGlycans[glycanBox.ModIds[0]].Target.ToString();
            TestContext.WriteLine($"Box target motif: {boxTargetMotif}");
            var withDecoys = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
            var decoySites = GlycoSpectralMatch.AddDecoyModSites(withDecoys, peptide, new string[] { "A" }, boxTargetMotif);
            var withDecoysProbabilities = RunGraph(withDecoys);

            TestContext.WriteLine($"Peptide {peptide.BaseSequence}");
            TestContext.WriteLine($"Real sites: {string.Join(",", targetsOnly.Select(p => p.Key + p.Value))}");
            TestContext.WriteLine($"Decoy sites: {string.Join(",", decoySites.OrderBy(p => p))}");
            TestContext.WriteLine("");
            TestContext.WriteLine("Site  Kind    TargetsOnly  WithDecoys   Delta");

            foreach (var site in withDecoysProbabilities.Keys.OrderBy(p => p))
            {
                bool isDecoy = decoySites.Contains(site);
                targetsOnlyProbabilities.TryGetValue(site, out double before);
                double after = withDecoysProbabilities[site];
                string kind = isDecoy ? "DECOY " : "target";
                string beforeText = isDecoy ? "     -" : before.ToString("F4");
                TestContext.WriteLine($"{site,4}  {kind}  {beforeText,11}  {after,10:F4}  {(isDecoy ? "" : (after - before).ToString("+0.0000;-0.0000"))}");
            }

            double decoyMass = withDecoysProbabilities.Where(p => decoySites.Contains(p.Key)).Sum(p => p.Value);
            double maxTargetShift = targetsOnlyProbabilities.Keys
                .Where(withDecoysProbabilities.ContainsKey)
                .Select(k => Math.Abs(withDecoysProbabilities[k] - targetsOnlyProbabilities[k]))
                .DefaultIfEmpty(0).Max();

            TestContext.WriteLine("");
            TestContext.WriteLine($"Total probability landing on decoy sites : {decoyMass:F4}");
            TestContext.WriteLine($"Largest shift in any real site           : {maxTargetShift:F4}");
        }

        /// <summary>
        /// Separates the two readings of M4. A decoy can take probability for two different reasons:
        /// because the fragment evidence genuinely cannot separate it from a real site (signal), or because
        /// inserting it between two ModPos entries re-partitions the fragment windows in GetLocalFragment
        /// (artifact). Those were confounded in M4 because the decoy that took mass happened to be adjacent
        /// to the real site.
        /// <para>
        /// This sweeps a single decoy across every non-candidate position of one peptide, one at a time, and
        /// reports the decoy's probability against the number of matched backbone ions that actually separate
        /// it from the real site. If probability tracks separating evidence, the effect is signal; if it
        /// tracks mere insertion, it is artifact.
        /// </para>
        /// </summary>
        [Test]
        [Explicit("Measurement harness, not a correctness test. Run deliberately.")]
        public static void Measure_DecoyProbability_VersusSeparatingEvidence()
        {
            var commonParameters = new CommonParameters(dissociationType: DissociationType.ETD, trimMsMsPeaks: false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"GlycoTestData\181217_Fusion_(LC2)_NewObj_Serum_deSA_Jacalin_HRM_4h_ETD_HCD_DDA_mz(400_1200)_21707.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scan = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).First();

            var protein = new Protein("AATVGSLAGQPLQER", "P16150");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            var products = new List<Product>();
            peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

            var glycanBox = OGlycanBoxes[1];
            var childBoxes = GlycanBox.BuildChildOGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds).ToArray();
            string boxTargetMotif = GlycanBox.GlobalOGlycans[glycanBox.ModIds[0]].Target.ToString();

            var targetsOnly = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
            // The site this box can actually reach, i.e. the one carrying the box's target motif.
            int realSite = targetsOnly.First(p => p.Value == boxTargetMotif).Key;

            // Matched backbone ions of the unmodified peptide. Used only to ask which backbone cleavages are
            // observed at all; it is a proxy, since a fragment spanning the glycosite carries the glycan mass.
            var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, products, commonParameters);
            var matchedResiduePositions = matchedIons
                .Select(p => p.NeutralTheoreticalProduct.ResiduePosition)
                .Distinct().OrderBy(p => p).ToList();
            TestContext.WriteLine($"Peptide {peptide.BaseSequence}, box target {boxTargetMotif}, real site {realSite}");
            TestContext.WriteLine($"Matched ion residue positions: {string.Join(",", matchedResiduePositions)}");
            TestContext.WriteLine("");
            TestContext.WriteLine("DecoyPos Residue Dist SeparatingIons  DecoyProb  RealSiteProb");

            var rows = new List<string> { "DecoyPos,Residue,Distance,SeparatingIons,DecoyProb,RealSiteProb" };

            for (int r = 0; r < peptide.BaseSequence.Length; r++)
            {
                int siteKey = r + 2;
                if (targetsOnly.ContainsKey(siteKey))
                {
                    continue; // real candidate site, not a decoy position
                }

                var modPos = new SortedDictionary<int, string>(targetsOnly.ToDictionary(p => p.Key, p => p.Value));
                modPos[siteKey] = boxTargetMotif;

                var graph = new LocalizationGraph(modPos, glycanBox, childBoxes, -1);
                LocalizationGraph.LocalizeOGlycan(graph, scan, commonParameters.ProductMassTolerance, products);
                var routes = LocalizationGraph.GetAllPaths_CalP(graph, 0.1, products.Count);
                var pairs = routes.SelectMany(p => p.ModSitePairs).Distinct().ToList();
                LocalizationGraph.CalProbabilityForModSitePair(routes, pairs);
                var bySite = pairs.GroupBy(p => p.SiteIndex).ToDictionary(g => g.Key, g => g.Sum(p => p.Probability));

                bySite.TryGetValue(siteKey, out double decoyProbability);
                bySite.TryGetValue(realSite, out double realProbability);

                // A matched ion separates the decoy from the real site when its cleavage falls between them.
                int low = Math.Min(siteKey, realSite);
                int high = Math.Max(siteKey, realSite);
                int separatingIons = matchedResiduePositions.Count(p => p >= low - 1 && p < high - 1);

                TestContext.WriteLine($"{siteKey,8} {peptide.BaseSequence[r],7} {Math.Abs(siteKey - realSite),4} {separatingIons,14}  {decoyProbability,9:F4}  {realProbability,12:F4}");
                rows.Add($"{siteKey},{peptide.BaseSequence[r]},{Math.Abs(siteKey - realSite)},{separatingIons},{decoyProbability:F4},{realProbability:F4}");
            }

            File.WriteAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "decoy_bracketing_sweep.csv"), rows);
        }

        /// <summary>
        /// Replicates the M5 bracketing sweep across every glycopeptide a real search identifies, rather than
        /// the single hand-built peptide M5 used. Runs a glyco search, reads the identifications back out of
        /// the written psmtsv, and for each one sweeps a single decoy across every non-candidate position.
        /// <para>
        /// Aggregating decoy probability against separating-ion count over many peptides is what turns M5
        /// from a striking case into a result.
        /// </para>
        /// </summary>
        [Test]
        [Explicit("Measurement harness, not a correctness test. Run deliberately.")]
        public static void Measure_BracketingSweep_AcrossIdentifiedGlycopeptides()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TESTGlycoSweep");
            if (Directory.Exists(outputFolder))
            {
                Directory.Delete(outputFolder, true);
            }
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Nett.Toml.ReadFile<GlycoSearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSnip.toml"),
                MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.WriteContaminants = false;

            var db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoProteinFASTA_7proteins.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoPepMix_snip.mzML");

            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) },
                new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder).Run();

            string psmtsv = Directory.GetFiles(outputFolder, "*.psmtsv", SearchOption.AllDirectories)
                .OrderByDescending(p => new FileInfo(p).Length).First();
            var lines = File.ReadAllLines(psmtsv);
            var header = lines[0].Split('\t');
            int iBaseSeq = Array.IndexOf(header, "Base Sequence");
            int iScan = Array.IndexOf(header, "Scan Number");
            int iGlycanMass = Array.IndexOf(header, "GlycanMass");
            Assert.That(Math.Min(iBaseSeq, Math.Min(iScan, iGlycanMass)), Is.GreaterThanOrEqualTo(0),
                "Required psmtsv columns not found.");

            // Load the spectra once and index by scan number.
            var commonParameters = new CommonParameters(dissociationType: DissociationType.EThcD, trimMsMsPeaks: false);
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scansByNumber = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters)
                .GroupBy(p => p.OneBasedScanNumber).ToDictionary(g => g.Key, g => g.First());

            var rows = new List<string> { "Peptide,Scan,DecoyPos,Residue,Distance,SeparatingIons,DecoyProb,RealSiteProb" };
            int peptidesSwept = 0;

            foreach (var line in lines.Skip(1).Where(l => !string.IsNullOrWhiteSpace(l)))
            {
                var fields = line.Split('\t');
                if (fields.Length <= Math.Max(iBaseSeq, Math.Max(iScan, iGlycanMass))) continue;

                string baseSequence = fields[iBaseSeq].Trim();
                if (string.IsNullOrWhiteSpace(baseSequence) || baseSequence.Contains("|")) continue; // ambiguous
                if (!int.TryParse(fields[iScan].Trim(), out int scanNumber)) continue;
                if (!double.TryParse(fields[iGlycanMass].Trim(), out double glycanMass)) continue;
                if (!scansByNumber.TryGetValue(scanNumber, out var scan)) continue;

                var peptide = new Protein(baseSequence, "sweep")
                    .Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>()).FirstOrDefault();
                if (peptide == null) continue;

                // The box whose mass matches what the search assigned.
                var glycanBox = OGlycanBoxes.FirstOrDefault(p => Math.Abs(p.Mass - glycanMass) < 0.02);
                if (glycanBox == null) continue;

                var products = new List<Product>();
                peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
                var childBoxes = GlycanBox.BuildChildOGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds).ToArray();
                string boxTargetMotif = GlycanBox.GlobalOGlycans[glycanBox.ModIds[0]].Target.ToString();

                var targetsOnly = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
                // Only peptides with a genuine localization choice are informative here.
                if (targetsOnly.Count(p => p.Value == boxTargetMotif) < 2) continue;
                int realSite = targetsOnly.First(p => p.Value == boxTargetMotif).Key;

                // Localization is scored on c/zDot ions, so it must run against the electron-based CHILD scan.
                // This data is HCD-pd-EThcD, so GetMs2Scans hands back the HCD parent, which carries no c/z at
                // all -- passing it produces a graph score of exactly zero and a uniform posterior that looks
                // like "decoys always win". Same trap as the localization-scan selection bug in #2692.
                var localizationScan = scan.ChildScans.FirstOrDefault(c =>
                    c.TheScan.DissociationType.HasValue
                    && c.TheScan.DissociationType != DissociationType.Autodetect
                    && GlycoPeptides.DissociationTypeContainETD(c.TheScan.DissociationType.Value, commonParameters.CustomIons))
                    ?? scan;

                var matchedPositions = MetaMorpheusEngine.MatchFragmentIons(localizationScan, products, commonParameters)
                    .Select(p => p.NeutralTheoreticalProduct.ResiduePosition).Distinct().OrderBy(p => p).ToList();

                peptidesSwept++;

                // Diagnostic: if the peptide/scan/box triple is not a genuine match, every cost is zero and
                // the posterior collapses to uniform 1/n, which would look like "decoys always take mass".
                {
                    var baseGraph = new LocalizationGraph(targetsOnly, glycanBox, childBoxes, -1);
                    LocalizationGraph.LocalizeOGlycan(baseGraph, localizationScan, commonParameters.ProductMassTolerance, products);
                    var baseRoutes = LocalizationGraph.GetAllPaths_CalP(baseGraph, 0.1, products.Count);
                    var basePairs = baseRoutes.SelectMany(p => p.ModSitePairs).Distinct().ToList();
                    LocalizationGraph.CalProbabilityForModSitePair(baseRoutes, basePairs);
                    string dist = string.Join(" ", basePairs.GroupBy(p => p.SiteIndex).OrderBy(g => g.Key)
                        .Select(g => $"{g.Key}:{g.Sum(p => p.Probability):F3}"));
                    TestContext.WriteLine($"DIAG {baseSequence} scan={scanNumber} sites={targetsOnly.Count} matchedIons={matchedPositions.Count} graphTotalScore={baseGraph.TotalScore:F3} targetsOnly[{dist}]");
                }

                for (int r = 0; r < peptide.BaseSequence.Length; r++)
                {
                    int siteKey = r + 2;
                    if (targetsOnly.ContainsKey(siteKey)) continue;

                    var modPos = new SortedDictionary<int, string>(targetsOnly.ToDictionary(p => p.Key, p => p.Value));
                    modPos[siteKey] = boxTargetMotif;

                    var graph = new LocalizationGraph(modPos, glycanBox, childBoxes, -1);
                    LocalizationGraph.LocalizeOGlycan(graph, localizationScan, commonParameters.ProductMassTolerance, products);
                    var routes = LocalizationGraph.GetAllPaths_CalP(graph, 0.1, products.Count);
                    if (routes.Count == 0) continue;
                    var pairs = routes.SelectMany(p => p.ModSitePairs).Distinct().ToList();
                    LocalizationGraph.CalProbabilityForModSitePair(routes, pairs);
                    var bySite = pairs.GroupBy(p => p.SiteIndex).ToDictionary(g => g.Key, g => g.Sum(p => p.Probability));

                    bySite.TryGetValue(siteKey, out double decoyProbability);
                    bySite.TryGetValue(realSite, out double realProbability);

                    int low = Math.Min(siteKey, realSite);
                    int high = Math.Max(siteKey, realSite);
                    int separating = matchedPositions.Count(p => p >= low - 1 && p < high - 1);

                    rows.Add($"{baseSequence},{scanNumber},{siteKey},{peptide.BaseSequence[r]},{Math.Abs(siteKey - realSite)},{separating},{decoyProbability:F4},{realProbability:F4}");
                }
            }

            string csv = Path.Combine(TestContext.CurrentContext.TestDirectory, "bracketing_sweep_all.csv");
            File.WriteAllLines(csv, rows);
            TestContext.WriteLine($"Peptides swept: {peptidesSwept}; decoy arms: {rows.Count - 1}");
            TestContext.WriteLine("Wrote " + csv);

            // Aggregate: does decoy probability track separating evidence across all peptides?
            var parsed = rows.Skip(1).Select(r => r.Split(',')).ToList();
            TestContext.WriteLine("");
            TestContext.WriteLine("SeparatingIons  Arms  MeanDecoyProb  MaxDecoyProb  ArmsWithProb>0.01");
            foreach (var group in parsed.GroupBy(f => int.Parse(f[5])).OrderBy(g => g.Key))
            {
                var probabilities = group.Select(f => double.Parse(f[6])).ToList();
                TestContext.WriteLine($"{group.Key,14}  {probabilities.Count,4}  {probabilities.Average(),13:F4}  {probabilities.Max(),12:F4}  {probabilities.Count(p => p > 0.01),17}");
            }

            Directory.Delete(outputFolder, true);
        }

        /// <summary>
        /// How much of a real glycopeptide corpus each shipped O-glycan database can actually represent.
        /// <para>
        /// M8 dropped 45 of 114 multi-site identifications because Byonic's glycan composition had no
        /// mass-matching box in the default database. That is a coverage bias, not a random sample, and it
        /// would propagate straight into an FLR: whichever glycans the database cannot express are silently
        /// excluded from the error estimate.
        /// </para>
        /// </summary>
        [Test]
        [Explicit("Measurement harness. Requires the PXD017646 glycoPSM table locally.")]
        public static void Measure_GlycanDatabaseCoverage_AgainstRealCorpus()
        {
            string tablePath = Path.Combine(@"E:\CodeReview\localization\data\raw", "StcEmix_35trig_EThcD25_rep1_GlycoPSMs.txt");
            if (!File.Exists(tablePath))
            {
                Assert.Ignore("PXD017646 glycoPSM table not present locally.");
            }

            var lines = File.ReadAllLines(tablePath);
            var header = lines[0].Split('\t');
            int iFrag = Array.IndexOf(header, "Fragmentation");
            int iCalcMH = Array.IndexOf(header, "calcMH");
            int iPepMass = Array.IndexOf(header, "PepMassNoMod");
            int iSites = Array.IndexOf(header, "#GlycoSitesOnPep");
            int iGlycans = Array.IndexOf(header, "Glycans");
            const double protonMass = 1.00727646677;

            // The multi-site ETD identifications an FLR would actually be computed over.
            var corpus = new List<(double Mass, string Composition)>();
            foreach (var line in lines.Skip(1).Where(l => !string.IsNullOrWhiteSpace(l)))
            {
                var f = line.Split('\t');
                if (f.Length <= Math.Max(iSites, iGlycans)) continue;
                if (!f[iFrag].Contains("ETD", StringComparison.OrdinalIgnoreCase)) continue;
                if (!int.TryParse(f[iSites].Trim(), out int sites) || sites < 2) continue;
                if (!double.TryParse(f[iCalcMH].Trim(), out double calcMH)) continue;
                if (!double.TryParse(f[iPepMass].Trim(), out double pepMass)) continue;
                corpus.Add((calcMH - protonMass - pepMass, f[iGlycans].Trim()));
            }
            TestContext.WriteLine($"Multi-site ETD identifications in corpus: {corpus.Count}");
            TestContext.WriteLine("");

            var originalGlycans = GlycanBox.GlobalOGlycans;
            var originalBoxes = GlycanBox.OGlycanBoxes;
            var databaseNames = new[]
            {
                "OGlycan.gdb",
                "OGlycan_withIsobaric.gdb",
                "Olgycan Database 28 glycans.txt",
                "Olgycan Database 32 glycans.txt",
                "Olgycan Database 36 glycans with Mann.txt",
            };

            try
            {
                TestContext.WriteLine("Database                                    Glycans  Boxes   Covered  Coverage");
                foreach (var name in databaseNames)
                {
                    string path = GlobalVariables.OGlycanDatabasePaths.FirstOrDefault(p => Path.GetFileName(p) == name);
                    if (path == null)
                    {
                        TestContext.WriteLine($"{name,-42}  (not registered)");
                        continue;
                    }

                    GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(path, true, true).ToArray();
                    var boxes = GlycanBox.BuildOGlycanBoxes(3).OrderBy(p => p.Mass).ToArray();
                    var boxMasses = boxes.Select(p => p.Mass).OrderBy(m => m).ToArray();

                    int covered = corpus.Count(entry => boxMasses.Any(m => Math.Abs(m - entry.Mass) < 0.02));
                    TestContext.WriteLine($"{name,-42} {GlycanBox.GlobalOGlycans.Length,8} {boxes.Length,6} {covered,9} {covered / (double)corpus.Count,9:P1}");
                }

                // Is the limit the database size, or the number of glycans allowed per peptide?
                TestContext.WriteLine("");
                TestContext.WriteLine("Database                                    maxOGlycanNum  Boxes   Covered  Coverage");
                foreach (var name in new[] { "OGlycan.gdb", "Olgycan Database 36 glycans with Mann.txt" })
                {
                    string path = GlobalVariables.OGlycanDatabasePaths.FirstOrDefault(p => Path.GetFileName(p) == name);
                    if (path == null) continue;
                    GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(path, true, true).ToArray();

                    foreach (int maxNum in new[] { 2, 3, 4 })
                    {
                        var boxes = GlycanBox.BuildOGlycanBoxes(maxNum).ToArray();
                        var masses = boxes.Select(p => p.Mass).ToArray();
                        int covered = corpus.Count(e => masses.Any(m => Math.Abs(m - e.Mass) < 0.02));
                        TestContext.WriteLine($"{name,-42} {maxNum,14} {boxes.Length,6} {covered,9} {covered / (double)corpus.Count,9:P1}");
                    }
                }

                // What the default database misses, by composition.
                string defaultPath = GlobalVariables.OGlycanDatabasePaths.First(p => Path.GetFileName(p) == "OGlycan.gdb");
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(defaultPath, true, true).ToArray();
                var defaultMasses = GlycanBox.BuildOGlycanBoxes(3).Select(p => p.Mass).ToArray();
                var missing = corpus.Where(e => !defaultMasses.Any(m => Math.Abs(m - e.Mass) < 0.02))
                    .GroupBy(e => e.Composition).OrderByDescending(g => g.Count()).ToList();

                TestContext.WriteLine("");
                TestContext.WriteLine($"Compositions the default database cannot express ({missing.Sum(g => g.Count())} identifications):");
                foreach (var group in missing.Take(15))
                {
                    TestContext.WriteLine($"  {group.Count(),4}  {group.Key}");
                }
            }
            finally
            {
                GlycanBox.GlobalOGlycans = originalGlycans;
                GlycanBox.OGlycanBoxes = originalBoxes;
            }
        }

        /// <summary>
        /// Phase C: a real false localization rate, rather than the per-arm rates of M8-M12.
        /// <para>
        /// All decoy positions compete at once, as they would in a search, and the question asked per
        /// identification is the one an FLR actually asks: <b>does the top-scoring placement land on a site
        /// that cannot carry the glycan?</b> Each such win is evidence that some proportion of target-site
        /// wins are equally unfounded, scaled by how many target sites there were to hit relative to decoy
        /// sites.
        /// </para>
        /// <para>
        /// Estimator, following the decoy-amino-acid method: a decoy win on an identification offering
        /// T target sites and D decoy sites implies T/D expected false target wins, so
        /// <c>FLR = sum(T_i/D_i over decoy wins) / (target wins)</c>. Computed over the multi-site
        /// population only, because single-candidate-site identifications localize by construction and
        /// including them dilutes the rate (M6: they are about half the corpus).
        /// </para>
        /// </summary>
        [Test]
        [Explicit("Measurement harness. Requires the PXD017646 raw files locally.")]
        public static void Measure_FalseLocalizationRate_AllDecoysAtOnce()
        {
            string dataDirectory = @"E:\CodeReview\localization\data\raw";
            var runs = Directory.GetFiles(dataDirectory, "*StcEmix_35trig_*.raw")
                .Select(raw =>
                {
                    string stem = Path.GetFileNameWithoutExtension(raw);
                    int idx = stem.IndexOf("StcEmix_35trig_", StringComparison.Ordinal);
                    return (Raw: raw, Table: Path.Combine(dataDirectory, stem.Substring(idx) + "_GlycoPSMs.txt"));
                })
                .Where(p => File.Exists(p.Table)).OrderBy(p => p.Raw).ToList();
            if (runs.Count == 0) Assert.Ignore("No PXD017646 raw/glycoPSM pairs present locally.");

            var commonParameters = new CommonParameters(dissociationType: DissociationType.EThcD, trimMsMsPeaks: false);
            const double protonMass = 1.00727646677;
            // Ala first, then Leu and Gly, following the decoy-amino-acid method. None can carry an
            // O-glycan. Using all three roughly triples the decoy-site count, which both stabilises the
            // estimate and stops peptides being dropped merely for lacking an alanine.
            string[] decoyResidues = { "A" };
            var rows = new List<string> { "Run,Activation,Peptide,Scan,TargetSites,DecoySites,WinnerIsDecoy,WinnerProb,TargetRatio,WinnerResidue,BestTargetProbNoDecoys" };

            foreach (var run in runs)
            {
                string runName = Path.GetFileNameWithoutExtension(run.Raw).Replace("2019_09_16_StcEmix_35trig_", "");
                string activation = System.Text.RegularExpressions.Regex.Replace(runName, "_rep\\d", "");

                var file = new MyFileManager(true).LoadFile(run.Raw, commonParameters);
                var scansByNumber = MetaMorpheusTask.GetMs2Scans(file, run.Raw, commonParameters)
                    .GroupBy(p => p.OneBasedScanNumber).ToDictionary(g => g.Key, g => g.First());

                var lines = File.ReadAllLines(run.Table);
                var header = lines[0].Split('\t');
                int iSeq = Array.IndexOf(header, "Sequence"), iScan = Array.IndexOf(header, "ScanNumber");
                int iFrag = Array.IndexOf(header, "Fragmentation"), iCalcMH = Array.IndexOf(header, "calcMH");
                int iPepMass = Array.IndexOf(header, "PepMassNoMod"), iSites = Array.IndexOf(header, "#GlycoSitesOnPep");

                foreach (var line in lines.Skip(1).Where(l => !string.IsNullOrWhiteSpace(l)))
                {
                    var f = line.Split('\t');
                    if (f.Length <= iSites) continue;
                    if (!f[iFrag].Contains("ETD", StringComparison.OrdinalIgnoreCase)) continue;
                    if (!int.TryParse(f[iSites].Trim(), out int byonicSites) || byonicSites < 2) continue;
                    if (!int.TryParse(f[iScan].Trim(), out int scanNumber)) continue;
                    if (!scansByNumber.TryGetValue(scanNumber, out var scan)) continue;
                    if (!double.TryParse(f[iCalcMH].Trim(), out double calcMH)) continue;
                    if (!double.TryParse(f[iPepMass].Trim(), out double pepMassNoMod)) continue;

                    string sequence = f[iSeq].Trim();
                    if (string.IsNullOrWhiteSpace(sequence) || !sequence.All(char.IsLetter)) continue;

                    var glycanBox = OGlycanBoxes.FirstOrDefault(p => Math.Abs(p.Mass - (calcMH - protonMass - pepMassNoMod)) < 0.02);
                    if (glycanBox == null) continue;

                    var peptide = new Protein(sequence, "byonic")
                        .Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>()).FirstOrDefault();
                    if (peptide == null) continue;

                    var products = new List<Product>();
                    peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
                    var childBoxes = GlycanBox.BuildChildOGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds).ToArray();
                    string boxTargetMotif = GlycanBox.GlobalOGlycans[glycanBox.ModIds[0]].Target.ToString();

                    var modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
                    int targetSites = modPos.Count(p => p.Value == boxTargetMotif);
                    if (targetSites < 2) continue;

                    // Every decoy residue competes, all at once.
                    var decoySites = GlycoSpectralMatch.AddDecoyModSites(modPos, peptide, decoyResidues, boxTargetMotif);
                    if (decoySites.Count == 0) continue;
                    if (!GraphCheck(modPos, glycanBox)) continue;

                    // Paired control: the same identification scored WITHOUT decoys. If adding decoy
                    // positions is merely revealing ambiguity, the best target site should keep roughly the
                    // probability it had here. If it collapses as decoy density rises, then inserting
                    // positions is re-partitioning the fragment windows GetLocalFragment brackets by, and the
                    // decoys are destroying the evidence rather than competing for it.
                    var targetsOnlyPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
                    double bestTargetProbNoDecoys = double.NaN;
                    if (GraphCheck(targetsOnlyPos, glycanBox))
                    {
                        var baseGraph = new LocalizationGraph(targetsOnlyPos, glycanBox, childBoxes, -1);
                        LocalizationGraph.LocalizeOGlycan(baseGraph, scan, commonParameters.ProductMassTolerance, products);
                        var baseRoutes = LocalizationGraph.GetAllPaths_CalP(baseGraph, 0.1, products.Count);
                        if (baseRoutes.Count > 0)
                        {
                            var basePairs = baseRoutes.SelectMany(p => p.ModSitePairs).Distinct().ToList();
                            LocalizationGraph.CalProbabilityForModSitePair(baseRoutes, basePairs);
                            bestTargetProbNoDecoys = basePairs.GroupBy(p => p.SiteIndex)
                                .Select(g => g.Sum(p => p.Probability)).DefaultIfEmpty(0).Max();
                        }
                    }

                    var graph = new LocalizationGraph(modPos, glycanBox, childBoxes, -1);
                    LocalizationGraph.LocalizeOGlycan(graph, scan, commonParameters.ProductMassTolerance, products);
                    if (graph.TotalScore < 2.0) continue;

                    var routes = LocalizationGraph.GetAllPaths_CalP(graph, 0.1, products.Count);
                    if (routes.Count == 0) continue;
                    var pairs = routes.SelectMany(p => p.ModSitePairs).Distinct().ToList();
                    LocalizationGraph.CalProbabilityForModSitePair(routes, pairs);

                    var bySite = pairs.GroupBy(p => p.SiteIndex)
                        .Select(g => new { Site = g.Key, Probability = g.Sum(p => p.Probability) })
                        .OrderByDescending(x => x.Probability).ToList();
                    if (bySite.Count == 0) continue;

                    var winner = bySite[0];
                    bool winnerIsDecoy = decoySites.Contains(winner.Site);
                    double ratio = targetSites / (double)decoySites.Count;
                    // Site keys are r+2, so the residue is at index key-2.
                    char winnerResidue = peptide.BaseSequence[winner.Site - 2];

                    rows.Add($"{runName},{activation},{sequence},{scanNumber},{targetSites},{decoySites.Count},{winnerIsDecoy},{winner.Probability:F4},{ratio:F4},{winnerResidue},{bestTargetProbNoDecoys:F4}");
                }
            }

            string csv = Path.Combine(TestContext.CurrentContext.TestDirectory, "flr_all_decoys.csv");
            File.WriteAllLines(csv, rows);
            TestContext.WriteLine($"Identifications scored with all decoys competing: {rows.Count - 1}");
            TestContext.WriteLine("Wrote " + csv);

            var parsed = rows.Skip(1).Select(r => r.Split(',')).ToList();
            if (parsed.Count == 0) return;

            TestContext.WriteLine("");
            TestContext.WriteLine("Activation      n   DecoyWins   RawDecoyRate   RatioNormalisedFLR");
            foreach (var group in parsed.GroupBy(x => x[1]).OrderBy(g => g.Key))
            {
                var decoyWins = group.Where(x => bool.Parse(x[6])).ToList();
                int targetWins = group.Count() - decoyWins.Count;
                double expectedFalseTargets = decoyWins.Sum(x => double.Parse(x[8]));
                double flr = targetWins > 0 ? expectedFalseTargets / targetWins : double.NaN;
                TestContext.WriteLine($"{group.Key,-12} {group.Count(),4}   {decoyWins.Count,9}   {decoyWins.Count / (double)group.Count(),12:P1}   {flr,18:P1}");
            }
        }

        /// <summary>
        /// Scans every plausible decoy residue, one at a time, and reports how often each wins per site it
        /// offers -- alongside how close its sites sit to genuine glycosites.
        /// <para>
        /// M14 showed alanine, leucine and glycine behave very differently as decoys (3.0%, 10.4%, 15.5% win
        /// rate per site) and attributed the difference to positional correlation with real sites rather than
        /// to chemistry. This tests that directly across the whole residue alphabet: if win rate tracks mean
        /// distance to the nearest real site, the requirement on a decoy residue is positional, and the right
        /// choice is measurable per dataset rather than inherited from the phospho literature.
        /// </para>
        /// </summary>
        [Test]
        [Explicit("Measurement harness. Requires the PXD017646 raw files locally.")]
        public static void Measure_DecoyResidueScan()
        {
            string dataDirectory = @"E:\CodeReview\localization\data\raw";
            var runs = Directory.GetFiles(dataDirectory, "*35trig_*.raw")
                .Select(raw =>
                {
                    // Raw stems are <date>_<sample>_35trig_<method>[_repN]; tables drop the date.
                    string stem = Path.GetFileNameWithoutExtension(raw);
                    var m = System.Text.RegularExpressions.Regex.Match(stem, @"^\d{4}_\d{2}_\d{2}_(.+)$");
                    string key = m.Success ? m.Groups[1].Value : stem;
                    return (Raw: raw, Table: Path.Combine(dataDirectory, key + "_GlycoPSMs.txt"), Key: key);
                })
                .Where(p => File.Exists(p.Table)).OrderBy(p => p.Key).ToList();
            if (runs.Count == 0) Assert.Ignore("No raw/glycoPSM pairs present locally.");
            TestContext.WriteLine($"Runs: {runs.Count}");

            // Everything that cannot carry an O-glycan. S and T are targets; Y is excluded because rare
            // Y-linked O-glycans exist and it would not be a clean null.
            string[] candidateResidues = { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "V", "W" };

            var commonParameters = new CommonParameters(dissociationType: DissociationType.EThcD, trimMsMsPeaks: false);
            const double protonMass = 1.00727646677;

            // "<sample>|<residue>" -> counts. Sample is tracked so the ranking can be tested for
            // transferability: a decoy residue that behaves well only on mucins is not a general choice.
            var sitesOffered = new Dictionary<string, int>();
            var wins = new Dictionary<string, int>();
            var distanceSum = new Dictionary<string, double>();
            var identifications = new Dictionary<string, int>();
            void Bump(Dictionary<string, int> d, string k, int by) { d.TryGetValue(k, out int v); d[k] = v + by; }

            foreach (var run in runs)
            {
                Dictionary<int, Ms2ScanWithSpecificMass> scansByNumber;
                try
                {
                    // A raw file still being downloaded is unreadable; skip it rather than abort the scan.
                    var file = new MyFileManager(true).LoadFile(run.Raw, commonParameters);
                    scansByNumber = MetaMorpheusTask.GetMs2Scans(file, run.Raw, commonParameters)
                        .GroupBy(p => p.OneBasedScanNumber).ToDictionary(g => g.Key, g => g.First());
                }
                catch (Exception e)
                {
                    TestContext.WriteLine($"SKIPPED {Path.GetFileName(run.Raw)}: {e.Message}");
                    continue;
                }

                var lines = File.ReadAllLines(run.Table);
                var header = lines[0].Split('\t');
                int iSeq = Array.IndexOf(header, "Sequence"), iScan = Array.IndexOf(header, "ScanNumber");
                int iFrag = Array.IndexOf(header, "Fragmentation"), iCalcMH = Array.IndexOf(header, "calcMH");
                int iPepMass = Array.IndexOf(header, "PepMassNoMod"), iSites = Array.IndexOf(header, "#GlycoSitesOnPep");
                if (Math.Min(iSeq, Math.Min(iScan, Math.Min(iFrag, Math.Min(iCalcMH, Math.Min(iPepMass, iSites))))) < 0) continue;

                foreach (var line in lines.Skip(1).Where(l => !string.IsNullOrWhiteSpace(l)))
                {
                    var f = line.Split('\t');
                    if (f.Length <= iSites) continue;
                    if (!f[iFrag].Contains("ETD", StringComparison.OrdinalIgnoreCase)) continue;
                    if (!int.TryParse(f[iSites].Trim(), out int bs) || bs < 2) continue;
                    if (!int.TryParse(f[iScan].Trim(), out int scanNumber)) continue;
                    if (!scansByNumber.TryGetValue(scanNumber, out var scan)) continue;
                    if (!double.TryParse(f[iCalcMH].Trim(), out double calcMH)) continue;
                    if (!double.TryParse(f[iPepMass].Trim(), out double pepMassNoMod)) continue;

                    string sequence = f[iSeq].Trim();
                    if (string.IsNullOrWhiteSpace(sequence) || !sequence.All(char.IsLetter)) continue;

                    var glycanBox = OGlycanBoxes.FirstOrDefault(p => Math.Abs(p.Mass - (calcMH - protonMass - pepMassNoMod)) < 0.02);
                    if (glycanBox == null) continue;

                    var peptide = new Protein(sequence, "byonic")
                        .Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>()).FirstOrDefault();
                    if (peptide == null) continue;

                    var products = new List<Product>();
                    peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
                    var childBoxes = GlycanBox.BuildChildOGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds).ToArray();
                    string boxTargetMotif = GlycanBox.GlobalOGlycans[glycanBox.ModIds[0]].Target.ToString();

                    var targetsOnly = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
                    var targetKeys = targetsOnly.Where(p => p.Value == boxTargetMotif).Select(p => p.Key).ToList();
                    if (targetKeys.Count < 2) continue;

                    foreach (var residue in candidateResidues)
                    {
                        var modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
                        var decoySites = GlycoSpectralMatch.AddDecoyModSites(modPos, peptide, new[] { residue }, boxTargetMotif);
                        if (decoySites.Count == 0) continue;
                        if (!GraphCheck(modPos, glycanBox)) continue;

                        var graph = new LocalizationGraph(modPos, glycanBox, childBoxes, -1);
                        LocalizationGraph.LocalizeOGlycan(graph, scan, commonParameters.ProductMassTolerance, products);
                        if (graph.TotalScore < 2.0) continue;

                        var routes = LocalizationGraph.GetAllPaths_CalP(graph, 0.1, products.Count);
                        if (routes.Count == 0) continue;
                        var pairs = routes.SelectMany(p => p.ModSitePairs).Distinct().ToList();
                        LocalizationGraph.CalProbabilityForModSitePair(routes, pairs);
                        var winner = pairs.GroupBy(p => p.SiteIndex)
                            .Select(g => new { Site = g.Key, P = g.Sum(x => x.Probability) })
                            .OrderByDescending(x => x.P).First();

                        string sample = run.Key.Contains("StcEmix") ? "StcEmix"
                            : run.Key.Contains("HEK293") ? "HEK293" : "GlycoPepMix";
                        string key = sample + "|" + residue;

                        Bump(identifications, key, 1);
                        Bump(sitesOffered, key, decoySites.Count);
                        if (decoySites.Contains(winner.Site)) Bump(wins, key, 1);
                        // Chance expectation D/(T+D): win rate per SITE is confounded by how many sites a
                        // residue contributes, so the comparable quantity is wins against chance.
                        distanceSum.TryGetValue(key, out double dv);
                        distanceSum[key] = dv + decoySites.Count / (double)(modPos.Count);
                    }
                }
            }

            var samples = new[] { "StcEmix", "GlycoPepMix", "HEK293" };
            var csvRows = new List<string> { "Sample,Residue,Identifications,SitesOffered,Wins,WinRatePerSite,MeanDistanceToNearestRealSite" };

            double? Rate(string sample, string residue)
            {
                string k = sample + "|" + residue;
                if (!sitesOffered.TryGetValue(k, out int s) || s == 0) return null;
                wins.TryGetValue(k, out int w);
                return w / (double)s;
            }

            TestContext.WriteLine("");
            TestContext.WriteLine("Observed decoy wins / chance expectation, by sample (idents in brackets).");
            TestContext.WriteLine(">1 means the decoy wins MORE often than chance alone would give.");
            TestContext.WriteLine("Residue        StcEmix            GlycoPepMix              HEK293");
            foreach (var r in candidateResidues.OrderBy(r => { distanceSum.TryGetValue("StcEmix|" + r, out double c); wins.TryGetValue("StcEmix|" + r, out int w); return c > 0 ? w / c : 99; }))
            {
                var cells = samples.Select(s =>
                {
                    string k = s + "|" + r;
                    if (!identifications.TryGetValue(k, out int id) || id == 0) return "        -      ";
                    wins.TryGetValue(k, out int w);
                    distanceSum.TryGetValue(k, out double chanceSum);
                    double enrich = chanceSum > 0 ? w / chanceSum : double.NaN;
                    return string.Format("{0,7:F2} ({1,4})", enrich, id);
                });
                TestContext.WriteLine($"{r,-7}  {string.Join("  ", cells)}");

                foreach (var s in samples)
                {
                    string k = s + "|" + r;
                    if (!sitesOffered.TryGetValue(k, out int n) || n == 0) continue;
                    wins.TryGetValue(k, out int w);
                    identifications.TryGetValue(k, out int id);
                    distanceSum.TryGetValue(k, out double dsum);
                    csvRows.Add($"{s},{r},{id},{n},{w},{w / (double)n:F4},{dsum / Math.Max(1, id):F3}");
                }
            }

            File.WriteAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "decoy_residue_scan.csv"), csvRows);
        }

        /// <summary>
        /// Reports the acquisition structure of a raw file: MS2 count, dissociation types, and whether
        /// precursors carry child scans. Needed before configuring a search on unfamiliar data -- O-Pair
        /// consumes a collision/electron pair, and the localization scan must be the electron-based child.
        /// </summary>
        [Test]
        [Explicit("Diagnostic. Point it at a local raw file.")]
        public static void Measure_AcquisitionStructure()
        {
            string rawPath = @"E:\CodeReview\localization\data\raw\MSV000083070_170919_11.raw";
            if (!File.Exists(rawPath)) Assert.Ignore("Raw file not present locally.");

            var commonParameters = new CommonParameters(dissociationType: DissociationType.EThcD, trimMsMsPeaks: false);
            var file = new MyFileManager(true).LoadFile(rawPath, commonParameters);
            var scans = MetaMorpheusTask.GetMs2Scans(file, rawPath, commonParameters).ToList();

            TestContext.WriteLine($"File: {Path.GetFileName(rawPath)}");
            TestContext.WriteLine($"MS2 scans (precursor-grouped): {scans.Count}");
            TestContext.WriteLine("");
            TestContext.WriteLine("Parent dissociation types:");
            foreach (var g in scans.GroupBy(s => s.TheScan.DissociationType?.ToString() ?? "(null)").OrderByDescending(g => g.Count()))
            {
                TestContext.WriteLine($"  {g.Key,-20} {g.Count(),7}");
            }

            TestContext.WriteLine("");
            TestContext.WriteLine("Child-scan counts per precursor:");
            foreach (var g in scans.GroupBy(s => s.ChildScans.Count).OrderBy(g => g.Key))
            {
                TestContext.WriteLine($"  {g.Key} child scan(s): {g.Count(),7}");
            }

            var withChildren = scans.Where(s => s.ChildScans.Count > 0).ToList();
            if (withChildren.Count > 0)
            {
                TestContext.WriteLine("");
                TestContext.WriteLine("Child dissociation types:");
                foreach (var g in withChildren.SelectMany(s => s.ChildScans)
                    .GroupBy(c => c.TheScan.DissociationType?.ToString() ?? "(null)").OrderByDescending(g => g.Count()))
                {
                    TestContext.WriteLine($"  {g.Key,-20} {g.Count(),7}");
                }
            }

            TestContext.WriteLine("");
            TestContext.WriteLine($"Precursor charge range: {scans.Min(s => s.PrecursorCharge)} - {scans.Max(s => s.PrecursorCharge)}");
            TestContext.WriteLine($"Precursor mass range  : {scans.Min(s => s.PrecursorMass):F0} - {scans.Max(s => s.PrecursorMass):F0}");
        }

        /// <summary>
        /// Repeats the decoy-residue scan on an independent O-glyco dataset (MSV000083070, human urine),
        /// driven by a MetaMorpheus glyco search rather than a deposited third-party table.
        /// <para>
        /// This is the transferability test M15 could not run: PXD017646 contains only one O-glyco sample,
        /// so the Q/I/F/A ranking there might be specific to StcE-digested mucins. Urine O-glycopeptides are
        /// a different sample, a different lab and a different sequence composition.
        /// </para>
        /// <para>
        /// Identifications come from MetaMorpheus itself here, which would be circular for an absolute FLR
        /// but is sound for a residue <i>ranking</i>: every residue is compared on the identical
        /// identification set, so any self-consistency bias applies equally to all of them and cannot
        /// manufacture an ordering.
        /// </para>
        /// </summary>
        [Test]
        [Explicit("Measurement harness. Requires the MSV000083070 pilot search output.")]
        public static void Measure_DecoyResidueScan_UrineIndependentDataset()
        {
            string rawPath = @"E:\CodeReview\localization\data\raw\MSV000083070_170919_11.raw";
            string psmtsv = @"E:\CodeReview\localization\results\MSV000083070_pilot\Task1GlycoSearchTask\oglyco.psmtsv";
            if (!File.Exists(rawPath) || !File.Exists(psmtsv)) Assert.Ignore("MSV000083070 pilot inputs not present.");

            var commonParameters = new CommonParameters(dissociationType: DissociationType.EThcD, trimMsMsPeaks: false);
            var file = new MyFileManager(true).LoadFile(rawPath, commonParameters);
            var scansByNumber = MetaMorpheusTask.GetMs2Scans(file, rawPath, commonParameters)
                .GroupBy(p => p.OneBasedScanNumber).ToDictionary(g => g.Key, g => g.First());

            var lines = File.ReadAllLines(psmtsv);
            var header = lines[0].Split('\t');
            int iSeq = Array.IndexOf(header, "Base Sequence");
            int iScan = Array.IndexOf(header, "Scan Number");
            int iGlycanMass = Array.IndexOf(header, "GlycanMass");
            int iQ = Array.IndexOf(header, "QValue");
            int iDecoy = Array.IndexOf(header, "Decoy/Contaminant/Target");
            Assert.That(Math.Min(iSeq, Math.Min(iScan, Math.Min(iGlycanMass, Math.Min(iQ, iDecoy)))), Is.GreaterThanOrEqualTo(0),
                "Required oglyco.psmtsv columns not found.");

            string[] candidateResidues = { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "V", "W" };
            var sitesOffered = candidateResidues.ToDictionary(r => r, r => 0);
            var wins = candidateResidues.ToDictionary(r => r, r => 0);
            var identifications = candidateResidues.ToDictionary(r => r, r => 0);
            // Win rate per SITE is confounded by how many sites a residue offers: a residue contributing five
            // decoy positions has each one diluted relative to a residue contributing one. The comparable
            // quantity is the per-identification win rate against what chance alone would give,
            // D/(T+D), summed per identification.
            var chanceExpectation = candidateResidues.ToDictionary(r => r, r => 0.0);
            int considered = 0, usable = 0;

            foreach (var line in lines.Skip(1).Where(l => !string.IsNullOrWhiteSpace(l)))
            {
                var f = line.Split('\t');
                if (f.Length <= Math.Max(iSeq, Math.Max(iScan, iGlycanMass))) continue;
                if (f[iDecoy].Trim() != "T") continue;
                if (!double.TryParse(f[iQ].Trim(), out double q) || q > 0.01) continue;
                considered++;

                string sequence = f[iSeq].Trim();
                if (string.IsNullOrWhiteSpace(sequence) || sequence.Contains("|") || !sequence.All(char.IsLetter)) continue;
                if (!int.TryParse(f[iScan].Trim(), out int scanNumber)) continue;
                if (!scansByNumber.TryGetValue(scanNumber, out var scan)) continue;
                if (!double.TryParse(f[iGlycanMass].Trim(), out double glycanMass)) continue;

                var glycanBox = OGlycanBoxes.FirstOrDefault(p => Math.Abs(p.Mass - glycanMass) < 0.02);
                if (glycanBox == null) continue;

                var peptide = new Protein(sequence, "mm")
                    .Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>()).FirstOrDefault();
                if (peptide == null) continue;

                var products = new List<Product>();
                peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
                var childBoxes = GlycanBox.BuildChildOGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds).ToArray();
                string boxTargetMotif = GlycanBox.GlobalOGlycans[glycanBox.ModIds[0]].Target.ToString();

                var targetsOnly = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
                if (targetsOnly.Count(p => p.Value == boxTargetMotif) < 2) continue;
                usable++;

                foreach (var residue in candidateResidues)
                {
                    var modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
                    var decoySites = GlycoSpectralMatch.AddDecoyModSites(modPos, peptide, new[] { residue }, boxTargetMotif);
                    if (decoySites.Count == 0) continue;
                    if (!GraphCheck(modPos, glycanBox)) continue;

                    var graph = new LocalizationGraph(modPos, glycanBox, childBoxes, -1);
                    LocalizationGraph.LocalizeOGlycan(graph, scan, commonParameters.ProductMassTolerance, products);
                    if (graph.TotalScore < 2.0) continue;

                    var routes = LocalizationGraph.GetAllPaths_CalP(graph, 0.1, products.Count);
                    if (routes.Count == 0) continue;
                    var pairs = routes.SelectMany(p => p.ModSitePairs).Distinct().ToList();
                    LocalizationGraph.CalProbabilityForModSitePair(routes, pairs);
                    var winner = pairs.GroupBy(p => p.SiteIndex)
                        .Select(g => new { Site = g.Key, P = g.Sum(x => x.Probability) })
                        .OrderByDescending(x => x.P).First();

                    identifications[residue]++;
                    sitesOffered[residue] += decoySites.Count;
                    if (decoySites.Contains(winner.Site)) wins[residue]++;
                    int targetsHere = modPos.Count - decoySites.Count;
                    chanceExpectation[residue] += decoySites.Count / (double)(targetsHere + decoySites.Count);
                }
            }

            TestContext.WriteLine($"GSMs at 1% FDR considered: {considered}; usable (multi-site, box matched): {usable}");
            TestContext.WriteLine("");
            TestContext.WriteLine("Residue  Idents  Sites  Wins  Win%/ident  Chance%  Obs/Chance");
            var rows = new List<string> { "Residue,Identifications,SitesOffered,Wins,WinRatePerIdent,ChanceRate,ObservedOverChance" };
            foreach (var r in candidateResidues.Where(r => identifications[r] > 0)
                .OrderBy(r => wins[r] / Math.Max(1e-9, chanceExpectation[r])))
            {
                double perIdent = wins[r] / (double)identifications[r];
                double chance = chanceExpectation[r] / identifications[r];
                double enrich = wins[r] / Math.Max(1e-9, chanceExpectation[r]);
                TestContext.WriteLine($"{r,-7} {identifications[r],7} {sitesOffered[r],6} {wins[r],5} {perIdent,11:P1} {chance,8:P1} {enrich,11:F2}");
                rows.Add($"{r},{identifications[r]},{sitesOffered[r]},{wins[r]},{perIdent:F4},{chance:F4},{enrich:F3}");
            }
            File.WriteAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "decoy_residue_scan_urine.csv"), rows);
        }

        /// <summary>
        /// Mirrors GlycoSearchEngine.GraphCheck: the peptide's candidate-site motifs must cover what the box
        /// requires. LocalizationGraph.LocalizeOGlycan dereferences its terminal node unguarded (line ~140),
        /// so it throws NullReferenceException if this is not checked first. The engine always checks; any
        /// harness driving the graph directly has to check too.
        /// </summary>
        private static bool GraphCheck(SortedDictionary<int, string> modPos, GlycanBox glycanBox)
        {
            if (modPos.Count < glycanBox.NumberOfMods)
            {
                return false;
            }

            var required = glycanBox.ModIds
                .Select(id => GlycanBox.GlobalOGlycans[id].Target.ToString())
                .GroupBy(m => m).ToDictionary(g => g.Key, g => g.Count());
            var available = modPos.Values.GroupBy(m => m).ToDictionary(g => g.Key, g => g.Count());

            return required.All(kv => available.TryGetValue(kv.Key, out int n) && n >= kv.Value);
        }

        /// <summary>
        /// The M5 replication, on the full StcEmix EThcD run from PXD017646 rather than a sliced fixture.
        /// <para>
        /// Peptide identifications come from the deposited Byonic glycoPSM table, not from a MetaMorpheus
        /// search, so the peptide/scan/glycan assignments are independent of O-Pair. For each identification
        /// with more than one candidate site, a single decoy is swept across every non-candidate position and
        /// its probability recorded against the number of matched backbone ions separating it from the real
        /// site. Aggregating over hundreds of identifications is what turns M5 into a result.
        /// </para>
        /// Requires the 395 MB raw file and the glycoPSM table; skips cleanly if either is absent.
        /// </summary>
        [Test]
        [Explicit("Measurement harness. Requires the PXD017646 raw file to be present locally.")]
        public static void Measure_BracketingSweep_FullRun_ByonicIdentifications()
        {
            string dataDirectory = @"E:\CodeReview\localization\data\raw";

            // Every StcEmix electron-based run present locally, each paired with its Byonic table.
            var runs = Directory.GetFiles(dataDirectory, "*StcEmix_35trig_*.raw")
                .Select(raw =>
                {
                    string stem = Path.GetFileNameWithoutExtension(raw);
                    int idx = stem.IndexOf("StcEmix_35trig_", StringComparison.Ordinal);
                    string table = Path.Combine(dataDirectory, stem.Substring(idx) + "_GlycoPSMs.txt");
                    return (Raw: raw, Table: table);
                })
                .Where(p => File.Exists(p.Table))
                .OrderBy(p => p.Raw).ToList();

            if (runs.Count == 0)
            {
                Assert.Ignore("No PXD017646 raw/glycoPSM pairs present locally.");
            }
            TestContext.WriteLine($"Runs available: {runs.Count}");

            var commonParameters = new CommonParameters(dissociationType: DissociationType.EThcD, trimMsMsPeaks: false);
            var allRows = new List<string> { "Run,Peptide,Scan,DecoyPos,Residue,Distance,SeparatingIons,DecoyProb,RealSiteProb,GraphScore" };
            int totalConsidered = 0, totalSwept = 0, totalNoBox = 0, totalNoScan = 0, totalNoEvidence = 0;

            foreach (var run in runs)
            {
            string rawPath = run.Raw;
            string tablePath = run.Table;
            string runName = Path.GetFileNameWithoutExtension(rawPath);

            var file = new MyFileManager(true).LoadFile(rawPath, commonParameters);
            var scansByNumber = MetaMorpheusTask.GetMs2Scans(file, rawPath, commonParameters)
                .GroupBy(p => p.OneBasedScanNumber).ToDictionary(g => g.Key, g => g.First());
            TestContext.WriteLine($"Loaded {scansByNumber.Count} MS2 scans from {Path.GetFileName(rawPath)}");

            var lines = File.ReadAllLines(tablePath);
            var header = lines[0].Split('\t');
            int iSeq = Array.IndexOf(header, "Sequence");
            int iScan = Array.IndexOf(header, "ScanNumber");
            int iFrag = Array.IndexOf(header, "Fragmentation");
            int iCalcMH = Array.IndexOf(header, "calcMH");
            int iPepMass = Array.IndexOf(header, "PepMassNoMod");
            int iSites = Array.IndexOf(header, "#GlycoSitesOnPep");

            const double protonMass = 1.00727646677;
            var rows = allRows;
            int considered = 0, swept = 0, skippedNoBox = 0, skippedNoScan = 0, skippedNoEvidence = 0;

            foreach (var line in lines.Skip(1).Where(l => !string.IsNullOrWhiteSpace(l)))
            {
                var f = line.Split('\t');
                if (f.Length <= iSites) continue;

                // Electron-based scans only: localization is scored on c/zDot ions.
                string fragmentation = f[iFrag].Trim();
                if (!fragmentation.Contains("ETD", StringComparison.OrdinalIgnoreCase)
                    && !fragmentation.Contains("EThcD", StringComparison.OrdinalIgnoreCase)) continue;

                // A single candidate site has nothing to localize (about half the corpus -- see M6).
                if (!int.TryParse(f[iSites].Trim(), out int byonicSites) || byonicSites < 2) continue;
                considered++;

                if (!int.TryParse(f[iScan].Trim(), out int scanNumber)) continue;
                if (!scansByNumber.TryGetValue(scanNumber, out var scan)) { skippedNoScan++; continue; }
                if (!double.TryParse(f[iCalcMH].Trim(), out double calcMH)) continue;
                if (!double.TryParse(f[iPepMass].Trim(), out double pepMassNoMod)) continue;

                string sequence = f[iSeq].Trim();
                if (string.IsNullOrWhiteSpace(sequence) || !sequence.All(char.IsLetter)) continue;

                double glycanMass = calcMH - protonMass - pepMassNoMod;
                var glycanBox = OGlycanBoxes.FirstOrDefault(p => Math.Abs(p.Mass - glycanMass) < 0.02);
                if (glycanBox == null) { skippedNoBox++; continue; }

                var peptide = new Protein(sequence, "byonic")
                    .Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>()).FirstOrDefault();
                if (peptide == null) continue;

                var products = new List<Product>();
                peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
                var childBoxes = GlycanBox.BuildChildOGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds).ToArray();
                string boxTargetMotif = GlycanBox.GlobalOGlycans[glycanBox.ModIds[0]].Target.ToString();

                var targetsOnly = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
                if (targetsOnly.Count(p => p.Value == boxTargetMotif) < 2) continue;
                if (!GraphCheck(targetsOnly, glycanBox)) { skippedNoBox++; continue; }
                int realSite = targetsOnly.First(p => p.Value == boxTargetMotif).Key;

                // Require the peptide/scan/box triple to carry real evidence, or the posterior is uniform and
                // the arm measures nothing (the M7 failure mode).
                var baseGraph = new LocalizationGraph(targetsOnly, glycanBox, childBoxes, -1);
                LocalizationGraph.LocalizeOGlycan(baseGraph, scan, commonParameters.ProductMassTolerance, products);
                if (baseGraph.TotalScore < 2.0) { skippedNoEvidence++; continue; }

                var matchedPositions = MetaMorpheusEngine.MatchFragmentIons(scan, products, commonParameters)
                    .Select(p => p.NeutralTheoreticalProduct.ResiduePosition).Distinct().OrderBy(p => p).ToList();

                swept++;

                for (int r = 0; r < peptide.BaseSequence.Length; r++)
                {
                    int siteKey = r + 2;
                    if (targetsOnly.ContainsKey(siteKey)) continue;

                    var modPos = new SortedDictionary<int, string>(targetsOnly.ToDictionary(p => p.Key, p => p.Value));
                    modPos[siteKey] = boxTargetMotif;

                    var graph = new LocalizationGraph(modPos, glycanBox, childBoxes, -1);
                    LocalizationGraph.LocalizeOGlycan(graph, scan, commonParameters.ProductMassTolerance, products);
                    var routes = LocalizationGraph.GetAllPaths_CalP(graph, 0.1, products.Count);
                    if (routes.Count == 0) continue;
                    var pairs = routes.SelectMany(p => p.ModSitePairs).Distinct().ToList();
                    LocalizationGraph.CalProbabilityForModSitePair(routes, pairs);
                    var bySite = pairs.GroupBy(p => p.SiteIndex).ToDictionary(g => g.Key, g => g.Sum(p => p.Probability));

                    bySite.TryGetValue(siteKey, out double decoyProbability);
                    bySite.TryGetValue(realSite, out double realProbability);

                    int low = Math.Min(siteKey, realSite);
                    int high = Math.Max(siteKey, realSite);
                    int separating = matchedPositions.Count(p => p >= low - 1 && p < high - 1);

                    rows.Add($"{runName},{sequence},{scanNumber},{siteKey},{peptide.BaseSequence[r]},{Math.Abs(siteKey - realSite)},{separating},{decoyProbability:F4},{realProbability:F4},{graph.TotalScore:F3}");
                }
            }

            TestContext.WriteLine($"  {runName}: considered {considered}, swept {swept} (no-box {skippedNoBox}, no-scan {skippedNoScan}, no-evidence {skippedNoEvidence})");
            totalConsidered += considered; totalSwept += swept;
            totalNoBox += skippedNoBox; totalNoScan += skippedNoScan; totalNoEvidence += skippedNoEvidence;
            } // end foreach run

            string csv = Path.Combine(TestContext.CurrentContext.TestDirectory, "bracketing_sweep_fullrun.csv");
            File.WriteAllLines(csv, allRows);
            TestContext.WriteLine("");
            TestContext.WriteLine($"TOTAL multi-site ETD identifications considered: {totalConsidered}");
            TestContext.WriteLine($"  skipped, scan not in raw : {totalNoScan}");
            TestContext.WriteLine($"  skipped, no matching box : {totalNoBox}");
            TestContext.WriteLine($"  skipped, no evidence     : {totalNoEvidence}");
            TestContext.WriteLine($"  SWEPT                    : {totalSwept}   (decoy arms: {allRows.Count - 1})");
            TestContext.WriteLine("Wrote " + csv);

            var parsed = allRows.Skip(1).Select(r => r.Split(',')).ToList();
            if (parsed.Count == 0) return;

            TestContext.WriteLine($"Unique peptide sequences: {parsed.Select(x => x[1]).Distinct().Count()}");
            TestContext.WriteLine("");
            TestContext.WriteLine("SeparatingIons  Arms  MeanDecoyProb  MedianDecoyProb  MaxDecoyProb  Frac>0.05");
            foreach (var group in parsed.GroupBy(x => Math.Min(int.Parse(x[6]), 6)).OrderBy(g => g.Key))
            {
                var probabilities = group.Select(x => double.Parse(x[7])).OrderBy(p => p).ToList();
                double median = probabilities[probabilities.Count / 2];
                string label = group.Key == 6 ? "6+" : group.Key.ToString();
                TestContext.WriteLine($"{label,14}  {probabilities.Count,4}  {probabilities.Average(),13:F4}  {median,15:F4}  {probabilities.Max(),12:F4}  {probabilities.Count(p => p > 0.05) / (double)probabilities.Count,9:F3}");
            }

            int decoyBeatsReal = parsed.Count(x => double.Parse(x[7]) > double.Parse(x[8]));
            int decoyOverCutoff = parsed.Count(x => double.Parse(x[7]) > 0.75);
            TestContext.WriteLine("");
            TestContext.WriteLine($"Decoy outscores the real site : {decoyBeatsReal} / {parsed.Count} = {decoyBeatsReal / (double)parsed.Count:P1}");
            TestContext.WriteLine($"Decoy exceeds the 0.75 cutoff : {decoyOverCutoff} / {parsed.Count} = {decoyOverCutoff / (double)parsed.Count:P1}");
        }

        /// <summary>
        /// Runs a real O-glyco search and summarises the written results by localization level, recording how
        /// many rows carry a site-specific probability. This is the before/after evidence for A1: today
        /// GlycoSearchTask only computes probabilities for Level1 and Level2, so Level3 rows -- the ambiguous
        /// ones an FLR is actually for -- come out empty.
        /// </summary>
        [Test]
        [Explicit("Measurement harness, not a correctness test. Run deliberately.")]
        public static void Measure_ProbabilityCoverage_ByLocalizationLevel()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TESTGlycoProbCoverage");
            if (Directory.Exists(outputFolder))
            {
                Directory.Delete(outputFolder, true);
            }
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Nett.Toml.ReadFile<GlycoSearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSnip.toml"),
                MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.WriteContaminants = false;

            var db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoProteinFASTA_7proteins.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"GlycoTestData\GlycoPepMix_snip.mzML");

            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) },
                new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder).Run();

            var psmtsv = Directory.GetFiles(outputFolder, "*.psmtsv", SearchOption.AllDirectories)
                .OrderByDescending(p => new FileInfo(p).Length).FirstOrDefault();
            Assert.That(psmtsv, Is.Not.Null, "No .psmtsv written by the glyco search.");
            TestContext.WriteLine("Reading " + psmtsv);

            var lines = File.ReadAllLines(psmtsv);
            var header = lines[0].Split('\t');
            int levelColumn = Array.IndexOf(header, "GlycanLocalizationLevel");
            // "Localized Glycans with ..." emits only pairs flagged Confident, which by construction excludes
            // every Level3 match. "All SiteSpecific Localization Probability" emits the full posterior over
            // candidate sites, so that is the column that shows whether a probability was computed at all.
            int confidentColumn = Array.IndexOf(header, "Localized Glycans with Peptide Site Specific Probability");
            int probabilityColumn = Array.IndexOf(header, "All SiteSpecific Localization Probability");
            Assert.That(levelColumn, Is.GreaterThanOrEqualTo(0), "GlycanLocalizationLevel column not found.");
            Assert.That(confidentColumn, Is.GreaterThanOrEqualTo(0), "Confident-pair probability column not found.");
            Assert.That(probabilityColumn, Is.GreaterThanOrEqualTo(0), "All-site probability column not found.");

            var byLevel = new SortedDictionary<string, (int Total, int WithProbability)>();
            var perRow = new List<string> { "Level\tProbabilityField" };
            foreach (var line in lines.Skip(1).Where(l => !string.IsNullOrWhiteSpace(l)))
            {
                var fields = line.Split('\t');
                if (fields.Length <= Math.Max(levelColumn, Math.Max(confidentColumn, probabilityColumn)))
                {
                    continue;
                }

                string level = string.IsNullOrWhiteSpace(fields[levelColumn]) ? "(blank)" : fields[levelColumn].Trim();
                bool hasProbability = !string.IsNullOrWhiteSpace(fields[probabilityColumn]);

                byLevel.TryGetValue(level, out var counts);
                byLevel[level] = (counts.Total + 1, counts.WithProbability + (hasProbability ? 1 : 0));
                perRow.Add(level + "\t" + fields[confidentColumn].Trim() + "\t" + fields[probabilityColumn].Trim());
            }

            TestContext.WriteLine("Level        Rows   WithProbability");
            foreach (var entry in byLevel)
            {
                TestContext.WriteLine($"{entry.Key,-10} {entry.Value.Total,6} {entry.Value.WithProbability,17}");
            }

            // Written so the values themselves can be diffed across a code change, not just the counts.
            perRow.Sort(StringComparer.Ordinal);
            string perRowPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "probability_coverage.tsv");
            File.WriteAllLines(perRowPath, perRow);
            TestContext.WriteLine("Wrote " + perRowPath);

            Directory.Delete(outputFolder, true);
        }

        /// <summary>
        /// Per-site probabilities on an ambiguous peptide, to see what the posterior actually looks like when
        /// the evidence does not discriminate. Today these are only computed for Level1/Level2 GSMs, so this
        /// is a preview of the population the FLR has to cover once that gate is lifted.
        /// </summary>
        [Test]
        [Explicit("Measurement harness, not a correctness test. Run deliberately.")]
        public static void Measure_PerSiteProbabilities_OnAmbiguousPeptide()
        {
            var commonParameters = new CommonParameters(dissociationType: DissociationType.ETD, trimMsMsPeaks: false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"GlycoTestData\181217_Fusion_(LC2)_NewObj_Serum_deSA_Jacalin_HRM_4h_ETD_HCD_DDA_mz(400_1200)_21707.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scan = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).First();

            foreach (int siteCount in new[] { 2, 4, 6 })
            {
                var peptide = MakePeptideWithSites(siteCount);
                var products = new List<Product>();
                peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

                var modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });
                var glycanBox = OGlycanBoxes.First(p => p.NumberOfMods == 1);
                var childBoxes = GlycanBox.BuildChildOGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds).ToArray();

                var localizationGraph = new LocalizationGraph(modPos, glycanBox, childBoxes, -1);
                LocalizationGraph.LocalizeOGlycan(localizationGraph, scan, commonParameters.ProductMassTolerance, products);

                var routes = LocalizationGraph.GetAllPaths_CalP(localizationGraph, 0.1, products.Count);
                var allPairs = routes.SelectMany(p => p.ModSitePairs).Distinct().ToList();
                LocalizationGraph.CalProbabilityForModSitePair(routes, allPairs);

                TestContext.WriteLine($"--- {siteCount} candidate sites, {routes.Count} routes ---");
                foreach (var pair in allPairs.OrderBy(p => p.SiteIndex))
                {
                    TestContext.WriteLine($"  site {pair.SiteIndex,3}  modId {pair.ModId,3}  probability {pair.Probability:F4}");
                }
                TestContext.WriteLine($"  probability sum = {allPairs.Sum(p => p.Probability):F4}");
            }
        }
    }
}
