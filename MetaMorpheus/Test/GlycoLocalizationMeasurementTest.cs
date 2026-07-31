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
