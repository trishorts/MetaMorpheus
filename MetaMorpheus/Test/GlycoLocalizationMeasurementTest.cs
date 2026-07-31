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
