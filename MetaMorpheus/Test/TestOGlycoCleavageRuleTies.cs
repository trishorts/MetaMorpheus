using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.GlycoSearch;
using Nett;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// <see cref="GlycoSearchParameters.BreakLocalizationTiesByCleavageRule"/>: when the spectrum scores several glycan
    /// placements equally, keep those the glycoprotease's cleavage rule allows. A placement the spectrum prefers is never
    /// overridden, and identifications do not change.
    /// </summary>
    [TestFixture]
    public static class TestOGlycoCleavageRuleTies
    {
        private static PeptideWithSetModifications Peptide(string sequence, string protease, int start, int end) =>
            new(new Protein(sequence, "P1"), new DigestionParams(protease), start, end, CleavageSpecificity.Full, "", 0,
                new Dictionary<int, Modification>(), 0);

        // A route placing one glycan on each two-based site key.
        private static Route RouteOn(params int[] keys)
        {
            var route = new Route();
            foreach (int key in keys)
                route.AddPos(key, 0, false);
            return route;
        }

        private static List<int[]> Kept(List<Route> routes, string protease, params PeptideWithSetModifications[] peptides) =>
            GlycoSpectralMatch.RoutesMeetingCleavageConditions(routes, peptides, ProteaseDictionary.Dictionary[protease])
                .Select(r => r.ModSitePairs.Select(p => p.SiteIndex).ToArray()).ToList();

        /// <summary>TSTGGG starts where OpeRATOR cut before T5, so T5 (key 2) needs a glycan.</summary>
        [Test]
        public static void OfTiedRoutes_KeepsThoseTheRuleAllows()
        {
            var peptide = Peptide("GGGGTSTGGG", "OpeRATOR", 5, 10);
            var routes = new List<Route> { RouteOn(2), RouteOn(3), RouteOn(4) };
            Assert.That(Kept(routes, "OpeRATOR", peptide), Is.EqualTo(new[] { new[] { 2 } }));
        }

        [Test]
        public static void WhenNoTiedRouteMeetsTheRule_KeepsThemAll()
        {
            var peptide = Peptide("GGGGTSTGGG", "OpeRATOR", 5, 10);
            var routes = new List<Route> { RouteOn(3), RouteOn(4) };
            Assert.That(Kept(routes, "OpeRATOR", peptide), Is.EqualTo(new[] { new[] { 3 }, new[] { 4 } }));
        }

        /// <summary>After K4, trypsin explains the cut, so OpeRATOR-trypsin|P imposes nothing.</summary>
        [Test]
        public static void WithoutAConditionInsideThePeptide_KeepsThemAll()
        {
            var peptide = Peptide("GGGKTSTGGG", "OpeRATOR-trypsin|P", 5, 10);
            var routes = new List<Route> { RouteOn(2), RouteOn(3) };
            Assert.That(Kept(routes, "OpeRATOR-trypsin|P", peptide), Is.EqualTo(new[] { new[] { 2 }, new[] { 3 } }));
        }

        /// <summary>GSGS ends where IMPa cut before T5, which it does only without a glycan on S4 (key 5).</summary>
        [Test]
        public static void AForbiddenSite_RemovesTheRoutesThatUseIt()
        {
            var peptide = Peptide("GSGSTGG", "IMPa", 1, 4);
            var routes = new List<Route> { RouteOn(3), RouteOn(5) };
            Assert.That(Kept(routes, "IMPa", peptide), Is.EqualTo(new[] { new[] { 3 } }));
        }

        /// <summary>A PSM shared by two peptides whose rules differ is left alone.</summary>
        [Test]
        public static void WhenTheCandidatePeptidesDisagree_KeepsThemAll()
        {
            var obligated = Peptide("GGGGTSTGGG", "OpeRATOR-trypsin|P", 5, 10);
            var free = Peptide("GGGKTSTGGG", "OpeRATOR-trypsin|P", 5, 10);
            var routes = new List<Route> { RouteOn(2), RouteOn(3) };
            Assert.That(Kept(routes, "OpeRATOR-trypsin|P", obligated, free), Is.EqualTo(new[] { new[] { 2 }, new[] { 3 } }));
        }

        /// <summary>
        /// PEPTISSTIDE, uncut by IMPa ([D]|T,[D]|S; a glycan needed at P1', forbidden at P1): naked it has 4 sites, before T4,
        /// S6, S7 and T8. With glycans on S6 and S7 IMPa could have cut only before S6 (S7's P1, S6, carries one); with
        /// glycans on T4 and T8, only before those two.
        /// </summary>
        [Test]
        [TestCase(new[] { 7, 8 }, 1)]
        [TestCase(new[] { 5, 9 }, 2)]
        public static void MissedCleavages_CountOnlyTheSitesThePlacedGlycansAllow(int[] keys, int expected)
        {
            var peptide = new Protein("PEPTISSTIDE", "P1").Digest(new DigestionParams("IMPa", maxMissedCleavages: 10, minPeptideLength: 1),
                new List<Modification>(), new List<Modification>()).Single(p => p.OneBasedStartResidue == 1 && p.OneBasedEndResidue == 11);
            Assert.That(peptide.MissedCleavages, Is.EqualTo(4), "premise: the naked count");
            Assert.That(GlycoPeptides.MissedCleavagesGivenGlycans(peptide, RouteOn(keys), ProteaseDictionary.Dictionary["IMPa"]), Is.EqualTo(expected));
        }

        [Test]
        public static void MissedCleavages_ASiteWithoutConditionsAlwaysCounts()
        {
            var peptide = Peptide("PEPTKIDERSEQK", "trypsin", 1, 13);
            Assert.That(GlycoPeptides.MissedCleavagesGivenGlycans(peptide, RouteOn(5), ProteaseDictionary.Dictionary["trypsin"]), Is.EqualTo(2));
        }

        /// <summary>
        /// StcE (S/T-X|S/T, a glycan at P2) cutting after Q: when P2 is inside the peptide its placed glycan decides; when P2
        /// is the residue before the peptide, an S or T there is taken to carry one.
        /// </summary>
        [Test]
        public static void MissedCleavages_StcEsP2_InsideIsAsPlaced_OutsideIsAssumed()
        {
            var agent = ProteaseDictionary.Dictionary["StcE"];
            var inside = Peptide("GTQSGG", "StcE", 2, 6); // T2 Q3 | S4
            Assert.That(GlycoPeptides.MissedCleavagesGivenGlycans(inside, RouteOn(4), agent), Is.EqualTo(0), "T2 bare");
            Assert.That(GlycoPeptides.MissedCleavagesGivenGlycans(inside, RouteOn(2), agent), Is.EqualTo(1), "T2 glycosylated");
            var outside = Peptide("TQSGG", "StcE", 2, 5);  // T1 | Q2 | S3: P2 is T1, on the neighbour
            Assert.That(GlycoPeptides.MissedCleavagesGivenGlycans(outside, RouteOn(3), agent), Is.EqualTo(1));
        }

        [Test]
        public static void TheOption_IsOffByDefault_RoundTripsThroughToml_AndIsOffInAnOlderToml()
        {
            Assert.That(new GlycoSearchTask()._glycoSearchParameters.BreakLocalizationTiesByCleavageRule, Is.False);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"{nameof(TestOGlycoCleavageRuleTies)}.toml");
            try
            {
                var task = new GlycoSearchTask();
                task._glycoSearchParameters.BreakLocalizationTiesByCleavageRule = true;
                Toml.WriteFile(task, path, MetaMorpheusTask.tomlConfig);
                Assert.That(Toml.ReadFile<GlycoSearchTask>(path, MetaMorpheusTask.tomlConfig)._glycoSearchParameters.BreakLocalizationTiesByCleavageRule, Is.True);

                var lines = File.ReadAllLines(path);
                var older = lines.Where(line => !line.TrimStart().StartsWith(nameof(GlycoSearchParameters.BreakLocalizationTiesByCleavageRule))).ToArray();
                Assert.That(older.Length, Is.LessThan(lines.Length), "premise: the option was written");
                File.WriteAllLines(path, older);
                Assert.That(Toml.ReadFile<GlycoSearchTask>(path, MetaMorpheusTask.tomlConfig)._glycoSearchParameters.BreakLocalizationTiesByCleavageRule, Is.False);
            }
            finally
            {
                if (File.Exists(path))
                    File.Delete(path);
            }
        }

        /// <summary>
        /// End to end: the leukosialin spectrum (TTGSLEPSSGASGPQVSSVK, see <see cref="GlycoSemiSpecificSearchTests"/>)
        /// searched with OpeRATOR-trypsin|P, matching only b/y ions, so the spectrum cannot place its glycans and many
        /// routes tie. The peptide follows a Met (...GPPVTM|TTGS...), so OpeRATOR's rule requires a glycan on T1. The
        /// spectrum is from a StcE digest; it stands in for an OpeRATOR one here to exercise the mechanism.
        /// </summary>
        [Test]
        public static void OGlycoSearch_BreaksTiesByTheRule_AndReportsTheRulesPlacement()
        {
            var off = SearchLeukosialin(false);
            var on = SearchLeukosialin(true);

            Assert.That(off.Total, Is.GreaterThan(1), "premise: the spectrum alone ties several placements");
            Assert.That(off.Routes.Any(r => !r.Contains(1)), Is.True, "premise: some tied placements leave T1 bare");
            Assert.That(off.Placed, Does.Not.Contain(1), "premise: O-Pair reports a placement leaving T1 bare");
            Assert.That(off.MissedCleavages, Is.EqualTo(7), "premise: the naked count, every internal S/T an OpeRATOR site");

            Assert.That(on.Routes.All(r => r.Contains(1)), Is.True, "every kept placement has a glycan on T1");
            Assert.That(on.Total, Is.LessThan(off.Total));
            Assert.That(on.Placed, Does.Contain(1), "the reported placement is one the rule allows");
            Assert.That(on.ComposedOfPlaced, Is.True, "the composition column describes the reported placement");
            // OpeRATOR could have cut before each glycosylated S/T after the first residue, and nowhere else inside.
            Assert.That(on.MissedCleavages, Is.EqualTo(on.Placed.Count(position => position > 1)));
            Assert.That(on.QValue, Is.LessThanOrEqualTo(0.01), "still identified");
        }

        // "H2N2A1" + "H1N1" -> "A1H3N3": monosaccharide counts, so two compositions compare regardless of grouping.
        private static string SumOf(IEnumerable<string> compositions) => string.Concat(compositions
            .SelectMany(c => Regex.Matches(c, @"([A-Z])(\d+)").Select(m => (m.Groups[1].Value, int.Parse(m.Groups[2].Value))))
            .GroupBy(x => x.Item1).OrderBy(g => g.Key).Select(g => g.Key + g.Sum(x => x.Item2)));

        // The routes as written (the first 10 at most) and how many there were; the reported placement's one-based glycan
        // positions; whether the composition column is their sum; and the PSM's missed cleavages and q-value.
        private static (List<HashSet<int>> Routes, int Total, HashSet<int> Placed, bool ComposedOfPlaced, int MissedCleavages, double QValue)
            SearchLeukosialin(bool breakTies) => SearchLeukosialin(breakTies, 1).Single();

        /// <summary>
        /// Two copies of the spectrum are searched at once, which builds the engine at the task's other construction site (the
        /// one file above goes the sequential way).
        /// </summary>
        [Test]
        public static void OGlycoSearch_ReportsTheRulesPlacement_WhenFilesAreSearchedTogether()
        {
            var psms = SearchLeukosialin(true, 2);
            Assert.That(psms.Count, Is.EqualTo(2), "premise: one PSM per file");
            Assert.That(psms.All(p => p.Placed.Contains(1)), Is.True);
        }

        // One entry per PSM of the glycopeptide, over this many copies of the spectrum.
        private static List<(List<HashSet<int>> Routes, int Total, HashSet<int> Placed, bool ComposedOfPlaced, int MissedCleavages, double QValue)>
            SearchLeukosialin(bool breakTies, int files)
        {
            string dataFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData");
            string settings = File.ReadAllText(Path.Combine(dataFolder, "GlycoSearchTaskconfigOGlycoTest_Run.toml"));
            Assert.That(settings, Does.Contain("Protease = \"StcE-trypsin\"").And.Contain("DissociationType = \"EThcD\"")
                .And.Contain("MS2ChildScanDissociationType = \"Unknown\"").And.Contain("MaxMissedCleavages = 5"), "premise: the settings file this test edits");
            settings = settings.Replace("Protease = \"StcE-trypsin\"", "Protease = \"OpeRATOR-trypsin|P\"") // also SpecificProtease
                .Replace("\nDissociationType = \"EThcD\"", "\nDissociationType = \"HCD\"")
                .Replace("MaxMissedCleavages = 5", "MaxMissedCleavages = 10"); // OpeRATOR's |S,|T put 7 sites inside the peptide

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, $"{nameof(TestOGlycoCleavageRuleTies)}_{breakTies}_{files}");
            string settingsPath = outputFolder + ".toml";
            try
            {
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                Directory.CreateDirectory(outputFolder);
                File.WriteAllText(settingsPath, settings);
                var task = MetaMorpheusTask.ReadTaskTomlWithLowResFallback<GlycoSearchTask>(settingsPath);
                Assert.That(task.CommonParameters.DissociationType, Is.EqualTo(MassSpectrometry.DissociationType.HCD), "premise: b/y ions only");
                task._glycoSearchParameters.BreakLocalizationTiesByCleavageRule = breakTies;

                string spectra = Path.Combine(dataFolder, "2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
                var spectraFiles = files == 1 ? new List<string> { spectra } : Enumerable.Range(0, files).Select(i =>
                {
                    string copy = Path.Combine(outputFolder, $"copy{i}.mgf");
                    File.Copy(spectra, copy);
                    return copy;
                }).ToList();
                new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, spectraFiles,
                    new List<DbForTask> { new DbForTask(Path.Combine(dataFolder, "P16150.fasta"), false) },
                    Path.Combine(outputFolder, "out")).Run();

                string prose = File.ReadAllText(Directory.GetFiles(Path.Combine(outputFolder, "out"), "AutoGeneratedManuscriptProse.txt", SearchOption.AllDirectories).Single());
                Assert.That(prose, files == 1 ? Does.Contain("spectra files searched one after another") : Does.Contain("spectra files searched in parallel"),
                    "premise: which engine construction site the search used");
                string oGlycoPath = Path.Combine(outputFolder, "out", "Task", "oglyco.psmtsv");
                Assert.That(File.Exists(oGlycoPath), Is.True, "premise: the search identifies the glycopeptide");
                var psms = SpectrumMatchTsvReader.ReadGlycoPsmTsv(oGlycoPath, out var errors).Where(p => p.BaseSeq == "TTGSLEPSSGASGPQVSSVK").ToList();
                Assert.That(errors, Is.Empty);

                return psms.Select(psm =>
                {
                    // "{@box[1-3,2-3,8-5]}{@box[...]}... In Total:N Paths": each route's one-based peptide positions and glycan
                    // ids, the first 10 only, then the count when there were more (GlycoSpectralMatch.AllLocalizationInfo).
                    var routes = Regex.Matches(psm.AllPotentialGlycanLocalization, @"\[([^\]]*)\]")
                        .Select(m => m.Groups[1].Value.Split(',').Select(pair => int.Parse(pair.Split('-')[0])).ToHashSet()).ToList();
                    var total = Regex.Match(psm.AllPotentialGlycanLocalization, @"In Total:(\d+)");
                    var placed = new HashSet<int>();
                    var glycans = new List<string>();
                    int position = 0;
                    foreach (Match m in Regex.Matches(psm.FullSequence, @"\[(O-linked glycosylation:([^ \]]+))?[^\]]*\]|[A-Z]"))
                    {
                        if (m.Groups[2].Success) { placed.Add(position); glycans.Add(m.Groups[2].Value); }
                        else if (!m.Value.StartsWith("[")) position++;
                    }
                    return (routes, total.Success ? int.Parse(total.Groups[1].Value) : routes.Count, placed,
                        SumOf(glycans) == SumOf(new[] { psm.GlycanComposition }), int.Parse(psm.MissedCleavage), psm.QValue);
                }).ToList();
            }
            finally
            {
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                if (File.Exists(settingsPath))
                    File.Delete(settingsPath);
            }
        }
    }
}
