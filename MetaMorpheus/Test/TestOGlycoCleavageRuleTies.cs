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
        public static void OGlycoSearch_BreaksTiesByTheRule_AndLeavesTheIdentificationAlone()
        {
            var off = SearchLeukosialin(false);
            var on = SearchLeukosialin(true);

            Assert.That(off.Total, Is.GreaterThan(1), "premise: the spectrum alone ties several placements");
            Assert.That(off.Routes.Any(r => !r.Contains(1)), Is.True, "premise: some tied placements leave T1 bare");
            Assert.That(on.Routes.All(r => r.Contains(1)), Is.True, "every kept placement has a glycan on T1");
            Assert.That(on.Total, Is.LessThan(off.Total));
            Assert.That(on.Score, Is.EqualTo(off.Score), "the identification does not change");
            Assert.That(on.QValue, Is.EqualTo(off.QValue));
        }

        // The routes as written (the first 10 at most), and how many there were.
        private static (List<HashSet<int>> Routes, int Total, double Score, double QValue) SearchLeukosialin(bool breakTies)
        {
            string dataFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData");
            string settings = File.ReadAllText(Path.Combine(dataFolder, "GlycoSearchTaskconfigOGlycoTest_Run.toml"));
            Assert.That(settings, Does.Contain("Protease = \"StcE-trypsin\"").And.Contain("DissociationType = \"EThcD\"")
                .And.Contain("MS2ChildScanDissociationType = \"Unknown\"").And.Contain("MaxMissedCleavages = 5"), "premise: the settings file this test edits");
            settings = settings.Replace("Protease = \"StcE-trypsin\"", "Protease = \"OpeRATOR-trypsin|P\"") // also SpecificProtease
                .Replace("\nDissociationType = \"EThcD\"", "\nDissociationType = \"HCD\"")
                .Replace("MaxMissedCleavages = 5", "MaxMissedCleavages = 10"); // OpeRATOR's |S,|T put 7 sites inside the peptide

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, $"{nameof(TestOGlycoCleavageRuleTies)}_{breakTies}");
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

                new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) },
                    new List<string> { Path.Combine(dataFolder, "2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf") },
                    new List<DbForTask> { new DbForTask(Path.Combine(dataFolder, "P16150.fasta"), false) },
                    outputFolder).Run();

                string oGlycoPath = Path.Combine(outputFolder, "Task", "oglyco.psmtsv");
                Assert.That(File.Exists(oGlycoPath), Is.True, "premise: the search identifies the glycopeptide");
                var psm = SpectrumMatchTsvReader.ReadGlycoPsmTsv(oGlycoPath, out var errors).Single(p => p.BaseSeq == "TTGSLEPSSGASGPQVSSVK");
                Assert.That(errors, Is.Empty);

                // "{@box[1-3,2-3,8-5]}{@box[...]}... In Total:N Paths": each route's one-based peptide positions and glycan
                // ids, the first 10 only, then the count when there were more (GlycoSpectralMatch.AllLocalizationInfo).
                var routes = Regex.Matches(psm.AllPotentialGlycanLocalization, @"\[([^\]]*)\]")
                    .Select(m => m.Groups[1].Value.Split(',').Select(pair => int.Parse(pair.Split('-')[0])).ToHashSet()).ToList();
                var total = Regex.Match(psm.AllPotentialGlycanLocalization, @"In Total:(\d+)");
                return (routes, total.Success ? int.Parse(total.Groups[1].Value) : routes.Count, psm.Score, psm.QValue);
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
