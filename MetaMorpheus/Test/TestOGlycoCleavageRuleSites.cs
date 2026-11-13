using EngineLayer.DatabaseLoading;
using EngineLayer.GlycoSearch;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// <see cref="CleavageRuleSite"/>: the residues where a protease's cuts say a glycan must, or must not, be, with what the
    /// Level 1 identifications say about each. Written as oglyco_cleavage_rule_sites.tsv when
    /// <see cref="GlycoSearchParameters.BreakLocalizationTiesByCleavageRule"/> is on.
    /// </summary>
    [TestFixture]
    public static class TestOGlycoCleavageRuleSites
    {
        // level1SiteKeys: the two-based site keys a Level 1 route places glycans on (key k is residue start + k - 2); none
        // means the PSM is not Level 1.
        private static (DigestionProduct, DigestionAgent, IEnumerable<int>) Psm(Protein protein, string protease, int start, int end, params int[] level1SiteKeys) =>
            (new PeptideWithSetModifications(protein, new DigestionParams(protease), start, end, CleavageSpecificity.Full, "", 0,
                new Dictionary<int, Modification>(), 0), ProteaseDictionary.Dictionary[protease],
             level1SiteKeys.Length == 0 ? null : level1SiteKeys);

        private static string Row(CleavageRuleSite s) =>
            $"{s.AminoAcid}{s.OneBasedResidue} {(s.IsForbidden ? "forbidden" : "required")} implied {s.ImpliedBy}, agreeing {s.ImpliersAgreeing}, " +
            $"against {s.ImpliersAgainst}, others with {s.OthersWithGlycan}, without {s.OthersWithout}";

        /// <summary>RPPITQSSL: StcE cuts T5-Q6|S7, so its P2, T5, needs a glycan, whichever side of the cut a peptide is on.</summary>
        [Test]
        public static void AStcECut_ImpliesP2_AndTheLevel1PsmsAreCounted()
        {
            var protein = new Protein("RPPITQSSL", "P1");
            var sites = CleavageRuleSite.Find(new[]
            {
                Psm(protein, "StcE", 1, 6, 6),    // ends at the cut, Level 1 with a glycan on T5 (key 6): agrees
                Psm(protein, "StcE", 1, 6, 5),    // ends at the cut, Level 1 with it on I4 instead: against
                Psm(protein, "StcE", 7, 9),       // starts at the cut, its T5 on the neighbour: implies, not Level 1
                Psm(protein, "StcE", 1, 9, 6),    // the whole protein, no cut: another Level 1 PSM with a glycan on T5
                Psm(protein, "StcE", 2, 9),       // no cut either, and not Level 1: not counted
            });
            Assert.That(sites.Select(Row), Is.EqualTo(new[] { "T5 required implied 3, agreeing 1, against 1, others with 1, without 0" }));
        }

        /// <summary>GSGSTGG: IMPa cuts before T5 only without a glycan on S4.</summary>
        [Test]
        public static void AnIMPaCut_ImpliesARequiredAndAForbiddenResidue()
        {
            var protein = new Protein("GSGSTGG", "P1");
            var sites = CleavageRuleSite.Find(new[] { Psm(protein, "IMPa", 1, 4, 5) }); // a glycan on S4 (key 5)
            Assert.That(sites.Select(Row), Is.EqualTo(new[]
            {
                "S4 forbidden implied 1, agreeing 0, against 1, others with 0, without 0",
                "T5 required implied 1, agreeing 0, against 0, others with 0, without 0",
            }));
        }

        [Test]
        public static void ACutTrypsinAlsoMakes_ImpliesNothing() =>
            Assert.That(CleavageRuleSite.Find(new[] { Psm(new Protein("GTKSGG", "P1"), "StcE-trypsin", 4, 6, 2) }), Is.Empty);

        /// <summary>
        /// End to end: the O-glyco search of the leukosialin spectrum (TTGSLEPSSGASGPQVSSVK, residues 208-227 of the test
        /// FASTA). StcE cut before it (...PPVTM|TT...), so StcE's P2, T206, on the neighbouring peptide, needs a glycan.
        /// </summary>
        [Test]
        [TestCase(true)]
        [TestCase(false)]
        public static void OGlycoSearch_WritesTheCleavageRuleSites_OnlyWithTheOptionOn(bool on)
        {
            string[] rows = SearchLeukosialin(on, settings => settings);
            Assert.That(rows, on ? Is.EqualTo(new[] { "P16150\tT206\tglycan required\t1\t0\t0\t0\t0" }) : Is.Null);
        }

        // The report's rows, or null when no report was written.
        private static string[] SearchLeukosialin(bool on, System.Func<string, string> editSettings)
        {
            string dataFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData");
            string settings = File.ReadAllText(Path.Combine(dataFolder, "GlycoSearchTaskconfigOGlycoTest_Run.toml"));
            Assert.That(settings, Does.Contain("Protease = \"StcE-trypsin\"").And.Contain("MaxMissedCleavages = 5"), "premise: the settings file this test edits");
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, $"{nameof(TestOGlycoCleavageRuleSites)}_{on}_{settings.GetHashCode() ^ editSettings(settings).GetHashCode()}");
            string settingsPath = outputFolder + ".toml";
            try
            {
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                Directory.CreateDirectory(outputFolder);
                File.WriteAllText(settingsPath, editSettings(settings));
                var task = MetaMorpheusTask.ReadTaskTomlWithLowResFallback<GlycoSearchTask>(settingsPath);
                task._glycoSearchParameters.BreakLocalizationTiesByCleavageRule = on;
                new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) },
                    new List<string> { Path.Combine(dataFolder, "2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf") },
                    new List<DbForTask> { new DbForTask(Path.Combine(dataFolder, "P16150.fasta"), false) },
                    outputFolder).Run();
                Assert.That(File.Exists(Path.Combine(outputFolder, "Task", "oglyco.psmtsv")), Is.True, "premise: the glycopeptide is identified");

                string report = Path.Combine(outputFolder, "Task", "oglyco_cleavage_rule_sites.tsv");
                return File.Exists(report) ? File.ReadAllLines(report).Skip(1).ToArray() : null;
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
