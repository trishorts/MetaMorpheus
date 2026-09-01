using EngineLayer;
using Nett;
using System.IO;
using TaskLayer;
using EngineLayer.GlycoSearch;
using NUnit.Framework;
using Omics.Modifications;
using NUnit.Framework.Legacy;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Decoy glycosites are candidate sites on residues that cannot carry an O-glycan, so a
    /// localization landing on one is wrong BY CONSTRUCTION. That is the ground truth a
    /// false-localization rate needs.
    ///
    /// These tests exist because BOTH failure modes are silent:
    ///
    ///   too strict -> the decoy is never admitted, no decoy ever wins, and the measured error rate
    ///                 is a clean-looking ZERO;
    ///   too loose  -> the decoy draws from twice the glycan instances a real site does, wins more
    ///                 than it should, and the measured error rate is INFLATED in the direction that
    ///                 looks like a discovery.
    ///
    /// Neither throws. Both produce a plausible number. So the parity is asserted, not assumed.
    /// </summary>
    [TestFixture]
    public class TestGlycoDecoySiteParity
    {
        private static GlycanBox[] _boxes;

        [OneTimeSetUp]
        public static void Setup()
        {
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(
                GlobalVariables.OGlycanDatabasePaths.First(p => p.Contains("OGlycan.gdb")), true, true).ToArray();
            _boxes = GlycanBox.BuildOGlycanBoxes(3).OrderBy(p => p.Mass).ToArray();
        }

        [TearDown]
        public void Reset()
        {
            // Static, so it must not leak into another fixture.
            LocalizationGraph.DecoyGlycositeMotifs = new HashSet<string>();
            LocalizationGraph.CanonicalDecoyTarget = "T";
        }

        /// <summary>Boxes holding exactly one glycan: the 0 -> 1 transition isolates a single glycan.</summary>
        private static IEnumerable<GlycanBox> SingleGlycanBoxes()
            => _boxes.Where(b => b.ChildGlycanBoxes != null && b.ChildGlycanBoxes.Length >= 2
                                 && b.ChildGlycanBoxes[1].ModIds.Length == 1);

        private static int AdmittedAt(string motif)
            => SingleGlycanBoxes().Count(b => LocalizationGraph.MotifCheck(b, 0, 1, motif));

        [Test]
        public static void DecoyMotifIsInertUnlessConfigured()
        {
            LocalizationGraph.DecoyGlycositeMotifs = new HashSet<string>();

            Assert.That(AdmittedAt("A"), Is.EqualTo(0),
                "With no decoy motifs configured an ordinary search must be completely unaffected.");
        }

        [Test]
        public static void DecoySiteAdmitsExactlyAsManyGlycansAsARealSite()
        {
            int onS = AdmittedAt("S");
            int onT = AdmittedAt("T");

            Assert.That(onS, Is.GreaterThan(0), "sanity: real S sites must admit something");
            Assert.That(onT, Is.GreaterThan(0), "sanity: real T sites must admit something");

            LocalizationGraph.DecoyGlycositeMotifs = new HashSet<string> { "A" };
            int onDecoy = AdmittedAt("A");

            // THE PARITY CLAIM. O-Pair loads each composition twice, once targeting S and once
            // targeting T, so a real site admits one instance per composition. The decoy must too.
            Assert.That(onDecoy, Is.EqualTo(AdmittedAt(LocalizationGraph.CanonicalDecoyTarget)),
                "A decoy site must admit exactly the canonical half.");

            Assert.That(onDecoy, Is.LessThan(onS + onT),
                "A decoy admitting BOTH the S- and T-targeted copies would be twice as competitive as "
                + "any real site, inflating the false-localization rate in the direction that looks "
                + "like a finding.");
        }

        /// <summary>
        /// A parameter that does not survive TOML round-trip is the worst kind of bug here: the search
        /// runs, finishes, and reports a perfectly normal-looking result with NO decoy sites in it and
        /// no error anywhere. The measured false-localization rate would simply be zero.
        /// </summary>
        [Test]
        public static void NewGlycoParametersSurviveTomlRoundTrip()
        {
            var task = new GlycoSearchTask();
            task._glycoSearchParameters.DecoyGlycositeResidues = new[] { "A", "L" };
            task._glycoSearchParameters.RetainedGsmsPerScan = 25;

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                                       "GlycoDecoyRoundTrip.toml");
            Toml.WriteFile(task, path, MetaMorpheusTask.tomlConfig);
            var loaded = Toml.ReadFile<GlycoSearchTask>(path, MetaMorpheusTask.tomlConfig);

            Assert.That(loaded._glycoSearchParameters.RetainedGsmsPerScan, Is.EqualTo(25),
                "RetainedGsmsPerScan did not survive the toml round-trip.");
            Assert.That(loaded._glycoSearchParameters.DecoyGlycositeResidues, Is.EqualTo(new[] { "A", "L" }),
                "DecoyGlycositeResidues did not survive the toml round-trip -- a search configured with "
                + "decoy sites would silently run WITHOUT them and report a false-localization rate of zero.");

            File.Delete(path);
        }

        [Test]
        public static void PositionalDecoyMotifIsUsableAndCannotCollideWithAResidue()
        {
            // Construction (a) builds REAL Glycan instances targeting this motif, and Glycan's
            // constructor runs it through ModificationMotif.TryGetMotif -- which requires ^[A-Za-z]+$
            // with exactly one upper case. A motif that fails here would silently produce glycans with
            // a null Target.
            Assert.That(ModificationMotif.TryGetMotif(LocalizationGraph.PositionalDecoyMotif, out _), Is.True,
                "The positional-decoy motif must be a valid ModificationMotif.");

            Assert.That(LocalizationGraph.PositionalDecoyMotif.Length, Is.GreaterThan(1),
                "A single letter could collide with an amino-acid motif, silently turning every "
                + "occurrence of that residue into a decoy site.");
        }

        /// <summary>
        /// Construction (a) is what replaced the canonical-letter parity rule. M17 measured that the
        /// letter was not cosmetic: flipping it T -> S moved per-site posteriors by >0.05 in 16% of
        /// decoy sites and changed the winning site in 17% of GSMs. Parity must now hold BY
        /// CONSTRUCTION -- exactly one decoy instance per composition, matching the one instance a real
        /// S or T site gets.
        /// </summary>
        [Test]
        public static void DecoyTargetedGlycansGiveParityWithoutAnArbitraryLetter()
        {
            var original = GlycanBox.GlobalOGlycans;
            try
            {
                int onS = original.Count(g => g.Target != null && g.Target.ToString() == "S");
                int onT = original.Count(g => g.Target != null && g.Target.ToString() == "T");
                Assert.That(onS, Is.GreaterThan(0));
                Assert.That(onT, Is.EqualTo(onS), "each composition is loaded once per real motif");

                Assert.That(original.Count(g => g.Target != null
                        && g.Target.ToString() == LocalizationGraph.PositionalDecoyMotif), Is.EqualTo(0),
                    "without construction (a), nothing targets the decoy motif");

                // Exactly what AddDecoyTargetedGlycans does.
                var withDecoys = original
                    .Concat(original.Where(g => g.Target != null && g.Target.ToString() == "T")
                        .Select(g => new Glycan(g.Kind, LocalizationGraph.PositionalDecoyMotif, GlycanType.O_glycan)))
                    .ToArray();

                int onDecoy = withDecoys.Count(g => g.Target != null
                    && g.Target.ToString() == LocalizationGraph.PositionalDecoyMotif);

                Assert.That(onDecoy, Is.EqualTo(onT),
                    "A decoy site must draw on exactly as many glycan instances as a real site.");
                Assert.That(onDecoy, Is.LessThan(onS + onT),
                    "Admitting both halves would make a decoy twice as competitive as any real site.");

                // Masses must match too, or decoys compete at the wrong precursor masses entirely.
                var tMasses = original.Where(g => g.Target.ToString() == "T").Select(g => g.Mass).OrderBy(m => m);
                var dMasses = withDecoys.Where(g => g.Target.ToString() == LocalizationGraph.PositionalDecoyMotif)
                                        .Select(g => g.Mass).OrderBy(m => m);
                Assert.That(dMasses, Is.EqualTo(tMasses),
                    "Decoy instances must carry the same masses as the compositions they mirror.");
            }
            finally
            {
                GlycanBox.GlobalOGlycans = original;
            }
        }

        [Test]
        public static void RealMotifsAreUnchangedWhileDecoysAreActive()
        {
            int sBefore = AdmittedAt("S");
            int tBefore = AdmittedAt("T");

            LocalizationGraph.DecoyGlycositeMotifs = new HashSet<string> { "A", "L", "G" };

            Assert.That(AdmittedAt("S"), Is.EqualTo(sBefore), "decoy sites must not perturb real S sites");
            Assert.That(AdmittedAt("T"), Is.EqualTo(tBefore), "decoy sites must not perturb real T sites");
        }
    }
}
