using EngineLayer;
using EngineLayer.GlycoSearch;
using NUnit.Framework;
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
