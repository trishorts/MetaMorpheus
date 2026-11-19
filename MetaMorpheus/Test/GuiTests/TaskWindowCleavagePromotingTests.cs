using GuiFunctions.Util;
using NUnit.Framework;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using Transcriptomics.Digestion;

namespace Test.GuiTests
{
    /// <summary>
    /// When a task window recommends ticking "Respect cleavage-promoting modifications" after the user picks a protease.
    /// </summary>
    [TestFixture]
    public static class TaskWindowCleavagePromotingTests
    {
        private static IEnumerable<TestCaseData> EveryProtease() =>
            ProteaseDictionary.Dictionary.Values.Select(protease =>
                new TestCaseData(protease).SetArgDisplayNames(protease.Name.Replace('.', '_')));

        /// <summary>
        /// Every protease the picker lists: the box is recommended exactly for those with a motif that cuts only beside a
        /// modification. The windows ask mzLib's motifs, not a list of names, so a glycoprotease a later mzLib adds is covered.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(EveryProtease))]
        public static void TheBoxIsRecommendedExactlyForAProteaseWithACleavageRequirement(Protease protease)
        {
            bool hasARequirement = protease.DigestionMotifs.Any(motif => motif.CleavageRequirements.Count > 0);
            Assert.That(TaskWindowCleavagePromoting.RecommendTickingTheBox(protease, false, CleavageSpecificity.Full),
                Is.EqualTo(hasARequirement));
        }

        /// <summary>StcE alone and in a co-digest, whose trypsin motifs carry no requirement.</summary>
        [Test]
        [TestCase("StcE")]
        [TestCase("StcE-trypsin")]
        public static void StcE_AndItsCoDigest_AreRecommendedTheBox(string protease) =>
            Assert.That(TaskWindowCleavagePromoting.RecommendTickingTheBox(ProteaseDictionary.Dictionary[protease], false,
                CleavageSpecificity.Full), Is.True);

        [Test]
        [TestCase("trypsin")]
        [TestCase("top-down")]
        [TestCase("non-specific")]
        public static void AProteaseThatCutsWhateverTheResidueCarries_IsNotRecommendedTheBox(string protease) =>
            Assert.That(TaskWindowCleavagePromoting.RecommendTickingTheBox(ProteaseDictionary.Dictionary[protease], false,
                CleavageSpecificity.Full), Is.False);

        [Test]
        public static void AnRnase_OrNoProtease_IsNotRecommendedTheBox()
        {
            Assert.That(TaskWindowCleavagePromoting.RecommendTickingTheBox(RnaseDictionary.Dictionary["RNase T1"], false,
                CleavageSpecificity.Full), Is.False);
            Assert.That(TaskWindowCleavagePromoting.RecommendTickingTheBox(null, false, CleavageSpecificity.Full), Is.False);
        }

        [Test]
        public static void AnAlreadyTickedBox_IsNotRecommendedAgain() =>
            Assert.That(TaskWindowCleavagePromoting.RecommendTickingTheBox(ProteaseDictionary.Dictionary["StcE"], true,
                CleavageSpecificity.Full), Is.False);

        /// <summary>mzLib applies the setting to fully specific digestion only, so a semi- or non-specific search is not asked.</summary>
        [Test]
        [TestCase(CleavageSpecificity.Semi)]
        [TestCase(CleavageSpecificity.None)]
        public static void ASearchThatIsNotFullySpecific_IsNotRecommendedTheBox(CleavageSpecificity searchMode) =>
            Assert.That(TaskWindowCleavagePromoting.RecommendTickingTheBox(ProteaseDictionary.Dictionary["StcE"], false, searchMode),
                Is.False);

        /// <summary>The prompt names the box by its label, and says what ticking it requires of the search.</summary>
        [Test]
        public static void ThePromptNamesTheBox_AndAsksForAnOGlycan()
        {
            Assert.That(TaskWindowCleavagePromoting.RecommendationTitle, Is.EqualTo("Use O-Glycoprotease Recommendations?"));
            Assert.That(TaskWindowCleavagePromoting.RecommendationMessage, Does.StartWith(
                "We recommend using the following parameters for O-glycoprotease searches:\n" +
                "\t-Check 'Respect cleavage-promoting modifications'\n"));
            Assert.That(TaskWindowCleavagePromoting.RecommendationMessage, Does.Contain("O-glycan as a variable modification"));
            Assert.That(TaskWindowCleavagePromoting.RecommendationMessage, Does.EndWith("Would you like to use these recommended settings?"));
        }
    }
}
