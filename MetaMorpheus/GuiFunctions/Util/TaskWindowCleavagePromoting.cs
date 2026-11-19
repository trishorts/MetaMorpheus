using System.Linq;
using Omics.Digestion;

namespace GuiFunctions.Util;

/// <summary>
/// When a task window recommends ticking "Respect cleavage-promoting modifications": the user picked a protease that cuts
/// only beside a modification (StcE needs an O-glycan at P2), the box is not ticked, and the search is fully specific.
/// </summary>
/// <remarks>
/// Without the box such a protease cuts wherever its sequence motif fits, as if the glycan were always there. mzLib applies
/// the setting to fully specific digestion only (<c>CleavageBlockingPolicy.For</c>), so the window does not ask otherwise.
/// </remarks>
public static class TaskWindowCleavagePromoting
{
    public const string RecommendationTitle = "Use O-Glycoprotease Recommendations?";

    public const string RecommendationMessage =
        "We recommend using the following parameters for O-glycoprotease searches:\n" +
        "\t-Check 'Respect cleavage-promoting modifications'\n" +
        "\t-Add an O-glycan as a variable modification, or use a database annotated with O-glycosites:\n" +
        "\t with the box checked, these enzymes cut only where an O-glycan can sit\n\n" +
        "Would you like to use these recommended settings?";

    /// <summary>True for a protease with a motif that cuts only beside a modification, the test mzLib's digestion uses.</summary>
    public static bool NeedsAPromotingModification(DigestionAgent protease) =>
        protease is not null && protease.DigestionMotifs.Any(motif => motif.CleavageRequirements.Count > 0);

    /// <summary>Whether the window should recommend ticking the box for the protease the user just picked.</summary>
    public static bool RecommendTickingTheBox(DigestionAgent protease, bool boxTicked, CleavageSpecificity searchMode) =>
        !boxTicked && searchMode == CleavageSpecificity.Full && NeedsAPromotingModification(protease);
}
