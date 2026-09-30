#nullable enable
using MzLibUtil;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// One library precursor scored against one DIA run: its apex, where that falls in iRT, and its score. Deliberately not a
/// <see cref="SpectralMatch"/>: a DIA identification is a chromatographic peak group, not one scan's hypothesis.
/// </summary>
/// <param name="PrecursorIndex">Index of the precursor in the library searched.</param>
/// <param name="Score">
/// The combined score the q-values are computed on: mzLib's semi-supervised TargetDecoyRescorer over
/// <paramref name="Features"/>. Before rescoring it holds the apex cosine × co-elution.
/// </param>
/// <param name="Features">The peak group's features, in <see cref="DiaPrecursorMatch.FeatureNames"/> order.</param>
/// <param name="QValue">
/// Target-decoy q-value among targets; NaN for a decoy, which competes but is not itself reported.
/// </param>
public sealed record DiaPrecursorMatch(
    int PrecursorIndex,
    string FullSequence,
    int Charge,
    double PrecursorMz,
    bool IsDecoy,
    Irt LibraryIrt,
    RtMinutes ApexRt,
    Irt ApexIrt,
    double Score,
    double[] Features,
    double QValue = double.NaN)
{
    /// <summary>What each entry of <see cref="Features"/> measures, in order.</summary>
    public static readonly string[] FeatureNames =
    [
        "ApexCosine",           // library-proportion match at the apex
        "CoElution",            // fragments rising and falling together around the apex
        "MatchedFraction",      // share of scored fragments seen at the apex
        "MeanAbsolutePpmError", // mass accuracy of the fragments seen at the apex
        "AbsoluteDeltaIrt",     // |apex iRT − library iRT| under this run's calibration
        "LogApexIntensity",     // log10(1 + summed fragment intensity at the apex)
        "PeakWidthScans",       // width of the summed fragment peak around the apex
        "FragmentCorrelation1", // each core fragment's correlation with the best-fragment profile, by library rank
        "FragmentCorrelation2",
        "FragmentCorrelation3",
        "FragmentCorrelation4",
        "FragmentCorrelation5",
        "FragmentCorrelation6",
        "TightCoElution",       // co-elution counting only peaks within 0.45 x the fragment tolerance
        "RemainingCoElution",   // correlation of the non-core fragments with the profile
        "WindowCosine",         // library cosine across the peak, weighted by the profile squared
    ];
}
