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
/// <param name="Quantity">
/// The fragments' summed peak areas (trapezoid over the peak, less each trace's baseline): the precursor's amount for
/// quantification. NaN when not measured.
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
    double QValue = double.NaN,
    double Quantity = double.NaN)
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
        "MedianCoElution",      // core fragments' correlation with their per-scan median profile
        "GoodFragments",        // core fragments with r >= 0.75 to the median profile
        "VeryGoodFragment",     // 1 when any core fragment has r >= 0.9
        "GaussianFit",          // median profile's correlation with a Gaussian at the apex
        "AreaSqrtCosine",       // cosine of square-root peak areas with square-root library intensities
        "AreaPearson",          // Pearson of peak areas with library intensities
        "AreaManhattan",        // Manhattan distance of sum-normalised areas and library intensities
        "ApexScoreDelta",       // (this apex's score - best competing scan outside its window) / this score
        "ApexScoreZ",           // z-score of this apex's score among the window's scored scans
        "PeakSignalFraction",   // share of the window's summed fragment signal inside this peak
        "WeightedCoElution",    // fragments' correlation with the profile, weighted by library intensity
        "WeightedMatchedFraction", // library-intensity-weighted share of fragments seen at the apex
        "Top3CoElution",        // correlation with the profile of the three most intense library fragments
        "Top3Cosine",           // apex cosine on the three most intense library fragments
        "Top1Present",          // 1 when the most intense library fragment is seen at the apex
        "Top2Present",          // 1 when the second most intense is
        "Top1PpmError",         // mass error of the most intense fragment at the apex (the tolerance when absent)
        "DetectableMatchedFraction", // share seen at the apex among fragments expected above 3x the scan's noise
        "DetectableCoElution",  // their correlation with the profile
        "DetectableFragments",  // how many fragments were expected to be detectable
        "Ms1Correlation",       // precursor MS1 monoisotopic trace vs the fragment profile (DIA-NN's use of MS1)
        "Ms1IsotopeCorrelation", // precursor M+1 MS1 trace vs the fragment profile
        "FragmentIsotopeOverlap", // weighted share of matched fragments with a bigger peak one isotope below (OpenSWATH)
        "FragmentM1Fraction",   // weighted share of matched fragments that show their own M+1 peak
        "MatchSurprisal",       // sum of -log10 chance-match probability over matched fragments (local peak density x window)
        "SurprisalFraction",    // matched surprisal over the total possible
        "SpecificityWeightedCoElution", // profile correlations weighted by each fragment's surprisal
        "PeptideLength",        // residues: a chance match is harder with more of them; a reversed decoy keeps its target's
        "Ms1EnvelopeCosine",    // M0-M3 at the apex's MS1 scan against the expected isotope pattern
        "Ms1MassErrorPpm",      // |M0 - library precursor m/z| in ppm; the full tolerance when no M0 is found
        "Ms1ApexShare",         // MS1 M0 at the apex (+-1 scan) over its maximum in the window: is the apex on an MS1 peak
        "SiblingCoElution",     // best co-elution of another charge state of the sequence with an apex within the co-elution window
        "SiblingApexDeltaMinutes", // RT gap to the nearest other-charge apex (1 when there is none, capped at 1)
        "PrecursorMz",          // library context (DIA-NN): lets the classifier learn when to trust a prediction; a decoy shares its target's
        "PrecursorCharge",
        "LibraryFragmentCount",
        "ExtraFragmentCoElution",       // fragments beyond the scored top N: mean correlation with the best-fragment profile (0 when none are read)
        "ExtraFragmentMatchedFraction", // their share seen at the apex
    ];
}
