#nullable enable
using System;

namespace EngineLayer.DiaLibrarySearch;

/// <param name="FragmentTolerancePpm">How far an observed peak may sit from a library fragment and still be it.</param>
/// <param name="IrtHalfWindow">
/// How far, in iRT, a precursor's apex may fall from its library iRT. Wide enough here to absorb a run whose true
/// RT-to-iRT curve bows away from the linear map it was searched with.
/// </param>
/// <param name="TopFragmentCount">How many of a library entry's most intense fragments are read and scored.</param>
/// <param name="ApexHalfWidthScans">Scans on either side of an apex over which fragment co-elution is measured.</param>
/// <param name="PrecursorSampleStride">
/// Score only precursors whose library index is a multiple of this; 1 scores them all. Targets and decoys are sampled
/// alike, so a sampled search's q-values stay fair. Calibration's first pass uses it to spend less time.
/// </param>
/// <param name="MaxApexCandidates">
/// Candidate peaks scored per precursor. The classifier's score picks among them, so a precursor whose true peak scores
/// second to interference on the apex score can still be found.
/// </param>
/// <param name="ClassifierTrainingQValue">
/// Targets at or below this q-value train the classifier as positives: 0.01, as Percolator. A looser cutoff (0.15)
/// once helped on a library of DIA-NN's own IDs, but on a whole-proteome library its positives are largely false targets and the network ranked hundreds of decoys among the best targets (4,896 against 20,966 precursors at 1%). With interference removal, 0.01 is as good or better on both libraries (results/entrapment/README.md).
/// </param>
/// <param name="InterferenceExplainedFragments">
/// Interference removal (<see cref="DiaInterferenceRemoval"/>): a match is dropped when a better, co-eluting match explains
/// at least this many of its six most intense library fragments. 0 turns removal off.
/// </param>
/// <param name="ClassifierNetworkMembers">
/// Networks in the classifier's ensemble. 12, as DIA-NN, trained 10 epochs each: at a matched entrapment FDP of 1%, +1.6%
/// (PXD005573 1 h) and +3.9% (HF-X) over 5 networks on the whole-proteome library, at about 1.5-2x the search time.
/// </param>
/// <param name="ClassifierNetworkEpochs">Training epochs per network (DIA-NN uses 1).</param>
/// <param name="ClassifierNetworkTrainingSample">
/// Which rows train the network when a fold has more than <see cref="MaxNetworkTrainingRows"/>: a random sample, or the
/// highest-ranked targets and decoys (half each), as DIA-NN trains after removing low-confidence identifications. Confident
/// by default: at a matched entrapment FDP of 1%, +9.4% (PXD022589 HF-X) and +9.9% (PXD005573 1 h) over a random sample,
/// at no extra time. A random sample holds about one real target in twenty on a whole-proteome library.
/// </param>
/// <param name="SiblingTopCandidateOnly">
/// Sibling support (SiblingCoElution, SiblingApexDeltaMinutes) from each other charge state's top candidate peak only,
/// rather than from any of its candidates near this apex, where a decoy gets several chances for noise to line up. On by
/// default: at a matched entrapment FDP of 1%, HF-X 57,200 -> 59,190 (+3.5%), PXD005573 1 h unchanged (39,280 -> 39,170).
/// </param>
/// <param name="MinimumApexFragments">
/// A cheap gate before the full features, as DIA-NN's: a scan can be a candidate apex only if at least this many of the six
/// most intense library fragments are seen there. 0 (default) turns it off.
/// </param>
/// <param name="FragmentApexCandidate">
/// One more candidate apex by a different rule: the apex of the best single fragment's smoothed trace over the whole window,
/// when no candidate already lies within <see cref="ApexHalfWidthScans"/> of it.
/// </param>
/// <param name="ClassifierNormalizationGroups">
/// When set, each classifier fold normalises its scores on a random sample of this many training groups instead of scoring
/// every training row (most of the classifier's scoring time). Null scores every row.
/// </param>
/// <param name="ClassifierNetworkPasses">
/// Network training passes: the second learns from the candidate peaks the first network picked, as DIA-NN trains twice.
/// </param>
/// <param name="ClassifierSeed">
/// The classifier's random seed. 0 is the default; other values give the run-to-run noise a result must exceed.
/// </param>
/// <param name="InterferenceSameMzOnly">
/// As DIA-NN, interference removal drops a match only for a better match with the same precursor m/z (or with it as the +1
/// isotope); false lets any co-eluting match in the same or a neighbouring window explain it.
/// </param>
/// <param name="MaxNetworkTrainingRows">
/// When set, the classifier's network trains on a random subsample of at most this many rows per fold (every row is still
/// scored). Null trains on all. Default 250,000, about DIA-NN's training size: on the whole-proteome library it took the
/// classifier from 10:10 to 7:03 and changed precursors at 1% from 23,200 to 23,721; smaller libraries never reach it.
/// </param>
/// <param name="ClassifierModel">The model the rescorer fits in each fold: a linear discriminant, or a small neural-network ensemble (DIA-NN's approach).</param>
public sealed record DiaLibrarySearchParameters(double FragmentTolerancePpm = 20, double IrtHalfWindow = 20,
    int TopFragmentCount = 12, int ApexHalfWidthScans = 3, int PrecursorSampleStride = 1, int MaxApexCandidates = 3, double ClassifierTrainingQValue = 0.01,
    StatisticalModels.RescoreModel ClassifierModel = StatisticalModels.RescoreModel.NeuralNetworkEnsemble,
    int InterferenceExplainedFragments = 4,
    int? MaxNetworkTrainingRows = 250_000,
    bool InterferenceSameMzOnly = false,
    int ClassifierSeed = 0,
    int ClassifierNetworkMembers = 12,
    int ClassifierNetworkEpochs = StatisticalModels.TargetDecoyRescorer.NetworkEpochs,
    int ClassifierNetworkPasses = 1,
    StatisticalModels.NetworkTrainingSample ClassifierNetworkTrainingSample = StatisticalModels.NetworkTrainingSample.Confident,
    bool SiblingTopCandidateOnly = true,
    int MinimumApexFragments = 0,
    bool FragmentApexCandidate = false,
    int? ClassifierNormalizationGroups = null)
{
    private readonly int _minimumApexFragments = MinimumApexFragments is >= 0 and <= 6 ? MinimumApexFragments
        : throw new ArgumentOutOfRangeException(nameof(MinimumApexFragments), MinimumApexFragments, "The gate counts the six most intense fragments, so 0 to 6.");

    /// <exception cref="ArgumentOutOfRangeException">Outside 0 to 6.</exception>
    public int MinimumApexFragments
    {
        get => _minimumApexFragments;
        init => _minimumApexFragments = value is >= 0 and <= 6 ? value
            : throw new ArgumentOutOfRangeException(nameof(MinimumApexFragments), value, "The gate counts the six most intense fragments, so 0 to 6.");
    }

    private readonly int _classifierNetworkPasses = ClassifierNetworkPasses >= 1 ? ClassifierNetworkPasses
        : throw new ArgumentOutOfRangeException(nameof(ClassifierNetworkPasses), ClassifierNetworkPasses, "The network needs at least one training pass.");

    /// <exception cref="ArgumentOutOfRangeException">Fewer than one pass.</exception>
    public int ClassifierNetworkPasses
    {
        get => _classifierNetworkPasses;
        init => _classifierNetworkPasses = value >= 1 ? value
            : throw new ArgumentOutOfRangeException(nameof(ClassifierNetworkPasses), value, "The network needs at least one training pass.");
    }

    private readonly int _precursorSampleStride = Positive(PrecursorSampleStride);

    /// <exception cref="ArgumentOutOfRangeException">The stride is less than 1.</exception>
    public int PrecursorSampleStride
    {
        get => _precursorSampleStride;
        init => _precursorSampleStride = Positive(value);
    }

    private readonly int _precursorSampleOffset;

    /// <summary>
    /// Which of the stride's samples is scored: precursors whose library index leaves this remainder. 0 by default. A sampled
    /// search with another offset draws a disjoint sample, so calibration can be repeated to see how much its sample matters.
    /// A search at stride 1 scores every precursor whatever the offset.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">The offset is negative or not less than the stride.</exception>
    public int PrecursorSampleOffset
    {
        get => _precursorSampleOffset;
        init => _precursorSampleOffset = value >= 0 && value < PrecursorSampleStride ? value
            : throw new ArgumentOutOfRangeException(nameof(PrecursorSampleOffset), value, "The sample offset must be at least 0 and less than the stride.");
    }

    private static int Positive(int stride) => stride >= 1 ? stride
        : throw new ArgumentOutOfRangeException(nameof(PrecursorSampleStride), stride, "The precursor sample stride must be at least 1.");
}