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
/// every training row (most of the classifier's scoring time). Null scores every row. 300,000 by default: IDs unchanged
/// within 0.2% at a matched paired entrapment FDP of 1%, and the classifier about a sixth faster.
/// </param>
/// <param name="MostIntenseFragmentPeak">
/// Read each fragment from the most intense peak within tolerance (DIA-NN 1.8) rather than the nearest one. On by default:
/// at a matched paired entrapment FDP of 1%, +0.2% (HF-X) and +1.4% (PXD005573 1 h).
/// </param>
/// <param name="DiaNnPeakFinding">
/// Candidate apexes as DIA-NN 1.8 picks them (<see cref="DiaLibrarySearchEngine.DiaNnCandidateApexes"/>) rather than by the
/// apex score, still at most <see cref="MaxApexCandidates"/>. On by default with 6 candidates: neutral at 3, +1.3% (HF-X) and
/// +2.8% (PXD005573 1 h) at 6, at a matched paired entrapment FDP of 1%; with <see cref="DiaNnScores"/>, +1.4% and +3.6%.
/// </param>
/// <param name="DiaNnPeakFindingMs1">With <see cref="DiaNnPeakFinding"/>, add the MS1 trace's correlation with the best fragment to the candidate score, as DIA-NN does.</param>
/// <param name="DiaNnSignalShare">DIA-NN's pSig: the six core fragments' shares of their window signal, as SignalShare1-6.</param>
/// <param name="ClassifierLinearIterations">
/// Rankings of the training rows before the network: the first by the best single feature, each later one by a refit
/// linear discriminant, which also picks each precursor's top candidate peak for the network (DIA-NN refits about 8 times).
/// </param>
/// <param name="FragmentMinimumMz">With <see cref="DiaNnFragmentFilter"/>, the lowest fragment m/z scored (DIA-NN: 200).</param>
/// <param name="FragmentMinimumResidues">With <see cref="DiaNnFragmentFilter"/>, the fewest residues a scored fragment spans (DIA-NN: 3).</param>
/// <param name="FragmentMaximumCharge">With <see cref="DiaNnFragmentFilter"/>, the highest fragment charge scored; 0 (default) sets no limit.</param>
/// <param name="FragmentChargeBelowPrecursor">With <see cref="DiaNnFragmentFilter"/>, score only fragments of lower charge than their precursor (a 2+ precursor's 1+ fragments, a 3+ precursor's 1+ and 2+).</param>
/// <param name="Ms1Offset">The run's MS1 offset over retention time (fitted in calibration); every MS1 lookup is shifted by it. Null shifts nothing.</param>
/// <param name="Ms1PeakFeatures">Features: the M0 mass error over the co-elution window (intensity-weighted) and the M+2 isotope trace's co-elution.</param>
/// <param name="ClassifierNetworkPositiveQValue">With the confident training sample, positives are only targets passing this q-value in the linear ranking (with as many top decoys); null keeps the top half-cap of targets.</param>
/// <param name="ClassifierNetworkLayers">The classifier network's hidden layers, input side first; null keeps the rescorer's default (DIA-NN 2020's 25-20-15-10-5).</param>
/// <param name="SqrtCoElution">A feature: the core fragments' co-elution on square-root traces, which damp the apex so the peak's flanks weigh more.</param>
/// <param name="Ms1ToleranceSpreadMultiple">When set, the MS1 tolerance follows the run's own mass accuracy: this multiple of the calibrated offset's residual spread (<see cref="EffectiveMs1TolerancePpm"/>).</param>
/// <param name="Ms1PeakEnvelope">A feature: the MS1 M-1 to M3 envelope summed over the co-elution window's MS1 scans, against the expected pattern with nothing at M-1.</param>
/// <param name="Ms1TightCorrelation">A feature: MS1 M0 co-elution from peaks within half the MS1 tolerance, as DIA-NN scores MS1 at several tolerances.</param>
/// <param name="Ms1TolerancePpm">
/// MS1 tolerance for the precursor's traces, isotope envelope and mass error, around the offset-corrected m/z
/// (<see cref="Ms1Offset"/>). 5 ppm by default, with the calibrated offset applied: against 10 ppm without it, two-seed means
/// at a matched paired entrapment FDP of 1% gain 1.6% (HF-X) and 1.3% (PXD005573 1 h). 10 ppm against 20 gained 1.1% and
/// 1.6%. Without the offset, 5 ppm clipped real precursors as the run's error drifted up to 3 ppm off zero.
/// </param>
/// <param name="DiaNnFragmentFilter">
/// Score only fragments DIA-NN would (at least 3 residues, 200-1800 m/z), when at least 3 of a precursor's qualify. On by
/// default: two-seed means at a matched paired entrapment FDP of 1%, HF-X +0.8%, PXD005573 1 h +3.9%.
/// </param>
/// <param name="MaxToleranceCoElution">
/// A feature as DIA-NN scores co-elution: each core fragment's best correlation with the profile among its traces at the full
/// tolerance, 0.45x and 0.2x (peaks beyond the fraction dropped), averaged.
/// </param>
/// <param name="ExtraFragmentCount">
/// Library fragments after the top <see cref="TopFragmentCount"/> (by library intensity) read as their own features
/// (ExtraCoElution, ExtraMatchedFraction, ExtraWeightedCoElution), as DIA-NN adds its remaining fragments: averaged into the
/// core scores they dilute them, since most are faint. 12 by default (fragments 13-24): +0.9% HF-X, +0.5% PXD005573 1 h
/// (two-seed means, matched paired entrapment FDP 1%). 0 reads none.
/// </param>
/// <param name="DiaNnScores">
/// DIA-NN 1.8 scores we otherwise lack: MinCorr, NFCorr and ShadowCorr (see <see cref="DiaPrecursorMatch.FeatureNames"/>).
/// Off, they are 0. On (default), seven more traces are read per precursor: +1.0% (HF-X) and +1.4% (PXD005573 1 h) at a
/// matched paired entrapment FDP of 1%.
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
    int TopFragmentCount = 12, int ApexHalfWidthScans = 3, int PrecursorSampleStride = 1, int MaxApexCandidates = 10, double ClassifierTrainingQValue = 0.01,
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
    int? ClassifierNormalizationGroups = 300_000,
    bool MostIntenseFragmentPeak = true,
    bool DiaNnPeakFinding = true,
    bool DiaNnScores = true,
    bool DiaNnPeakFindingMs1 = false,
    bool DiaNnSignalShare = false,
    int ClassifierLinearIterations = 3,
    int ExtraFragmentCount = 12,
    bool MaxToleranceCoElution = false,
    bool DiaNnFragmentFilter = true,
    double FragmentMinimumMz = 200,
    int FragmentMinimumResidues = 3,
    int FragmentMaximumCharge = 0,
    bool FragmentChargeBelowPrecursor = false,
    double Ms1TolerancePpm = 5,
    bool Ms1TightCorrelation = false,
    int[]? ClassifierNetworkLayers = null,
    double? ClassifierNetworkPositiveQValue = null,
    bool Ms1PeakFeatures = false,
    Ms1OffsetModel? Ms1Offset = null,
    bool SqrtCoElution = false,
    bool Ms1PeakEnvelope = true,
    double? Ms1ToleranceSpreadMultiple = null)
{
    /// <summary>
    /// The MS1 tolerance a search uses: <see cref="Ms1ToleranceSpreadMultiple"/> times the calibrated offset's residual
    /// spread (<see cref="Ms1OffsetModel.ResidualSpreadPpm"/>) when both are known, else <see cref="Ms1TolerancePpm"/>.
    /// </summary>
    public double EffectiveMs1TolerancePpm =>
        Ms1ToleranceSpreadMultiple is { } multiple && Ms1Offset is { ResidualSpreadPpm: var spread } && double.IsFinite(spread) && spread > 0
            ? multiple * spread
            : Ms1TolerancePpm;

    private readonly int _minimumApexFragments = MinimumApexFragments is >= 0 and <= 6 ? MinimumApexFragments
        : throw new ArgumentOutOfRangeException(nameof(MinimumApexFragments), MinimumApexFragments, "The gate counts the six most intense fragments, so 0 to 6.");

    /// <exception cref="ArgumentOutOfRangeException">Outside 0 to 6.</exception>
    public int MinimumApexFragments
    {
        get => _minimumApexFragments;
        init => _minimumApexFragments = value is >= 0 and <= 6 ? value
            : throw new ArgumentOutOfRangeException(nameof(MinimumApexFragments), value, "The gate counts the six most intense fragments, so 0 to 6.");
    }

    private readonly int _extraFragmentCount = ExtraFragmentCount >= 0 ? ExtraFragmentCount
        : throw new ArgumentOutOfRangeException(nameof(ExtraFragmentCount), ExtraFragmentCount, "The extra fragment count cannot be negative.");

    /// <exception cref="ArgumentOutOfRangeException">Negative.</exception>
    public int ExtraFragmentCount
    {
        get => _extraFragmentCount;
        init => _extraFragmentCount = value >= 0 ? value
            : throw new ArgumentOutOfRangeException(nameof(ExtraFragmentCount), value, "The extra fragment count cannot be negative.");
    }

    private readonly int _classifierLinearIterations = ClassifierLinearIterations >= 1 ? ClassifierLinearIterations
        : throw new ArgumentOutOfRangeException(nameof(ClassifierLinearIterations), ClassifierLinearIterations, "The linear model needs at least one iteration.");

    /// <exception cref="ArgumentOutOfRangeException">Fewer than one iteration.</exception>
    public int ClassifierLinearIterations
    {
        get => _classifierLinearIterations;
        init => _classifierLinearIterations = value >= 1 ? value
            : throw new ArgumentOutOfRangeException(nameof(ClassifierLinearIterations), value, "The linear model needs at least one iteration.");
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