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
/// Targets at or below this q-value train the classifier as positives (pyProphet trains at 0.15 then 0.05; our
/// default is 0.01). A looser cutoff teaches it the marginal cases where extra features help.
/// </param>
/// <param name="InterferenceExplainedFragments">
/// Interference removal (<see cref="DiaInterferenceRemoval"/>): a match is dropped when a better, co-eluting match explains
/// at least this many of its six most intense library fragments. 0 turns removal off.
/// </param>
/// <param name="ClassifierModel">The model the rescorer fits in each fold: a linear discriminant, or a small neural-network ensemble (DIA-NN's approach).</param>
public sealed record DiaLibrarySearchParameters(double FragmentTolerancePpm = 20, double IrtHalfWindow = 20,
    int TopFragmentCount = 12, int ApexHalfWidthScans = 3, int PrecursorSampleStride = 1, int MaxApexCandidates = 1, double ClassifierTrainingQValue = 0.15,
    StatisticalModels.RescoreModel ClassifierModel = StatisticalModels.RescoreModel.NeuralNetworkEnsemble,
    int InterferenceExplainedFragments = 4)
{
    private readonly int _precursorSampleStride = Positive(PrecursorSampleStride);

    /// <exception cref="ArgumentOutOfRangeException">The stride is less than 1.</exception>
    public int PrecursorSampleStride
    {
        get => _precursorSampleStride;
        init => _precursorSampleStride = Positive(value);
    }

    private static int Positive(int stride) => stride >= 1 ? stride
        : throw new ArgumentOutOfRangeException(nameof(PrecursorSampleStride), stride, "The precursor sample stride must be at least 1.");
}