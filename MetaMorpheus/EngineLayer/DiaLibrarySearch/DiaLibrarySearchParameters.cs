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
public sealed record DiaLibrarySearchParameters(double FragmentTolerancePpm = 20, double IrtHalfWindow = 20,
    int TopFragmentCount = 12, int ApexHalfWidthScans = 3, int PrecursorSampleStride = 1)
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