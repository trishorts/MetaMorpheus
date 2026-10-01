#nullable enable
using System;
using System.Linq;
using Chromatography.RetentionTimeCalibration;
using MassSpectrometry;
using MzLibUtil;
using Readers.SpectralLibrary;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>A run's retention-time calibration, and the iRT window a search should use with it.</summary>
/// <param name="AnchorCount">Confident first-pass identifications the calibration was fitted on.</param>
public sealed record DiaIrtCalibration(IrtCalibrationModel Model, double IrtHalfWindow, int AnchorCount);

/// <summary>
/// Calibrates a DIA run onto its library's iRT scale from the run alone, with nothing borrowed from another search.
/// A first pass searches with a provisional line (the library's iRT range laid over the run's MS2 time span) and a wide
/// window. Its confident targets become anchors for mzLib's robust <see cref="IrtCalibration"/>. The fitting lives in
/// mzLib; this class only orchestrates the passes.
/// </summary>
public static class DiaIrtSelfCalibration
{
    /// <summary>The most confident first-pass targets used as anchors; more add time, not accuracy.</summary>
    public const int MaximumAnchors = 2000;

    /// <summary>Anchors must pass this target-decoy q-value.</summary>
    public const double AnchorQValue = 0.01;

    /// <summary>
    /// A first pass with fewer anchors than this samples twice as many precursors and searches again, until it has them
    /// or has searched the whole library. A sparse library (a whole proteome, of which a run holds a few percent) needs this.
    /// </summary>
    public const int DesiredAnchors = 500;

    /// <summary>The first pass's iRT half-window, as a fraction of the library's central iRT range.</summary>
    public const double FirstPassWindowFraction = 0.25;

    /// <summary>
    /// Targets the first pass aims to search. It finds several times <see cref="MaximumAnchors"/> at 1% from these, so a
    /// larger library is sampled down to about this size.
    /// </summary>
    public const int FirstPassTargetCount = 10_000;

    /// <summary>The first pass's <see cref="DiaLibrarySearchParameters.PrecursorSampleStride"/> for a library of this many targets.</summary>
    public static int FirstPassStride(int targetCount, int firstPassTargetCount = FirstPassTargetCount) =>
        Math.Max(1, (int)Math.Ceiling((double)targetCount / firstPassTargetCount));

    /// <exception cref="MetaMorpheusException">
    /// The run has no DIA MS2 scans, the library holds no targets, or too few confident first-pass identifications to
    /// calibrate on.
    /// </exception>
    public static DiaIrtCalibration Calibrate(MsDataScan[] scans, MslLibrary library, DiaLibrarySearchParameters parameters,
        CommonParameters commonParameters, IrtCalibrationOptions? options = null, int firstPassTargetCount = FirstPassTargetCount)
    {
        ArgumentNullException.ThrowIfNull(scans);
        ArgumentNullException.ThrowIfNull(library);
        ArgumentNullException.ThrowIfNull(parameters);
        options ??= new IrtCalibrationOptions();

        double[] ms2Rts = scans.Where(s => s.MsnOrder == 2 && s.IsolationRange is not null).Select(s => s.RetentionTime).ToArray();
        if (ms2Rts.Length == 0)
            throw new MetaMorpheusException("Could not calibrate retention time: the run has no DIA MS2 scans.");

        double[] targetIrts = library.QueryMzWindow(0, float.MaxValue).ToArray()
            .Where(e => e.IsDecoy == 0).Select(e => (double)e.Irt).Order().ToArray();
        if (targetIrts.Length == 0)
            throw new MetaMorpheusException("Could not calibrate retention time: the library holds no target precursors.");

        double lowIrt = targetIrts[(int)(0.02 * (targetIrts.Length - 1))];
        double highIrt = targetIrts[(int)(0.98 * (targetIrts.Length - 1))];
        var provisional = IrtCalibration.Line((new RtMinutes(ms2Rts.Min()), new Irt(lowIrt)), (new RtMinutes(ms2Rts.Max()), new Irt(highIrt)));

        int stride = FirstPassStride(targetIrts.Length, firstPassTargetCount);
        while (true)
        {
            var firstPass = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(scans, library, provisional,
                parameters with
                {
                    IrtHalfWindow = FirstPassWindowFraction * (highIrt - lowIrt),
                    PrecursorSampleStride = stride,
                }, commonParameters, [], []).Run();

            var anchors = firstPass.Matches
                .Where(m => !m.IsDecoy && m.QValue <= AnchorQValue)
                .OrderByDescending(m => m.Score)
                .Take(MaximumAnchors)
                .Select(m => (m.ApexRt, m.LibraryIrt))
                .ToList();
            if (anchors.Count >= DesiredAnchors || stride == 1)
                return Fit(anchors, options);
            stride = Math.Max(1, stride / 2);
        }
    }

    private static DiaIrtCalibration Fit(System.Collections.Generic.List<(RtMinutes ApexRt, Irt LibraryIrt)> anchors, IrtCalibrationOptions options)
    {
        if (anchors.Count < options.MinimumAnchors)
            throw new MetaMorpheusException($"Could not calibrate retention time: the first pass found only {anchors.Count} confident " +
                $"identifications, and at least {options.MinimumAnchors} are needed.");

        var model = IrtCalibration.Fit(anchors, options);
        return new DiaIrtCalibration(model, Math.Max(5, 4 * model.ResidualSd), anchors.Count);
    }
}
