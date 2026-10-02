#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Chromatography.RetentionTimeCalibration;
using MassSpectrometry;
using MzLibUtil;
using Readers.SpectralLibrary;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>A run's retention-time calibration, and the iRT window a search should use with it.</summary>
/// <param name="AnchorCount">Confident first-pass identifications the calibration was fitted on.</param>
/// <param name="Rounds">Calibration rounds run: the first with the provisional line, each later one with the previous fit.</param>
public sealed record DiaIrtCalibration(IrtCalibrationModel Model, double IrtHalfWindow, int AnchorCount, int Rounds = 1);

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

    /// <summary>The search window's half-width, in residual SDs of the calibration (at least 5 iRT).</summary>
    public const double DefaultWindowSds = 4;

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
        CommonParameters commonParameters, IrtCalibrationOptions? options = null, int firstPassTargetCount = FirstPassTargetCount,
        int desiredAnchors = DesiredAnchors, int rounds = 1, double windowSds = DefaultWindowSds)
    {
        ArgumentNullException.ThrowIfNull(scans);
        if (rounds < 1)
            throw new ArgumentOutOfRangeException(nameof(rounds), rounds, "At least one calibration round is needed.");
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
            if (anchors.Count >= Math.Min(desiredAnchors, MaximumAnchors) || stride == 1)
                return LaterRounds(Fit(anchors, options, windowSds), scans, library, parameters, commonParameters, options, stride, rounds, windowSds);
            stride = Math.Max(1, stride / 2);
        }
    }

    /// <summary>
    /// Iterative calibration, as DIA-NN calibrates in rounds: the same sample searched again with the previous round's fit and
    /// window, so anchors come from the whole gradient rather than where the provisional straight line happened to be right.
    /// A round that cannot fit keeps the previous calibration.
    /// </summary>
    private static DiaIrtCalibration LaterRounds(DiaIrtCalibration calibration, MsDataScan[] scans, MslLibrary library,
        DiaLibrarySearchParameters parameters, CommonParameters commonParameters, IrtCalibrationOptions options, int stride, int rounds, double windowSds)
    {
        for (int round = 2; round <= rounds; round++)
        {
            var pass = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(scans, library, calibration.Model,
                parameters with { IrtHalfWindow = calibration.IrtHalfWindow, PrecursorSampleStride = stride }, commonParameters, [], []).Run();
            var anchors = pass.Matches
                .Where(m => !m.IsDecoy && m.QValue <= AnchorQValue)
                .OrderByDescending(m => m.Score)
                .Take(MaximumAnchors)
                .Select(m => (m.ApexRt, m.LibraryIrt))
                .ToList();
            if (anchors.Count < options.MinimumAnchors)
                break;
            calibration = Fit(anchors, options, windowSds) with { Rounds = round };
        }
        return calibration;
    }

    /// <summary>The most confident main-search targets a <see cref="Refine"/> fit uses.</summary>
    public const int MaximumRefinementAnchors = 20_000;

    /// <summary>
    /// The second pass, as DIA-NN refits after its first search: the run's RT->iRT map refitted on the main search's targets at
    /// <see cref="AnchorQValue"/> (the best <see cref="MaximumRefinementAnchors"/> by score), many more than the first pass's.
    /// Decoys and targets above the cutoff never take part.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="matches"/> is null.</exception>
    /// <exception cref="MetaMorpheusException">Too few confident targets to fit.</exception>
    public static DiaIrtCalibration Refine(IReadOnlyList<DiaPrecursorMatch> matches, IrtCalibrationOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(matches);
        var anchors = matches
            .Where(m => !m.IsDecoy && m.QValue <= AnchorQValue)
            .OrderByDescending(m => m.Score).ThenBy(m => m.PrecursorIndex)
            .Take(MaximumRefinementAnchors)
            .Select(m => (m.ApexRt, m.LibraryIrt))
            .ToList();
        return Fit(anchors, options ?? new IrtCalibrationOptions(), DefaultWindowSds);
    }

    private static DiaIrtCalibration Fit(System.Collections.Generic.List<(RtMinutes ApexRt, Irt LibraryIrt)> anchors, IrtCalibrationOptions options, double windowSds)
    {
        if (anchors.Count < options.MinimumAnchors)
            throw new MetaMorpheusException($"Could not calibrate retention time: the first pass found only {anchors.Count} confident " +
                $"identifications, and at least {options.MinimumAnchors} are needed.");

        var model = IrtCalibration.Fit(anchors, options);
        return new DiaIrtCalibration(model, Math.Max(5, windowSds * model.ResidualSd), anchors.Count);
    }
}
