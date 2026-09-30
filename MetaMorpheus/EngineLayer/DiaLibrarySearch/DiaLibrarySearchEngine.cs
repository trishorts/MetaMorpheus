#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MassSpectrometry;
using Chromatography.RetentionTimeCalibration;
using FlashLFQ;
using MathNet.Numerics.Statistics;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using StatisticalModels;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// Library-based DIA search. It scores every library precursor isolated by each DIA window whose library iRT falls
/// within reach of the run, then assigns target-decoy q-values.
/// </summary>
/// <remarks>
/// For each isolation window, the window's retention-time span is mapped into iRT, and the library is queried there
/// with decoys included. For each candidate, its fragments are read from every scan in the window whose iRT is within
/// <see cref="DiaLibrarySearchParameters.IrtHalfWindow"/> of the library iRT. The apex is the scan with the most summed
/// fragment intensity, and the score is the cosine between the apex intensities and the library's. A precursor isolated
/// by several windows keeps its best score. q-values are (D+1)/T among targets, from mzLib's
/// <see cref="DeconvolutionQValueCalculator.AssignQValues"/>.
/// <para>
/// Everything is compared in iRT. The library is never rewritten in minutes.
/// </para>
/// </remarks>
public class DiaLibrarySearchEngine : MetaMorpheusEngine
{
    private readonly MsDataScan[] _scans;
    private readonly MslLibrary _library;
    private readonly IrtCalibrationModel _irtMap;
    private readonly DiaLibrarySearchParameters _parameters;

    // The run's MS1 peaks, indexed once with FlashLFQ's indexer; null when the run has no MS1 scans
    private PeakIndexingEngine? _ms1Index;
    private double[] _ms1Rts = [];

    /// <param name="scans">The run's scans. Only MS2 scans with an isolation range are searched.</param>
    /// <param name="library">A loaded library holding targets and decoys. The engine does not dispose it.</param>
    /// <param name="irtMap">This run's calibration from minutes to library iRT (mzLib <see cref="IrtCalibrationModel"/>).</param>
    public DiaLibrarySearchEngine(MsDataScan[] scans, MslLibrary library, IrtCalibrationModel irtMap, DiaLibrarySearchParameters parameters,
        CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
        List<string> nestedIds)
        : base(commonParameters, fileSpecificParameters, nestedIds)
    {
        _scans = scans ?? throw new ArgumentNullException(nameof(scans));
        _library = library ?? throw new ArgumentNullException(nameof(library));
        _irtMap = irtMap ?? throw new ArgumentNullException(nameof(irtMap));
        _parameters = parameters ?? throw new ArgumentNullException(nameof(parameters));
    }

    /// <exception cref="MetaMorpheusException">The library holds no decoys, so no q-value could be estimated.</exception>
    protected override MetaMorpheusEngineResults RunSpecific()
    {
        if (_library.DecoyCount == 0)
            throw new MetaMorpheusException("The spectral library contains no decoys, so a DIA search cannot estimate its FDR. " +
                "Search a library that includes decoy precursors.");

        Status("Running DIA library search...");
        var ms1Scans = _scans.Where(scan => scan.MsnOrder == 1).OrderBy(scan => scan.RetentionTime).ToArray();
        _ms1Index = ms1Scans.Length > 0 ? PeakIndexingEngine.InitializeIndexingEngine(ms1Scans) : null;
        _ms1Rts = ms1Scans.Select(scan => scan.RetentionTime).ToArray();
        var tolerance = new PpmTolerance(_parameters.FragmentTolerancePpm);
        var windows = _scans
            .Where(scan => scan.MsnOrder == 2 && scan.IsolationRange is not null)
            .GroupBy(scan => (scan.IsolationRange.Minimum, scan.IsolationRange.Maximum))
            .OrderBy(window => window.Key.Minimum)
            .ToList();

        var best = new Dictionary<int, DiaPrecursorMatch>();
        for (int w = 0; w < windows.Count; w++)
        {
            if (GlobalVariables.StopLoops)
                break;
            var scans = windows[w].OrderBy(scan => scan.RetentionTime).ToArray();
            foreach (var match in SearchWindow(windows[w].Key, scans, tolerance))
                if (!best.TryGetValue(match.PrecursorIndex, out var existing) || match.Score > existing.Score)
                    best[match.PrecursorIndex] = match;
            ReportProgress(new ProgressEventArgs((int)(100.0 * (w + 1) / windows.Count), "Searching DIA windows...", NestedIds));
        }

        var matches = AssignQValues(best.Values.ToList());
        Status("Done.");
        return new DiaLibrarySearchResults(this, matches);
    }

    private List<DiaPrecursorMatch> SearchWindow((double Minimum, double Maximum) window, MsDataScan[] scans, PpmTolerance tolerance)
    {
        var matches = new List<DiaPrecursorMatch>();
        if (scans.Length == 0)
            return matches;

        double[] scanIrts = scans.Select(scan => _irtMap.ToIrt(new RtMinutes(scan.RetentionTime)).Value).ToArray();
        double irtLow = scanIrts.Min() - _parameters.IrtHalfWindow;
        double irtHigh = scanIrts.Max() + _parameters.IrtHalfWindow;

        MslPrecursorIndexEntry[] candidates;
        using (var hits = _library.QueryWindow((float)window.Minimum, (float)window.Maximum, (float)irtLow, (float)irtHigh, includeDecoys: true))
            candidates = hits.Entries.ToArray();

        foreach (var candidate in candidates)
        {
            if (candidate.PrecursorIdx % _parameters.PrecursorSampleStride != 0)
                continue;
            var entry = _library.GetEntry(candidate.PrecursorIdx);
            if (entry is null || entry.MatchedFragmentIons.Count == 0)
                continue;
            var match = ScoreCandidate(candidate, entry, scans, scanIrts, tolerance);
            if (match is not null)
                matches.Add(match);
        }
        return matches;
    }

    /// <summary>
    /// Reads the entry's most intense fragments across the scans within reach of its library iRT, then scores the peak
    /// group with mzLib's <see cref="FragmentCoElution"/>: the apex is where the fragments co-elute in library
    /// proportions, and the score is the apex cosine times the co-elution around it.
    /// </summary>
    private DiaPrecursorMatch? ScoreCandidate(MslPrecursorIndexEntry candidate, MslLibraryEntry entry, MsDataScan[] scans,
        double[] scanIrts, PpmTolerance tolerance)
    {
        var allIntensities = entry.MatchedFragmentIons.Select(f => (double)f.Intensity).ToArray();
        var fragments = FragmentCoElution.TopIndices(allIntensities, _parameters.TopFragmentCount)
            .Select(i => entry.MatchedFragmentIons[i]).ToList();
        double[] libraryIntensities = fragments.Select(f => (double)f.Intensity).ToArray();

        var reachable = Enumerable.Range(0, scans.Length)
            .Where(s => Math.Abs(scanIrts[s] - candidate.Irt) <= _parameters.IrtHalfWindow)
            .ToArray();
        if (reachable.Length == 0)
            return null;

        var traces = fragments.Select(_ => new double[reachable.Length]).ToArray();
        for (int k = 0; k < reachable.Length; k++)
        {
            var spectrum = scans[reachable[k]].MassSpectrum;
            if (spectrum.Size == 0)
                continue;
            for (int f = 0; f < fragments.Count; f++)
            {
                int i = spectrum.GetClosestPeakIndex(fragments[f].Mz);
                if (tolerance.Within(spectrum.XArray[i], fragments[f].Mz))
                    traces[f][k] = spectrum.YArray[i];
            }
        }

        int apex = FragmentCoElution.FindApex(traces, libraryIntensities, _parameters.ApexHalfWidthScans);
        if (apex < 0)
            return null;

        double[] apexIntensities = traces.Select(trace => trace[apex]).ToArray();
        double coElution = FragmentCoElution.Score(traces,
            Math.Max(0, apex - _parameters.ApexHalfWidthScans), Math.Min(reachable.Length - 1, apex + _parameters.ApexHalfWidthScans));
        double cosine = SpectralSimilarity.CosineOfAlignedVectors(apexIntensities, libraryIntensities);

        // Mass accuracy of the fragments seen at the apex
        var apexSpectrum = scans[reachable[apex]].MassSpectrum;
        var ppmErrors = new List<double>();
        for (int f = 0; f < fragments.Count; f++)
        {
            if (apexIntensities[f] <= 0)
                continue;
            double observed = apexSpectrum.XArray[apexSpectrum.GetClosestPeakIndex(fragments[f].Mz)];
            ppmErrors.Add(Math.Abs(observed - fragments[f].Mz) / fragments[f].Mz * 1e6);
        }

        double[] summed = Enumerable.Range(0, reachable.Length).Select(k => traces.Sum(trace => trace[k])).ToArray();
        var (peakStart, peakEnd) = FragmentCoElution.PeakBounds(summed, apex);
        var (ms1Correlation, ms1PpmError, ms1ApexIntensity) = Ms1Evidence(candidate.PrecursorMz, scans, reachable, summed, apex);

        var apexRt = new RtMinutes(scans[reachable[apex]].RetentionTime);
        var apexIrt = _irtMap.ToIrt(apexRt);

        double[] features =
        [
            cosine,
            coElution,
            (double)ppmErrors.Count / fragments.Count,
            ppmErrors.Count > 0 ? ppmErrors.Average() : _parameters.FragmentTolerancePpm,
            Math.Abs(apexIrt.Value - candidate.Irt),
            Math.Log10(1 + apexIntensities.Sum()),
            peakEnd - peakStart + 1,
            ms1Correlation,
            ms1PpmError,
            Math.Log10(1 + ms1ApexIntensity),
        ];

        return new DiaPrecursorMatch(
            candidate.PrecursorIdx,
            entry.FullSequence,
            candidate.Charge,
            candidate.PrecursorMz,
            candidate.IsDecoy != 0,
            new Irt(candidate.Irt),
            apexRt,
            apexIrt,
            cosine * coElution,
            features);
    }

    /// <summary>
    /// The precursor's own evidence in MS1: its monoisotopic peak in the MS1 scan nearest each reachable MS2 scan. The
    /// correlation is Pearson's between that trace and the summed fragment trace around the apex (0 when either is flat or
    /// the run has no MS1). The ppm error and intensity are taken at the apex, and a missing peak counts as the full
    /// precursor tolerance and zero intensity.
    /// </summary>
    private (double Correlation, double AbsolutePpmError, double ApexIntensity) Ms1Evidence(double precursorMz, MsDataScan[] scans,
        int[] reachable, double[] summedFragments, int apex)
    {
        double tolerancePpm = CommonParameters.PrecursorMassTolerance.Value;
        if (_ms1Index is null)
            return (0, tolerancePpm, 0);

        var trace = new double[reachable.Length];
        IIndexedPeak? apexPeak = null;
        for (int k = 0; k < reachable.Length; k++)
        {
            var peak = _ms1Index.GetIndexedPeak(precursorMz, NearestMs1(scans[reachable[k]].RetentionTime), CommonParameters.PrecursorMassTolerance);
            trace[k] = peak?.Intensity ?? 0;
            if (k == apex)
                apexPeak = peak;
        }

        int from = Math.Max(0, apex - _parameters.ApexHalfWidthScans);
        int length = Math.Min(reachable.Length - 1, apex + _parameters.ApexHalfWidthScans) - from + 1;
        double correlation = length < 3 ? 0 : Correlation.Pearson(trace.Skip(from).Take(length), summedFragments.Skip(from).Take(length));
        double ppm = apexPeak is null ? tolerancePpm : Math.Abs(apexPeak.M - precursorMz) / precursorMz * 1e6;
        return (double.IsFinite(correlation) ? correlation : 0, ppm, apexPeak?.Intensity ?? 0);
    }

    private int NearestMs1(double retentionTime)
    {
        int i = Array.BinarySearch(_ms1Rts, retentionTime);
        if (i >= 0)
            return i;
        i = ~i;
        if (i == 0)
            return 0;
        if (i == _ms1Rts.Length)
            return _ms1Rts.Length - 1;
        return retentionTime - _ms1Rts[i - 1] <= _ms1Rts[i] - retentionTime ? i - 1 : i;
    }

    /// <summary>
    /// Combines each match's features into one score with mzLib's semi-supervised <see cref="TargetDecoyRescorer"/>, with folds
    /// grouped by sequence so that no match is scored by a model trained on it. Then assigns (D+1)/T q-values among targets.
    /// If rescoring cannot run (for example too few matches), the pre-rescoring score stands.
    /// </summary>
    private static List<DiaPrecursorMatch> AssignQValues(List<DiaPrecursorMatch> matches)
    {
        if (matches.Count > 0)
        {
            var rescored = TargetDecoyRescorer.Score(
                matches.Select(m => m.Features).ToList(),
                matches.Select(m => m.IsDecoy).ToList(),
                matches.Select(m => m.FullSequence).ToList());
            if (rescored.Scores.All(double.IsFinite))
                matches = matches.Select((m, i) => m with { Score = rescored.Scores[i] }).ToList();
        }

        var targets = matches.Where(m => !m.IsDecoy).ToList();
        var decoyScores = matches.Where(m => m.IsDecoy).Select(m => m.Score).ToList();
        double[] qValues = DeconvolutionQValueCalculator.AssignQValues(targets.Select(m => m.Score).ToList(), decoyScores);

        var assigned = targets.Select((m, i) => m with { QValue = qValues[i] }).ToList();
        assigned.AddRange(matches.Where(m => m.IsDecoy));
        return assigned;
    }
}

/// <summary>The precursors a <see cref="DiaLibrarySearchEngine"/> scored, targets with their q-values.</summary>
public class DiaLibrarySearchResults(DiaLibrarySearchEngine engine, List<DiaPrecursorMatch> matches) : MetaMorpheusEngineResults(engine)
{
    public List<DiaPrecursorMatch> Matches { get; init; } = matches;

    public int TargetCount => Matches.Count(m => !m.IsDecoy);

    public int DecoyCount => Matches.Count(m => m.IsDecoy);

    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine(base.ToString());
        sb.AppendLine($"Target precursors: {TargetCount}");
        sb.AppendLine($"Decoy precursors: {DecoyCount}");
        sb.AppendLine($"Target precursors with q-value <= 0.01: {Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01)}");
        return sb.ToString();
    }
}
