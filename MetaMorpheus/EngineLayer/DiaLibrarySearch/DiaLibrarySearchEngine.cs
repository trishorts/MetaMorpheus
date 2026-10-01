#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MassSpectrometry;
using Chromatography.RetentionTimeCalibration;
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

    /// <summary>The most intense library fragments that define the elution profile.</summary>
    private const int CoreFragmentCount = 6;

    /// <summary>Tight co-elution counts only peaks within this fraction of the fragment tolerance (DIA-NN uses 0.45).</summary>
    private const double TightToleranceFraction = 0.45;

    /// <summary>MS1 tolerance for the precursor's elution-profile features (DIA-NN chose 17–22 ppm on PXD005573).</summary>
    private static readonly PpmTolerance Ms1Tolerance = new(20);

    // The run's MS1 peaks, indexed once with FlashLFQ's indexer; null when the run has no MS1 scans
    private FlashLFQ.PeakIndexingEngine? _ms1Index;
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
        var tolerance = new PpmTolerance(_parameters.FragmentTolerancePpm);
        var ms1Scans = _scans.Where(scan => scan.MsnOrder == 1).OrderBy(scan => scan.RetentionTime).ToArray();
        _ms1Index = ms1Scans.Length > 0 ? FlashLFQ.PeakIndexingEngine.InitializeIndexingEngine(ms1Scans) : null;
        _ms1Rts = ms1Scans.Select(scan => scan.RetentionTime).ToArray();
        var windows = _scans
            .Where(scan => scan.MsnOrder == 2 && scan.IsolationRange is not null)
            .GroupBy(scan => (scan.IsolationRange.Minimum, scan.IsolationRange.Maximum))
            .OrderBy(window => window.Key.Minimum)
            .ToList();

        var rows = new List<DiaPrecursorMatch>();
        for (int w = 0; w < windows.Count; w++)
        {
            if (GlobalVariables.StopLoops)
                break;
            var scans = windows[w].OrderBy(scan => scan.RetentionTime).ToArray();
            rows.AddRange(SearchWindow(windows[w].Key, scans, tolerance));
            ReportProgress(new ProgressEventArgs((int)(100.0 * (w + 1) / windows.Count), "Searching DIA windows...", NestedIds));
        }

        var matches = AssignQValues(rows);
        Status("Done.");
        return new DiaLibrarySearchResults(this, matches, DiaPeptideFdr.Assign(matches));
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
            matches.AddRange(ScoreCandidates(candidate, entry, scans, scanIrts, tolerance));
        }
        return matches;
    }

    /// <summary>
    /// Reads the entry's most intense fragments across the scans within reach of its library iRT, then scores the peak
    /// group with mzLib's <see cref="FragmentCoElution"/>: the apex is where the fragments co-elute in library
    /// proportions, and the score is the apex cosine times the co-elution around it. Up to
    /// <see cref="DiaLibrarySearchParameters.MaxApexCandidates"/> candidate apexes are scored, each as its own match; the
    /// classifier's score decides among them (<see cref="AssignQValues"/>).
    /// </summary>
    private IEnumerable<DiaPrecursorMatch> ScoreCandidates(MslPrecursorIndexEntry candidate, MslLibraryEntry entry, MsDataScan[] scans,
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
            yield break;

        var traces = fragments.Select(_ => new double[reachable.Length]).ToArray();
        var ppm = fragments.Select(_ => new double[reachable.Length]).ToArray();
        for (int k = 0; k < reachable.Length; k++)
        {
            var spectrum = scans[reachable[k]].MassSpectrum;
            if (spectrum.Size == 0)
                continue;
            for (int f = 0; f < fragments.Count; f++)
            {
                int i = spectrum.GetClosestPeakIndex(fragments[f].Mz);
                if (tolerance.Within(spectrum.XArray[i], fragments[f].Mz))
                {
                    traces[f][k] = spectrum.YArray[i];
                    ppm[f][k] = Math.Abs(spectrum.XArray[i] - fragments[f].Mz) / fragments[f].Mz * 1e6;
                }
            }
        }

        // The core: the most intense library fragments, in library rank; the rest enter as their own feature
        int[] coreIndices = FragmentCoElution.TopIndices(libraryIntensities, CoreFragmentCount)
            .OrderByDescending(i => libraryIntensities[i]).ThenBy(i => i).ToArray();
        var core = coreIndices.Select(i => traces[i]).ToList();
        var rest = Enumerable.Range(0, fragments.Count).Except(coreIndices).Select(i => traces[i]).ToList();
        // The core again, counting only peaks within a tight fraction of the tolerance
        double tightPpm = TightToleranceFraction * _parameters.FragmentTolerancePpm;
        var tightCore = coreIndices.Select(i => traces[i].Select((v, k) => ppm[i][k] <= tightPpm ? v : 0).ToArray()).ToList();
        double[] summed = Enumerable.Range(0, reachable.Length).Select(k => traces.Sum(trace => trace[k])).ToArray();
        // Every scan's apex score, to judge how far a candidate stands out from the rest of its window (PECAN, OpenSWATH)
        double[] apexScores = FragmentCoElution.ApexScores(traces, libraryIntensities, _parameters.ApexHalfWidthScans);
        double[] scored = apexScores.Where(v => v > 0).ToArray();
        double scoreMean = scored.Length > 0 ? scored.Average() : 0;
        double scoreSd = scored.Length > 1 ? Math.Sqrt(scored.Sum(v => (v - scoreMean) * (v - scoreMean)) / (scored.Length - 1)) : 0;
        double windowSignal = summed.Sum();

        // The precursor's MS1 traces, monoisotopic and M+1, in the MS1 scan nearest each reachable MS2 scan
        double[]? ms1Mono = null, ms1Isotope = null;
        if (_ms1Index is not null)
        {
            ms1Mono = new double[reachable.Length];
            ms1Isotope = new double[reachable.Length];
            double isotopeMz = candidate.PrecursorMz + 1.0033548 / Math.Max(1, (int)candidate.Charge);
            for (int k = 0; k < reachable.Length; k++)
            {
                int ms1 = NearestMs1(scans[reachable[k]].RetentionTime);
                ms1Mono[k] = _ms1Index.GetIndexedPeak(candidate.PrecursorMz, ms1, Ms1Tolerance)?.Intensity ?? 0;
                ms1Isotope[k] = _ms1Index.GetIndexedPeak(isotopeMz, ms1, Ms1Tolerance)?.Intensity ?? 0;
            }
        }
        foreach (int apex in FragmentCoElution.FindApexes(apexScores, _parameters.ApexHalfWidthScans, _parameters.MaxApexCandidates))
        {
            double[] apexIntensities = traces.Select(trace => trace[apex]).ToArray();
            // Co-elution against the best of the six most intense library fragments, smoothed: one reliable profile rather
            // than an average that an interfered fragment drags along (DIA-NN's approach, Demichev et al. 2020)
            int from = Math.Max(0, apex - _parameters.ApexHalfWidthScans);
            int to = Math.Min(reachable.Length - 1, apex + _parameters.ApexHalfWidthScans);
            double coElution = 0, tightCoElution = 0, remainingCoElution = 0;
            var fragmentCorrelations = new double[CoreFragmentCount];
            double[]? reference = null;
            if (to > from)
            {
                reference = FragmentCoElution.Smooth(core[FragmentCoElution.BestFragment(core, from, to)]);
                double[] correlations = FragmentCoElution.CorrelationsTo(core, reference, from, to);
                coElution = correlations.Average();
                Array.Copy(correlations, fragmentCorrelations, correlations.Length);
                if (rest.Count > 0)
                    remainingCoElution = FragmentCoElution.CorrelationsTo(rest, reference, from, to).Average();
                var tightReference = FragmentCoElution.Smooth(tightCore[FragmentCoElution.BestFragment(tightCore, from, to)]);
                tightCoElution = FragmentCoElution.CorrelationsTo(tightCore, tightReference, from, to).Average();
            }
            double cosine = SpectralSimilarity.CosineOfAlignedVectors(apexIntensities, libraryIntensities);

            // Library similarity across the peak, each scan weighted by the profile squared
            double windowCosine = cosine;
            if (reference is not null)
            {
                double weighted = 0, weights = 0;
                for (int s = from; s <= to; s++)
                {
                    double weight = reference[s] * reference[s];
                    if (weight <= 0)
                        continue;
                    weighted += weight * SpectralSimilarity.CosineOfAlignedVectors(traces.Select(trace => trace[s]).ToArray(), libraryIntensities);
                    weights += weight;
                }
                if (weights > 0)
                    windowCosine = weighted / weights;
            }

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

            var (peakStart, peakEnd) = FragmentCoElution.PeakBounds(summed, apex);

            // Median-profile co-elution (EncyclopeDIA, AlphaDIA): each core trace normalised to sum 1 over the window,
            // the per-scan median as a reference that no single fragment can hijack
            double medianCoElution = 0, goodFragments = 0, veryGoodFragment = 0, gaussianFit = 0;
            if (to > from)
            {
                double[] median = MedianProfile(core, from, to, reachable.Length);
                double[] toMedian = FragmentCoElution.CorrelationsTo(core, median, from, to);
                medianCoElution = toMedian.Average();
                goodFragments = toMedian.Count(r => r >= 0.75);
                veryGoodFragment = toMedian.Any(r => r >= 0.9) ? 1 : 0;
                double sigma = Math.Max(1, (to - from) / 4.0);
                double[] gaussian = Enumerable.Range(0, reachable.Length).Select(s => Math.Exp(-0.5 * Math.Pow((s - apex) / sigma, 2))).ToArray();
                gaussianFit = FragmentCoElution.CorrelationsTo([median], gaussian, from, to)[0];
            }

            // Library agreement on peak areas rather than one apex scan (OpenSWATH, Skyline): cosine of square-root areas,
            // Pearson of areas, and Manhattan distance of sum-normalised areas
            double[] areas = traces.Select(trace => Area(trace, peakStart, peakEnd)).ToArray();
            double areaSqrtCosine = SpectralSimilarity.CosineOfAlignedVectors(areas.Select(Math.Sqrt).ToArray(), libraryIntensities.Select(Math.Sqrt).ToArray());
            double areaPearson = areas.Any(a => a > 0) ? Math.Max(0, FragmentCoElution.CorrelationsTo([areas], libraryIntensities, 0, areas.Length - 1)[0]) : 0;
            double areaSum = areas.Sum(), librarySum = libraryIntensities.Sum();
            double areaManhattan = areaSum > 0 ? areas.Select((a, f) => Math.Abs(a / areaSum - libraryIntensities[f] / librarySum)).Sum() : 2;

            // Fewer fragments, weighted by expectation: a weak precursor shows only its strongest fragments above noise, so
            // score it on those rather than penalise it for fragments it could not show
            double[] allCorrelations = reference is null ? new double[fragments.Count] : FragmentCoElution.CorrelationsTo(traces, reference, from, to);
            double libraryTotal = libraryIntensities.Sum();
            double weightedCoElution = libraryIntensities.Select((w, f) => w * allCorrelations[f]).Sum() / libraryTotal;
            double weightedMatchedFraction = libraryIntensities.Select((w, f) => apexIntensities[f] > 0 ? w : 0).Sum() / libraryTotal;
            int[] byRank = Enumerable.Range(0, fragments.Count).OrderByDescending(f => libraryIntensities[f]).ThenBy(f => f).ToArray();
            int[] top3 = byRank.Take(3).ToArray();
            double top3CoElution = top3.Average(f => allCorrelations[f]);
            double top3Cosine = SpectralSimilarity.CosineOfAlignedVectors(top3.Select(f => apexIntensities[f]).ToArray(), top3.Select(f => libraryIntensities[f]).ToArray());
            double top1Present = apexIntensities[byRank[0]] > 0 ? 1 : 0;
            double top2Present = byRank.Length > 1 && apexIntensities[byRank[1]] > 0 ? 1 : 0;
            double top1Ppm = _parameters.FragmentTolerancePpm;
            if (top1Present > 0)
            {
                var top1 = fragments[byRank[0]];
                double observedTop1 = apexSpectrum.XArray[apexSpectrum.GetClosestPeakIndex(top1.Mz)];
                top1Ppm = Math.Abs(observedTop1 - top1.Mz) / top1.Mz * 1e6;
            }

            // Detectable fragments: the library pattern scaled to the apex by least squares, against the apex scan's noise
            // (its 25th-percentile intensity); fragments expected below 3x noise are not held against the precursor
            double libraryScale = libraryIntensities.Select((l, f) => l * apexIntensities[f]).Sum() / libraryIntensities.Sum(l => l * l);
            double noise = apexSpectrum.Size > 0 ? apexSpectrum.YArray.Order().ElementAt(apexSpectrum.Size / 4) : 0;
            int[] detectable = Enumerable.Range(0, fragments.Count).Where(f => libraryScale * libraryIntensities[f] >= 3 * noise).ToArray();
            if (detectable.Length == 0)
                detectable = [byRank[0]];
            double detectableMatchedFraction = detectable.Count(f => apexIntensities[f] > 0) / (double)detectable.Length;
            double detectableCoElution = detectable.Average(f => allCorrelations[f]);

            // Fragment isotopes at the apex (OpenSWATH): a matched peak with a bigger peak one isotope below is probably another
            // ion's isotope; a real fragment usually shows its own M+1
            double isotopeOverlap = 0, fragmentM1 = 0, matchedWeight = 0;
            for (int f = 0; f < fragments.Count; f++)
            {
                if (apexIntensities[f] <= 0)
                    continue;
                double w = libraryIntensities[f];
                matchedWeight += w;
                double spacing = 1.0033548 / Math.Max(1, fragments[f].Charge);
                int below = apexSpectrum.GetClosestPeakIndex(fragments[f].Mz - spacing);
                if (tolerance.Within(apexSpectrum.XArray[below], fragments[f].Mz - spacing) && apexSpectrum.YArray[below] > apexIntensities[f])
                    isotopeOverlap += w;
                int above = apexSpectrum.GetClosestPeakIndex(fragments[f].Mz + spacing);
                if (tolerance.Within(apexSpectrum.XArray[above], fragments[f].Mz + spacing))
                    fragmentM1 += w;
            }
            if (matchedWeight > 0)
            {
                isotopeOverlap /= matchedWeight;
                fragmentM1 /= matchedWeight;
            }

            // Specificity (PECAN's idea): a match where peaks are sparse is less likely by chance than one in the crowded
            // low-m/z region. Chance probability p = 1 - exp(-density x window), with density the apex spectrum's peaks per
            // Th within +-5 Th and window the tolerance width in Th; each fragment's evidence is its surprisal -log10 p
            double matchSurprisal = 0, possibleSurprisal = 0, specificityWeightedCoElution = 0;
            for (int f = 0; f < fragments.Count; f++)
            {
                double mz = fragments[f].Mz;
                int lo = apexSpectrum.Size == 0 ? 0 : LowerBound(apexSpectrum.XArray, mz - 5);
                int hi = apexSpectrum.Size == 0 ? 0 : LowerBound(apexSpectrum.XArray, mz + 5);
                double density = Math.Max(1, hi - lo) / 10.0;
                double window = 2 * mz * _parameters.FragmentTolerancePpm * 1e-6;
                double chance = Math.Clamp(1 - Math.Exp(-density * window), 1e-6, 1);
                double surprisal = -Math.Log10(chance);
                possibleSurprisal += surprisal;
                specificityWeightedCoElution += surprisal * allCorrelations[f];
                if (apexIntensities[f] > 0)
                    matchSurprisal += surprisal;
            }
            double surprisalFraction = possibleSurprisal > 0 ? matchSurprisal / possibleSurprisal : 0;
            specificityWeightedCoElution = possibleSurprisal > 0 ? specificityWeightedCoElution / possibleSurprisal : 0;

            // MS1 as DIA-NN uses it: the precursor's MS1 trace (monoisotopic and M+1) correlated with the fragment profile,
            // never its raw intensity
            double ms1Correlation = 0, ms1IsotopeCorrelation = 0;
            if (reference is not null && ms1Mono is not null && ms1Isotope is not null)
            {
                ms1Correlation = FragmentCoElution.CorrelationsTo([ms1Mono], reference, from, to)[0];
                ms1IsotopeCorrelation = FragmentCoElution.CorrelationsTo([ms1Isotope], reference, from, to)[0];
            }

            // Uniqueness: this apex against the best competing scan outside its co-elution window, its z-score among the
            // window's scans, and the share of the window's fragment signal inside its peak
            double apexScore = apexScores[apex];
            double competitor = Enumerable.Range(0, apexScores.Length).Where(s => Math.Abs(s - apex) > _parameters.ApexHalfWidthScans)
                .Select(s => apexScores[s]).DefaultIfEmpty(0).Max();
            double apexScoreDelta = apexScore > 0 ? Math.Clamp((apexScore - competitor) / apexScore, -1, 1) : 0;
            double apexScoreZ = scoreSd > 0 ? (apexScore - scoreMean) / scoreSd : 0;
            double peakSignalFraction = windowSignal > 0 ? Enumerable.Range(peakStart, peakEnd - peakStart + 1).Sum(s => summed[s]) / windowSignal : 0;
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
                .. fragmentCorrelations,
                tightCoElution,
                remainingCoElution,
                windowCosine,
                medianCoElution,
                goodFragments,
                veryGoodFragment,
                gaussianFit,
                areaSqrtCosine,
                areaPearson,
                areaManhattan,
                apexScoreDelta,
                apexScoreZ,
                peakSignalFraction,
                weightedCoElution,
                weightedMatchedFraction,
                top3CoElution,
                top3Cosine,
                top1Present,
                top2Present,
                top1Ppm,
                detectableMatchedFraction,
                detectableCoElution,
                detectable.Length,
                ms1Correlation,
                ms1IsotopeCorrelation,
                isotopeOverlap,
                fragmentM1,
                matchSurprisal,
                surprisalFraction,
                specificityWeightedCoElution,
            ];

            yield return new DiaPrecursorMatch(
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
    }

    /// <summary>The first index whose value is at least <paramref name="x"/> in an ascending array.</summary>
    private static int LowerBound(double[] sorted, double x)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) / 2;
            if (sorted[mid] < x) lo = mid + 1; else hi = mid;
        }
        return lo;
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

    /// <summary>Per-scan median of the traces, each normalised to sum 1 over [from, to]; zero outside the range.</summary>
    private static double[] MedianProfile(IReadOnlyList<double[]> traces, int from, int to, int length)
    {
        var median = new double[length];
        var sums = traces.Select(t => Enumerable.Range(from, to - from + 1).Sum(s => t[s])).ToArray();
        var values = new double[traces.Count];
        for (int s = from; s <= to; s++)
        {
            for (int f = 0; f < traces.Count; f++)
                values[f] = sums[f] > 0 ? traces[f][s] / sums[f] : 0;
            var sorted = values.Order().ToArray();
            median[s] = sorted.Length % 2 == 1 ? sorted[sorted.Length / 2] : 0.5 * (sorted[sorted.Length / 2 - 1] + sorted[sorted.Length / 2]);
        }
        return median;
    }

    /// <summary>Trapezoid area of the trace over [start, end], less the baseline at its lower end point.</summary>
    private static double Area(double[] trace, int start, int end)
    {
        double area = 0;
        for (int s = start; s < end; s++)
            area += 0.5 * (trace[s] + trace[s + 1]);
        double baseline = Math.Min(trace[start], trace[end]) * (end - start);
        return Math.Max(0, area - baseline);
    }

    /// <summary>
    /// Combines each match's features into one score with mzLib's semi-supervised <see cref="TargetDecoyRescorer"/>, with folds
    /// grouped by sequence so that no match is scored by a model trained on it. Each precursor then keeps its best-scoring
    /// candidate, and (D+1)/T q-values are assigned among targets.
    /// If rescoring cannot run (for example too few matches), the pre-rescoring score stands.
    /// </summary>
    private List<DiaPrecursorMatch> AssignQValues(List<DiaPrecursorMatch> matches)
    {
        if (matches.Count > 0)
        {
            var rescored = TargetDecoyRescorer.Score(
                matches.Select(m => m.Features).ToList(),
                matches.Select(m => m.IsDecoy).ToList(),
                matches.Select(m => m.FullSequence).ToList(),
                positiveQValue: _parameters.ClassifierTrainingQValue,
                model: _parameters.ClassifierModel,
                candidateGroups: matches.Select(m => m.PrecursorIndex).ToList());
            if (rescored.Scores.All(double.IsFinite))
                matches = matches.Select((m, i) => m with { Score = rescored.Scores[i] }).ToList();
        }

        // One match per precursor: its best-scoring candidate apex, across windows
        matches = matches.GroupBy(m => m.PrecursorIndex).Select(g => g.OrderByDescending(m => m.Score).ThenBy(m => m.ApexRt.Value).First()).ToList();

        var targets = matches.Where(m => !m.IsDecoy).ToList();
        var decoyScores = matches.Where(m => m.IsDecoy).Select(m => m.Score).ToList();
        double[] qValues = DeconvolutionQValueCalculator.AssignQValues(targets.Select(m => m.Score).ToList(), decoyScores);

        var assigned = targets.Select((m, i) => m with { QValue = qValues[i] }).ToList();
        assigned.AddRange(matches.Where(m => m.IsDecoy));
        return assigned;
    }
}

/// <summary>
/// The precursors a <see cref="DiaLibrarySearchEngine"/> scored, targets with their q-values, and the peptides they
/// collapse to, with peptide-level q-values from <see cref="DiaPeptideFdr"/>.
/// </summary>
public class DiaLibrarySearchResults(DiaLibrarySearchEngine engine, List<DiaPrecursorMatch> matches, List<DiaPeptideMatch> peptides)
    : MetaMorpheusEngineResults(engine)
{
    public List<DiaPrecursorMatch> Matches { get; init; } = matches;

    public List<DiaPeptideMatch> Peptides { get; init; } = peptides;

    public int TargetCount => Matches.Count(m => !m.IsDecoy);

    public int DecoyCount => Matches.Count(m => m.IsDecoy);

    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine(base.ToString());
        sb.AppendLine($"Target precursors: {TargetCount}");
        sb.AppendLine($"Decoy precursors: {DecoyCount}");
        sb.AppendLine($"Target precursors with q-value <= 0.01: {Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01)}");
        sb.AppendLine($"Target peptides with q-value <= 0.01: {Peptides.Count(p => !p.IsDecoy && p.QValue <= 0.01)}");
        return sb.ToString();
    }
}
