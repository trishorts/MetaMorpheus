#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using StatisticalModels;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// Removal of interfering precursors, after DIA-NN (remove_ifs): a match whose most intense library fragments are already
/// explained by a better-scoring match co-eluting in the same or a neighbouring isolation window is dropped. Such a match
/// rides on the other's signal. Targets and decoys are treated alike, so a false target riding on a real peptide goes the
/// same way as a decoy would, and the target-decoy estimate stays fair.
/// </summary>
public static class DiaInterferenceRemoval
{
    /// <summary>A match's most intense library fragments that are checked against a better match.</summary>
    public const int CheckedFragments = 6;

    /// <summary>Only matches scoring at least as well as the target at this provisional q-value take part; the rest
    /// cannot be reported anyway, and leaving them out keeps the pass cheap.</summary>
    public const double ProvisionalQValue = 0.10;

    /// <param name="matches">One match per precursor, already scored.</param>
    /// <param name="fragmentMzsByIntensity">A precursor's library fragment m/z, most intense first.</param>
    /// <param name="windowOf">The index of the isolation window a precursor m/z falls in.</param>
    /// <param name="rtToleranceMinutes">How far apart two apexes may be and still co-elute.</param>
    /// <param name="explainedFragments">How many of the checked fragments a better match must explain.</param>
    /// <returns>The matches that are not explained by a better one, in their input order.</returns>
    public static List<DiaPrecursorMatch> Remove(IReadOnlyList<DiaPrecursorMatch> matches, Func<int, float[]> fragmentMzsByIntensity,
        Func<double, int> windowOf, double rtToleranceMinutes, double tolerancePpm, int explainedFragments, Action<string>? report = null)
    {
        ArgumentNullException.ThrowIfNull(matches);
        ArgumentNullException.ThrowIfNull(fragmentMzsByIntensity);
        ArgumentNullException.ThrowIfNull(windowOf);
        if (matches.Count == 0)
            return [];

        double[] q = TargetDecoyQValues.Compute(matches.Select(m => m.Score).ToArray(), matches.Select(m => m.IsDecoy).ToArray());
        double floor = matches.Where((m, i) => q[i] <= ProvisionalQValue).Select(m => m.Score).DefaultIfEmpty(double.PositiveInfinity).Min();

        var clock = System.Diagnostics.Stopwatch.StartNew();
        long neighbours = 0;
        int considered = 0;
        var removed = new HashSet<int>();
        var kept = new Dictionary<int, List<(double Rt, float[] Mzs)>>();
        foreach (var match in matches.Where((m, i) => q[i] <= ProvisionalQValue).OrderByDescending(m => m.Score).ThenBy(m => m.PrecursorIndex))
        {
            considered++;
            int window = windowOf(match.PrecursorMz);
            float[] own = fragmentMzsByIntensity(match.PrecursorIndex);
            float[] checkedMzs = own.Take(CheckedFragments).ToArray();
            bool explained = false;
            for (int w = window - 1; w <= window + 1 && !explained; w++)
            {
                if (!kept.TryGetValue(w, out var better))
                    continue;
                foreach (var (rt, mzs) in better)
                {
                    if (Math.Abs(rt - match.ApexRt.Value) > rtToleranceMinutes)
                        continue;
                    neighbours++;
                    int count = checkedMzs.Count(mz => mzs.Any(other => Math.Abs(other - mz) <= mz * tolerancePpm * 1e-6));
                    if (count >= explainedFragments)
                    {
                        explained = true;
                        break;
                    }
                }
            }
            if (explained)
            {
                removed.Add(match.PrecursorIndex);
                continue;
            }
            if (!kept.TryGetValue(window, out var list))
                kept[window] = list = [];
            list.Add((match.ApexRt.Value, own));
        }
        report?.Invoke($"interference removal: {considered} of {matches.Count} matches above the floor {floor:G4}; RT tolerance {rtToleranceMinutes:F4} min; " +
            $"{neighbours} co-eluting pairs checked; {removed.Count} removed ({matches.Count(m => m.IsDecoy && removed.Contains(m.PrecursorIndex))} decoys)  [{clock.Elapsed:mm\\:ss}]");
        return matches.Where(m => !removed.Contains(m.PrecursorIndex)).ToList();
    }
}
