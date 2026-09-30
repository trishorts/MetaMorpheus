#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using StatisticalModels;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>A peptide found by the DIA search: its best-scoring precursor, and a peptide-level q-value.</summary>
/// <param name="BestPrecursorIndex">Library index of the precursor (charge state) whose score the peptide carries.</param>
/// <param name="QValue">Target-decoy q-value among peptides; for a decoy, the value at its rank.</param>
public sealed record DiaPeptideMatch(string FullSequence, bool IsDecoy, int BestPrecursorIndex, double Score, double QValue);

/// <summary>
/// Peptide-level FDR for the DIA search, following MetaMorpheus's DDA convention (FdrAnalysisEngine): collapse to
/// peptides by full sequence, keeping each peptide's best-scoring precursor, then compute q-values afresh among
/// peptides. A peptide seen at two charges counts once, never twice.
/// </summary>
public static class DiaPeptideFdr
{
    /// <summary>
    /// One peptide per full sequence and decoy flag, ordered by full sequence. Equal best scores go to the lower
    /// precursor index. q-values are mzLib's shared (D+1)/T (<see cref="TargetDecoyQValues"/>).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="precursors"/> is null.</exception>
    public static List<DiaPeptideMatch> Assign(IReadOnlyList<DiaPrecursorMatch> precursors)
    {
        ArgumentNullException.ThrowIfNull(precursors);

        var best = precursors
            .GroupBy(m => (m.FullSequence, m.IsDecoy))
            .Select(group => group.OrderByDescending(m => m.Score).ThenBy(m => m.PrecursorIndex).First())
            .OrderBy(m => m.FullSequence, StringComparer.Ordinal).ThenBy(m => m.IsDecoy)
            .ToList();
        double[] qValues = TargetDecoyQValues.Compute(best.Select(m => m.Score).ToArray(), best.Select(m => m.IsDecoy).ToArray());
        return best.Select((m, i) => new DiaPeptideMatch(m.FullSequence, m.IsDecoy, m.PrecursorIndex, m.Score, qValues[i])).ToList();
    }
}
