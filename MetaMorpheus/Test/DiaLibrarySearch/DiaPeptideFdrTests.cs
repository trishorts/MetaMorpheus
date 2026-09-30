#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer.DiaLibrarySearch;
using MzLibUtil;
using NUnit.Framework;
using StatisticalModels;

namespace Test.DiaLibrarySearch;

/// <summary>
/// Peptide-level FDR for the DIA search. Precursors (charge states) collapse to peptides, each keeping its best-scoring
/// precursor. q-values then come from a fresh target-decoy competition among peptides, never inherited from precursors,
/// as MetaMorpheus's FdrAnalysisEngine does for DDA.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaPeptideFdrTests
{
    private static DiaPrecursorMatch Match(int index, string sequence, int charge, bool isDecoy, double score) =>
        new(index, sequence, charge, 500, isDecoy, new Irt(0), new RtMinutes(10), new Irt(0), score, [score], QValue: 0.5);

    [Test]
    public void ChargeStatesCollapseToThePeptidesBestPrecursor()
    {
        var precursors = new List<DiaPrecursorMatch>
        {
            Match(0, "PEPTIDEK", 2, false, 3.0),
            Match(1, "PEPTIDEK", 3, false, 7.0),
            Match(2, "EDITPEPK", 2, true, 1.0),
        };

        var peptides = DiaPeptideFdr.Assign(precursors);

        Assert.That(peptides, Has.Count.EqualTo(2));
        var peptide = peptides.Single(p => p.FullSequence == "PEPTIDEK");
        Assert.That(peptide.Score, Is.EqualTo(7.0));
        Assert.That(peptide.BestPrecursorIndex, Is.EqualTo(1));
        Assert.That(peptide.IsDecoy, Is.False);
    }

    /// <summary>
    /// The q-values are those of a fresh competition among peptides: the shared (D+1)/T over the peptides' best scores. Two
    /// charge states of a target no longer count twice.
    /// </summary>
    [Test]
    public void PeptideQValuesComeFromAFreshCompetitionAmongPeptides()
    {
        var random = new Random(3);
        var precursors = new List<DiaPrecursorMatch>();
        for (int i = 0; i < 300; i++)
        {
            bool decoy = i % 3 == 0;
            string sequence = $"PEP{i:D4}K";
            precursors.Add(Match(2 * i, sequence, 2, decoy, random.NextDouble() + (decoy ? 0 : 0.5)));
            precursors.Add(Match(2 * i + 1, sequence, 3, decoy, random.NextDouble() + (decoy ? 0 : 0.5)));
        }

        var peptides = DiaPeptideFdr.Assign(precursors);

        var expected = TargetDecoyQValues.Compute(peptides.Select(p => p.Score).ToArray(), peptides.Select(p => p.IsDecoy).ToArray());
        Assert.That(peptides, Has.Count.EqualTo(300));
        for (int i = 0; i < peptides.Count; i++)
            Assert.That(peptides[i].QValue, Is.EqualTo(expected[i]).Within(1e-12), peptides[i].FullSequence);
        Assert.That(peptides.All(p => p.Score == precursors.Where(m => m.FullSequence == p.FullSequence).Max(m => m.Score)));
    }

    /// <summary>Equal best scores go to the lower precursor index, so the choice never depends on input order.</summary>
    [Test]
    public void TiedChargeStatesKeepTheLowerPrecursorIndex()
    {
        var precursors = new List<DiaPrecursorMatch> { Match(9, "PEPTIDEK", 3, false, 5.0), Match(4, "PEPTIDEK", 2, false, 5.0) };

        Assert.That(DiaPeptideFdr.Assign(precursors).Single().BestPrecursorIndex, Is.EqualTo(4));
        Assert.That(DiaPeptideFdr.Assign(precursors.AsEnumerable().Reverse().ToList()).Single().BestPrecursorIndex, Is.EqualTo(4));
    }

    [Test]
    public void NoPrecursorsMeansNoPeptides()
    {
        Assert.That(DiaPeptideFdr.Assign([]), Is.Empty);
        Assert.Throws<ArgumentNullException>(() => DiaPeptideFdr.Assign(null!));
    }
}
