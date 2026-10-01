#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer.DiaLibrarySearch;
using MzLibUtil;
using NUnit.Framework;

namespace Test.DiaLibrarySearch;

/// <summary>
/// Interference removal, after DIA-NN's remove_ifs: a match whose top library fragments a better-scoring match co-eluting
/// in the same or a neighbouring window already explains is dropped, target or decoy alike. On the whole-proteome HeLa
/// library such "riders" were most of the decoys scoring above the 2% target threshold.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaInterferenceRemovalTests
{
    private const double Ppm = 20;
    private const double RtTolerance = 0.1;
    private const int Explained = 4;

    // Precursor 0's fragments; the others copy some of them
    private static readonly float[] Base = [300f, 400f, 500f, 600f, 700f, 800f];

    private static float[] Sharing(int shared) => Base.Take(shared).Concat(Enumerable.Range(0, 6 - shared).Select(i => 1100f + 37 * i)).ToArray();

    private static DiaPrecursorMatch Match(int index, double score, bool isDecoy = false, double rt = 10, double mz = 500) =>
        new(index, $"PEPTIDE{index}K", 2, mz, isDecoy, new Irt(0), new RtMinutes(rt), new Irt(0), score, [score]);

    // 200 unrelated, well-separated targets keep the provisional floor low, so the cases under test take part
    private static List<DiaPrecursorMatch> Background() =>
        Enumerable.Range(100, 200).Select(i => Match(i, 1 + i / 1000.0, rt: 20 + i)).ToList();

    private static List<DiaPrecursorMatch> Remove(List<DiaPrecursorMatch> matches, Dictionary<int, float[]> fragments) =>
        DiaInterferenceRemoval.Remove(matches,
            index => fragments.TryGetValue(index, out var mzs) ? mzs : Enumerable.Range(0, 6).Select(i => 2000f + index * 7 + i).ToArray(),
            mz => (int)((mz - 400) / 25), RtTolerance, Ppm, Explained);

    [Test]
    public void ADecoyRidingOnABetterTargetIsRemovedAndTheTargetKept()
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9, isDecoy: true)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(5) });

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Contain(0).And.Not.Contain(1));
        Assert.That(kept, Has.Count.EqualTo(matches.Count - 1));
    }

    /// <summary>Removal is blind to the label: a target explained by a better decoy goes too.</summary>
    [Test]
    public void ATargetExplainedByABetterDecoyIsRemovedToo()
    {
        var matches = Background().Append(Match(0, 5, isDecoy: true)).Append(Match(1, 4.9)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(4) });

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Contain(0).And.Not.Contain(1));
    }

    [TestCase(3, true, TestName = "Three shared fragments are coincidence")]
    [TestCase(4, false, TestName = "Four shared fragments are interference")]
    public void TheThresholdIsOnSharedTopFragments(int shared, bool survives)
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(shared) });

        Assert.That(kept.Any(m => m.PrecursorIndex == 1), Is.EqualTo(survives));
    }

    [TestCase(10.05, 500, false, TestName = "Co-eluting in the same window")]
    [TestCase(10.5, 500, true, TestName = "Eluting apart")]
    [TestCase(10.05, 530, false, TestName = "In the neighbouring window")]
    [TestCase(10.05, 560, true, TestName = "Two windows away")]
    public void OnlyACoElutingMatchInTheSameOrANeighbouringWindowExplains(double rt, double mz, bool survives)
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9, rt: rt, mz: mz)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(6) });

        Assert.That(kept.Any(m => m.PrecursorIndex == 1), Is.EqualTo(survives));
    }

    /// <summary>The better match explains; the worse one, even when it copies the better, never removes it.</summary>
    [Test]
    public void TheBetterScoringMatchIsTheOneKept()
    {
        var matches = Background().Append(Match(0, 4.9)).Append(Match(1, 5)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(6) });

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Contain(1).And.Not.Contain(0));
    }

    /// <summary>Matches too weak to be reported take no part, either as explainers or as explained.</summary>
    [Test]
    public void MatchesBelowTheProvisionalFloorAreLeftAlone()
    {
        var decoys = Enumerable.Range(400, 200).Select(i => Match(i, -5 + i / 1000.0, isDecoy: true, rt: 50 + i)).ToList();
        var matches = Background().Concat(decoys).Append(Match(0, -9)).Append(Match(1, -9.1)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(6) });

        Assert.That(kept, Has.Count.EqualTo(matches.Count));
    }

    [Test]
    public void ArgumentsAreChecked()
    {
        Assert.Throws<ArgumentNullException>(() => DiaInterferenceRemoval.Remove(null!, _ => [], _ => 0, 1, 20, 4));
        Assert.That(DiaInterferenceRemoval.Remove([], _ => [], _ => 0, 1, 20, 4), Is.Empty);
    }
}
