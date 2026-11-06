#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer.DiaLibrarySearch;
using NUnit.Framework;

namespace Test.DiaLibrarySearch;

/// <summary>
/// What a run's confident targets say of the library's iRT, by sequence composition. A predicted iRT is fixed for some other column, and on both benchmark runs what is left after the run's RT->iRT
/// curve has a shape: K, R and H elute 1.5 to 3 iRT earlier than predicted. A fit on 1,000 anchors takes the median error of the targets a search scores too low from 2.47 to 1.59 and from 2.22 to 1.25.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaIrtCompositionModelTests
{
    private const string Acids = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>Random tryptic-looking peptides, each with the residual the rule gives it.</summary>
    private static List<DiaIrtCompositionModel.Anchor> Anchors(int count, Func<string, int, double> residual, int seed = 1)
    {
        var random = new Random(seed);
        return Enumerable.Range(0, count).Select(_ =>
        {
            int length = random.Next(7, 25);
            string sequence = new string(Enumerable.Range(0, length - 1).Select(_ => Acids[random.Next(Acids.Length)]).ToArray()) + (random.Next(2) == 0 ? 'K' : 'R');
            int charge = random.Next(2, 4);
            return new DiaIrtCompositionModel.Anchor(sequence, sequence, charge, residual(sequence, charge));
        }).ToList();
    }

    private static int Count(string sequence, char residue) => sequence.Count(c => c == residue);

    [Test]
    public void ItLearnsWhatEachResidueDoes()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(2000, (s, _) => -3 * Count(s, 'K') + 2 * Count(s, 'P')));

        Assert.That(model.Offset("AAAAAAAAKK", "AAAAAAAAKK", 2) - model.Offset("AAAAAAAAAK", "AAAAAAAAAK", 2), Is.EqualTo(-3).Within(0.2));
        Assert.That(model.Offset("APAAAAAAAK", "APAAAAAAAK", 2) - model.Offset("AAAAAAAAAK", "AAAAAAAAAK", 2), Is.EqualTo(2).Within(0.2));
    }

    /// <summary>A permuted decoy shares its target's composition, so it moves as its target does, and the target-decoy estimate stays fair.</summary>
    [Test]
    public void APermutationOfTheSameCompositionMovesAlike()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(2000, (s, _) => -3 * Count(s, 'H') + Count(s, 'V')));

        Assert.That(model.Offset("AVGHLTEK", "AVGHLTEK", 2), Is.EqualTo(model.Offset("AGVLHTEK", "AGVLHTEK", 2)).Within(1e-9));
    }

    [Test]
    public void ItLearnsWhatChargeDoes()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(2000, (_, z) => z == 3 ? 2 : 0));

        Assert.That(model.Offset("AAAAAAAAAK", "AAAAAAAAAK", 3) - model.Offset("AAAAAAAAAK", "AAAAAAAAAK", 2), Is.EqualTo(2).Within(0.2));
    }

    [Test]
    public void ItLearnsWhatLengthDoes()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(4000, (s, _) => s.Length <= 9 ? 3 : 0));

        Assert.That(model.Offset(new string('A', 7) + "K", new string('A', 7) + "K", 2), Is.GreaterThan(model.Offset(new string('A', 15) + "K", new string('A', 15) + "K", 2) + 1));
    }

    /// <summary>Basic and acidic residues at the N-terminus do not read as they do further in: the free amine is right beside them.</summary>
    [TestCase("KRH", 'K')]
    [TestCase("DE", 'D')]
    public void ItLearnsWhatTheFirstResidueDoes(string group, char residue)
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(4000, (s, _) => group.Contains(s[0]) ? 3 : 0));
        string first = residue + "AAAAAAAAK", later = "AAAAA" + residue + "AAAK";

        Assert.That(model.Offset(first, first, 2) - model.Offset(later, later, 2), Is.EqualTo(3).Within(0.6));
    }

    [Test]
    public void ItLearnsWhatABasicSecondResidueDoes()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(4000, (s, _) => "KRH".Contains(s[1]) ? 3 : 0));

        Assert.That(model.Offset("ARAAAAAAAK", "ARAAAAAAAK", 2) - model.Offset("AAAAARAAAK", "AAAAARAAAK", 2), Is.EqualTo(3).Within(0.6));
    }

    [Test]
    public void AnOxidationIsCountedFromTheModifiedSequence()
    {
        var anchors = Anchors(2000, (_, _) => 0).Select((a, i) => i % 2 == 0 ? a : a with { FullSequence = a.BaseSequence + "[Common Variable:Oxidation on M]", Residual = -4 }).ToList();
        var model = DiaIrtCompositionModel.Fit(anchors);

        Assert.That(model.Offset("AMAAAK", "AM[Common Variable:Oxidation on M]AAAK", 2) - model.Offset("AMAAAK", "AMAAAK", 2), Is.EqualTo(-4).Within(0.2));
    }

    /// <summary>It only says what is left after the run's curve, so a run with nothing left moves nothing.</summary>
    [Test]
    public void WithNothingLeftNothingMoves()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(1000, (_, _) => 0));

        Assert.That(model.Offset("AVGHLTEK", "AVGHLTEK", 2), Is.EqualTo(0).Within(1e-9));
    }

    [Test]
    public void TheMeanOfWhatIsLeftIsLearnedToo()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(1000, (_, _) => 1.5));

        Assert.That(model.Offset("AVGHLTEK", "AVGHLTEK", 2), Is.EqualTo(1.5).Within(0.1));
    }

    [Test]
    public void TooFewAnchorsMoveNothing()
    {
        int below = DiaIrtCompositionModel.MinimumAnchors - 1;

        Assert.That(DiaIrtCompositionModel.Fit(Anchors(below, (s, _) => 5 * Count(s, 'K'))).Offset("AAAAKK", "AAAAKK", 2), Is.EqualTo(0));
        Assert.That(DiaIrtCompositionModel.Fit(Anchors(below + 1, (s, _) => 5 * Count(s, 'K'))).Offset("AAAAKK", "AAAAKK", 2), Is.Not.EqualTo(0));
        Assert.That(DiaIrtCompositionModel.MinimumAnchors, Is.EqualTo(200));
    }

    /// <summary>A residue no anchor holds has nothing to say, and a library holding it must not be thrown far off.</summary>
    [Test]
    public void AResidueNoAnchorHoldsDoesNothing()
    {
        var anchors = Anchors(2000, (s, _) => -2 * Count(s, 'K')).Where(a => !a.BaseSequence.Contains('W')).ToList();
        var model = DiaIrtCompositionModel.Fit(anchors);

        Assert.That(model.Offset("AAAWAAAK", "AAAWAAAK", 2) - model.Offset("AAAAAAAK", "AAAAAAAK", 2), Is.EqualTo(0).Within(0.2), "it has one A fewer, and no term of its own for W");
    }

    [Test]
    public void AnUnknownResidueIsIgnored()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(2000, (s, _) => Count(s, 'K')));

        Assert.That(double.IsFinite(model.Offset("XXXXXXXK", "XXXXXXXK", 2)), Is.True);
        Assert.That(model.Offset("XXXXXXXK", "XXXXXXXK", 2) - model.Offset("XXXXXXXX", "XXXXXXXX", 2), Is.EqualTo(1).Within(0.2), "an unknown residue counts as nothing, and the K still does");
    }

    /// <summary>A wild anchor is one in 2,000: it must not move the fit much, since calibration's anchors are confident but not all right.</summary>
    [Test]
    public void OneWildAnchorDoesNotMoveItFar()
    {
        var anchors = Anchors(2000, (s, _) => -2 * Count(s, 'K'));
        anchors[0] = anchors[0] with { Residual = 500 };
        var model = DiaIrtCompositionModel.Fit(anchors);

        Assert.That(model.Offset("AAAAAAAAKK", "AAAAAAAAKK", 2) - model.Offset("AAAAAAAAAK", "AAAAAAAAAK", 2), Is.EqualTo(-2).Within(0.3));
    }

    /// <summary>The library is moved 5 iRT earlier per W after the run was made, so the run reads those precursors 5 iRT later than it says, and calibration's fit says so.</summary>
    [Test]
    public void CalibrationFitsItFromItsAnchors()
    {
        string directory = System.IO.Path.Combine(TestContext.CurrentContext.WorkDirectory, $"DiaIrtCompositionModelTests_{Guid.NewGuid():N}");
        System.IO.Directory.CreateDirectory(directory);
        try
        {
            var run = SyntheticDiaRun.Build(1000, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 60);
            foreach (var entry in run.Library.Where(e => e.BaseSequence.Contains('W')))
                entry.RetentionTime -= 5 * entry.BaseSequence.Count(c => c == 'W');
            using var library = Readers.SpectralLibrary.MslLibrary.Load(run.WriteLibrary(directory));

            var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new EngineLayer.CommonParameters());

            Assert.That(calibration.IrtComposition, Is.Not.Null);
            double offset = calibration.IrtComposition!.Offset("AAAAWAAAAK", "AAAAWAAAAK", 2) - calibration.IrtComposition.Offset("AAAAAAAAAK", "AAAAAAAAAK", 2);
            Assert.That(offset, Is.EqualTo(5).Within(1.5));

            var oneRound = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new EngineLayer.CommonParameters(), rounds: 1);
            Assert.That(oneRound.IrtComposition, Is.Null, "as with the MS1 offset, only a later round fits one");
        }
        finally
        {
            System.IO.Directory.Delete(directory, true);
        }
    }

    [Test]
    public void CalibrationIgnoresAModelPassedToIt()
    {
        var model = DiaIrtCompositionModel.Fit([]);

        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(IrtComposition: model)).IrtComposition, Is.Null);
    }

    [Test]
    public void AnEmptySequenceDoesNotThrow()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(1000, (s, _) => Count(s, 'K')));

        Assert.That(double.IsFinite(model.Offset("", "", 2)), Is.True);
        Assert.That(double.IsFinite(model.Offset("K", "K", 2)), Is.True);
    }
}
