#nullable enable
using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DiaLibrarySearch;
using MzLibUtil;
using NUnit.Framework;
using Readers.SpectralLibrary;

namespace Test.DiaLibrarySearch;

/// <summary>
/// A DIA run must be searchable with nothing but its scans and a library. It calibrates its own minutes onto the
/// library's iRT from a wide first pass, and nothing is borrowed from another search engine.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaIrtSelfCalibrationTests
{
    private string _directory = "";

    [OneTimeSetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DiaIrtSelfCalibrationTests");
        Directory.CreateDirectory(_directory);
    }

    [OneTimeTearDown]
    public void TearDown()
    {
        Directory.Delete(_directory, true);
    }

    /// <summary>
    /// Calibration scores every fragment whatever the search's fragment rule: with the rule in calibration too, the held-out
    /// run lost 2.1% (two-seed means at a matched paired FDP of 1%). Calibration's parameters are the search's with the rule off,
    /// and nothing else changed; and calibration gives the same result whichever way the search sets the rule.
    /// </summary>
    [Test]
    public void CalibrationIgnoresTheFragmentRule()
    {
        var search = new DiaLibrarySearchParameters(DiaNnFragmentFilter: true, TopFragmentCount: 10);
        var calibrationParameters = DiaIrtSelfCalibration.CalibrationParameters(search);
        Assert.That(calibrationParameters.DiaNnFragmentFilter, Is.False);
        Assert.That(calibrationParameters with { DiaNnFragmentFilter = true, Ms1TolerancePpm = search.Ms1TolerancePpm, DiaNnSignalShare = search.DiaNnSignalShare,
                InterferenceExplainedFragments = search.InterferenceExplainedFragments},
            Is.EqualTo(search), "nothing else changes but calibration's own MS1 tolerance, the signal share and the interference rule");

        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var on = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(DiaNnFragmentFilter: true), new CommonParameters());
        var off = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(DiaNnFragmentFilter: false), new CommonParameters());
        Assert.That(on.Anchors, Is.EqualTo(off.Anchors));
    }

    /// <summary>
    /// Calibration scores without DIA-NN's signal share, whatever the search does: the share was benchmarked in the main search
    /// only (+0.9% / +0.6%), and switching it on in calibration too moved the anchors and cost HF-X 1.5% on seed 0 (58,284
    /// against 59,170). Calibration gives the same anchors whichever way the search sets it.
    /// </summary>
    [Test]
    public void CalibrationIgnoresTheSignalShare()
    {
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(DiaNnSignalShare: true)).DiaNnSignalShare, Is.False);

        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var on = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(DiaNnSignalShare: true), new CommonParameters());
        var off = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(DiaNnSignalShare: false), new CommonParameters());
        Assert.That(on.Anchors, Is.EqualTo(off.Anchors));
    }

    /// <summary>
    /// Calibration removes interference at 4 explained fragments, whatever the search uses: the main search's 3 was
    /// benchmarked with calibration at 4. Calibration gives the same anchors whichever rule the search sets.
    /// </summary>
    [Test]
    public void CalibrationKeepsItsOwnInterferenceRule()
    {
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(InterferenceExplainedFragments: 3)).InterferenceExplainedFragments,
            Is.EqualTo(DiaIrtSelfCalibration.CalibrationInterferenceExplainedFragments));
        Assert.That(DiaIrtSelfCalibration.CalibrationInterferenceExplainedFragments, Is.EqualTo(4));
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(InterferenceSameMzOnly: true)).InterferenceSameMzOnly, Is.True,
            "the same-m/z pairing was benchmarked with calibration pairing the same way");

        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 3000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var three = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(InterferenceExplainedFragments: 3), new CommonParameters());
        var six = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(InterferenceExplainedFragments: 6), new CommonParameters());
        Assert.That(three.Anchors, Is.EqualTo(six.Anchors));
    }

    /// <summary>
    /// Calibration also fits the run's MS1 offset from its anchors' raw MS1 errors; the synthetic run has none, so the
    /// fitted offset is near 0. Calibration's own passes search without an offset.
    /// </summary>
    [Test]
    public void CalibrationFitsTheRunsMs1Offset()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());

        Assert.That(calibration.Ms1Offset, Is.Not.Null);
        Assert.That(Math.Abs(calibration.Ms1Offset!.OffsetPpm(5)), Is.LessThan(1));
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(Ms1Offset: Ms1OffsetModel.Constant(3))).Ms1Offset, Is.Null);
        // Calibration reads MS1 at its own fixed tolerance (the offset is not known yet), whatever the search's
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(Ms1TolerancePpm: 5)).Ms1TolerancePpm,
            Is.EqualTo(DiaIrtSelfCalibration.CalibrationMs1TolerancePpm));
        // Applying a calibration sets the search's iRT window and MS1 offset together
        var applied = calibration.ApplyTo(new DiaLibrarySearchParameters());
        Assert.That(applied.IrtHalfWindow, Is.EqualTo(calibration.IrtHalfWindow));
        Assert.That(applied.Ms1Offset, Is.SameAs(calibration.Ms1Offset));
    }

    /// <summary>The first pass's confident identifications recover the run's true, nonlinear RT(iRT).</summary>
    [Test]
    public void TheFirstPassRecoversTheRunsTrueCurve()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());

        for (double irt = -10; irt <= 110; irt += 10)
        {
            var observed = calibration.Model.ToIrt(new RtMinutes(SyntheticDiaRun.TrueRtMinutes(irt)));
            Assert.That(observed.Value, Is.EqualTo(irt).Within(2.0), $"at library iRT {irt}");
        }
        Assert.That(calibration.IrtHalfWindow, Is.LessThan(20), "a calibrated run needs a narrower window than the first pass");
        Assert.That(calibration.AnchorCount, Is.GreaterThanOrEqualTo(50));
        // The anchors are kept, so a bad fit can be traced to the identifications it was fitted on
        Assert.That(calibration.Anchors, Has.Count.EqualTo(calibration.AnchorCount));
        Assert.That(calibration.Anchors.Count(a => Math.Abs(calibration.Model.ToIrt(a.ApexRt).Value - a.LibraryIrt.Value) < 3), Is.GreaterThan(0.9 * calibration.AnchorCount));
    }

    /// <summary>Two passes find what one pass with a borrowed map found, with no map supplied at all.</summary>
    [Test]
    public void ASelfCalibratedSearchFindsThePlantedPrecursors()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, calibration.Model,
            new DiaLibrarySearchParameters(IrtHalfWindow: calibration.IrtHalfWindow), new CommonParameters(), [], []).Run();

        int found = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
    }

    /// <summary>
    /// The first pass needs only enough targets to yield its anchors. A large library is sampled down to about
    /// <see cref="DiaIrtSelfCalibration.FirstPassTargetCount"/> targets; a small one is searched whole.
    /// </summary>
    [TestCase(300, 1)]
    [TestCase(10_000, 1)]
    [TestCase(10_001, 2)]
    [TestCase(38_987, 4)]
    public void TheFirstPassSamplesALargeLibrary(int targets, int stride)
    {
        Assert.That(DiaIrtSelfCalibration.FirstPassStride(targets), Is.EqualTo(stride));
    }

    /// <summary>
    /// A whole-proteome library is sparse: few of its precursors are in any run. A sample then holds too few for any to
    /// reach 1% q, since (D+1)/T needs a hundred targets before the first decoy. The first pass samples twice as many
    /// precursors until it finds enough anchors, up to the whole library.
    /// </summary>
    [Test]
    public void ASparseLibraryIsSampledMoreUntilTheFirstPassFindsAnchors()
    {
        // 3000 targets of which about 5% elute: a 500-target sample holds about 25
        var run = SyntheticDiaRun.Build(3000, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 20) == 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(),
            firstPassTargetCount: 500);

        Assert.That(calibration.AnchorCount, Is.GreaterThanOrEqualTo(100));
        for (double irt = -10; irt <= 110; irt += 20)
            Assert.That(calibration.Model.ToIrt(new RtMinutes(SyntheticDiaRun.TrueRtMinutes(irt))).Value, Is.EqualTo(irt).Within(2.0), $"at library iRT {irt}");
    }

    /// <summary>Calibration can be repeated on a disjoint first-pass sample; each sample calibrates the run on its own.</summary>
    [TestCase(1)]
    [TestCase(3)]
    public void AnotherFirstPassSampleCalibratesToo(int offset)
    {
        var run = SyntheticDiaRun.Build(3000, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) == 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(),
            firstPassTargetCount: 500, sampleOffset: offset);

        Assert.That(calibration.AnchorCount, Is.GreaterThanOrEqualTo(100));
        for (double irt = -10; irt <= 110; irt += 20)
            Assert.That(calibration.Model.ToIrt(new RtMinutes(SyntheticDiaRun.TrueRtMinutes(irt))).Value, Is.EqualTo(irt).Within(2.0), $"at library iRT {irt}");
        Assert.Throws<ArgumentOutOfRangeException>(() => DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(),
            new CommonParameters(), sampleOffset: -1));
    }

    /// <summary>
    /// The first pass widens until it has the anchors asked for. A run whose first pass stopped at 843 anchors needed a second
    /// search to fix its RT model (PXD022589, +18%); asking for more anchors up front is the cheaper fix to try.
    /// </summary>
    [Test]
    public void AskingForMoreAnchorsWidensTheFirstPassFurther()
    {
        var run = SyntheticDiaRun.Build(3000, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) == 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var few = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(),
            firstPassTargetCount: 500, desiredAnchors: 100);
        var many = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(),
            firstPassTargetCount: 500, desiredAnchors: 600);

        Assert.That(few.AnchorCount, Is.GreaterThanOrEqualTo(100).And.LessThan(600));
        Assert.That(many.AnchorCount, Is.GreaterThanOrEqualTo(600));
    }

    /// <summary>
    /// The search window is a multiple of the calibration's residual SD (at least 5 iRT). On PXD005573, every change that
    /// narrowed it cost identifications, so the multiple is a parameter to measure rather than a constant.
    /// </summary>
    [TestCase(4.0)]
    [TestCase(5.0)]
    public void TheWindowIsTheChosenMultipleOfTheResidualSd(double sds)
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(), windowSds: sds);

        Assert.That(calibration.IrtHalfWindow, Is.EqualTo(Math.Max(5, sds * calibration.Model.ResidualSd)).Within(1e-9));
    }

    /// <summary>
    /// Calibration passes pick anchors with the linear discriminant, whatever model the main search uses: anchors only need to
    /// be confidently right. On the whole-proteome library this cut calibration from 12:18 to 3:51 (PXD005573 1 h) and gave
    /// slightly more precursors at matched FDP.
    /// </summary>
    [Test]
    public void CalibrationUsesTheLinearModelWhateverTheSearchUses()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var network = DiaIrtSelfCalibration.Calibrate(run.Scans, library,
            new DiaLibrarySearchParameters(ClassifierModel: StatisticalModels.RescoreModel.NeuralNetworkEnsemble), new CommonParameters());
        var linear = DiaIrtSelfCalibration.Calibrate(run.Scans, library,
            new DiaLibrarySearchParameters(ClassifierModel: StatisticalModels.RescoreModel.LinearDiscriminant), new CommonParameters());

        Assert.That(network.AnchorCount, Is.EqualTo(linear.AnchorCount));
        Assert.That(network.IrtHalfWindow, Is.EqualTo(linear.IrtHalfWindow));
        Assert.That(network.Model.ResidualSd, Is.EqualTo(linear.Model.ResidualSd));
    }

    /// <summary>
    /// Iterative calibration (DIA-NN calibrates in rounds): each later round searches the same sample again with the previous
    /// round's fitted curve and window instead of the provisional straight line, so anchors come from the whole gradient
    /// rather than where the straight line happened to be right. The curve stays accurate, and the window is the last round's.
    /// </summary>
    [Test]
    public void ALaterCalibrationRoundSearchesWithTheFittedCurve()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var one = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(), rounds: 1);
        var two = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(), rounds: 2);

        Assert.That(two.Rounds, Is.EqualTo(2));
        Assert.That(one.Rounds, Is.EqualTo(1));
        Assert.That(two.AnchorCount, Is.GreaterThanOrEqualTo(one.AnchorCount * 0.9));
        for (double irt = -10; irt <= 110; irt += 10)
            Assert.That(two.Model.ToIrt(new RtMinutes(SyntheticDiaRun.TrueRtMinutes(irt))).Value, Is.EqualTo(irt).Within(2.0), $"at library iRT {irt}");
    }

    /// <summary>
    /// The second pass (DIA-NN refits after its first search): the main search's confident targets refit the run's RT->iRT map.
    /// Only targets at 1% count. A far-off target above 1%, or any decoy, must not pull the fit.
    /// </summary>
    [Test]
    public void ARefitFromTheMainSearchUsesOnlyItsConfidentTargets()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, calibration.Model,
            new DiaLibrarySearchParameters(IrtHalfWindow: calibration.IrtHalfWindow), new CommonParameters(), [], []).Run();
        int confident = results.Matches.Count(m => !m.IsDecoy && m.QValue <= DiaIrtSelfCalibration.AnchorQValue);
        // Poison: wrong apexes far from the curve, as a weak target and as a decoy
        var poisoned = results.Matches.Concat(Enumerable.Range(0, 50).SelectMany(i => new[]
        {
            new DiaPrecursorMatch(100_000 + i, "POISONK", 2, 500, false, new Irt(100), new RtMinutes(1), new Irt(0), 1, [1.0], QValue: 0.5),
            new DiaPrecursorMatch(200_000 + i, "KNOSIOP", 2, 500, true, new Irt(100), new RtMinutes(1), new Irt(0), 99, [1.0]),
        })).ToList();

        var refined = DiaIrtSelfCalibration.Refine(poisoned);

        Assert.That(refined.AnchorCount, Is.EqualTo(Math.Min(confident, DiaIrtSelfCalibration.MaximumRefinementAnchors)));
        for (double irt = -10; irt <= 110; irt += 10)
            Assert.That(refined.Model.ToIrt(new RtMinutes(SyntheticDiaRun.TrueRtMinutes(irt))).Value, Is.EqualTo(irt).Within(2.0), $"at library iRT {irt}");
    }

    /// <summary>With nothing to find there is nothing to calibrate on, and that is an error, not a guess.</summary>
    [Test]
    public void ARunWithNothingToFindCannotBeCalibrated()
    {
        var run = SyntheticDiaRun.Build(200, _ => false);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var e = Assert.Throws<MetaMorpheusException>(() =>
            DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters()));
        Assert.That(e!.Message, Does.Contain("calibrat"));
    }
}
