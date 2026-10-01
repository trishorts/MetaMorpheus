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
