#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using Readers.SpectralLibrary;

namespace Test.DiaLibrarySearch;

/// <summary>
/// Mass calibration for DIA, by MetaMorpheus's own CalibrationEngine (the DDA calibration task's): the run's confident
/// identifications supply labelled points (fragments at the apex scan, the precursor's M0 in the nearest MS1 scan), and the
/// engine corrects every scan's m/z by its smoothed local error. Instruments other than PXD005573's drift by several ppm.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaMassCalibrationTests
{
    private string _directory = "";

    [OneTimeSetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DiaMassCalibrationTests");
        Directory.CreateDirectory(_directory);
    }

    [OneTimeTearDown]
    public void TearDown() => Directory.Delete(_directory, true);

    private static double MedianPlantedFragmentPpm(IEnumerable<MsDataScan> scans, SyntheticDiaRun run)
    {
        var planted = run.Library.Where(e => run.PlantedSequences.Contains(e.FullSequence)).ToList();
        var errors = new List<double>();
        foreach (var scan in scans.Where(s => s.MsnOrder == 2))
            foreach (var entry in planted.Where(e => e.PrecursorMz >= scan.IsolationRange.Minimum && e.PrecursorMz < scan.IsolationRange.Maximum
                         && Math.Abs(scan.RetentionTime - SyntheticDiaRun.TrueRtMinutes(e.RetentionTime)) < 0.02))
                foreach (var f in entry.MatchedFragmentIons)
                {
                    int? i = scan.MassSpectrum.GetClosestPeakIndex(f.Mz);
                    if (i is int k)
                    {
                        double ppm = (scan.MassSpectrum.XArray[k] - f.Mz) / f.Mz * 1e6;
                        if (Math.Abs(ppm) < 20)
                            errors.Add(ppm);
                    }
                }
        return errors.Order().ElementAt(errors.Count / 2);
    }

    [Test]
    public void ARunOffsetBy8PpmIsCalibratedFromItsOwnIdentifications()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, withMs1: true, ppmOffset: 8);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var rt = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, rt.Model,
            new DiaLibrarySearchParameters(IrtHalfWindow: rt.IrtHalfWindow), new CommonParameters(), [], []).Run();
        Assert.That(MedianPlantedFragmentPpm(run.Scans, run), Is.EqualTo(8).Within(0.5), "the fixture's offset");

        var calibration = DiaMassCalibration.Calibrate(run.Scans, results.Matches, library, new CommonParameters());

        Assert.That(calibration.Ms2Points, Is.GreaterThan(100));
        Assert.That(calibration.Ms1Points, Is.GreaterThan(50));
        Assert.That(calibration.Ms2MedianPpmBefore, Is.EqualTo(8).Within(1));
        Assert.That(calibration.Ms1MedianPpmBefore, Is.EqualTo(8).Within(1));
        Assert.That(MedianPlantedFragmentPpm(calibration.Scans, run), Is.EqualTo(0).Within(1), "corrected");
        Assert.That(calibration.Scans, Has.Length.EqualTo(run.Scans.Length));
    }

    /// <summary>With no confident identifications there is nothing to calibrate on: the scans come back unchanged.</summary>
    [Test]
    public void WithoutIdentificationsTheScansAreReturnedUnchanged()
    {
        var run = SyntheticDiaRun.Build(50, _ => false, withMs1: true, ppmOffset: 8);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaMassCalibration.Calibrate(run.Scans, [], library, new CommonParameters());

        Assert.That(calibration.Scans, Is.SameAs(run.Scans));
        Assert.That(calibration.Ms2Points, Is.EqualTo(0));
    }
}
