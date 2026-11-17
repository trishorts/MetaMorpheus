#nullable enable
using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;

namespace Test.DiaLibrarySearch;

/// <summary>
/// A run's fragment mass error, measured on calibration's anchors. DIA-NN sets its MS2 tolerance per run from such a
/// measurement (12 ppm on A1 and PXD005573 1 h, 20 on HF-X), and our fixed 20 ppm costs A1 10% against 10 ppm.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaFragmentMassErrorTests
{
    /// <summary>
    /// -2..+2 ppm around +1.5, twenty times each: the median is 1.5, the median absolute deviation 1 (a robust SD of 1.4826),
    /// and 95% of the deviations are within 2.
    /// </summary>
    [Test]
    public void TheSpreadIsARobustSdAroundTheMedian()
    {
        var errors = Enumerable.Range(0, 20).SelectMany(_ => new[] { -2.0, -1, 0, 1, 2 }).Select(e => e + 1.5).ToList();

        var fitted = DiaFragmentMassError.Fit(errors);

        Assert.That(fitted.MedianPpm, Is.EqualTo(1.5).Within(1e-12));
        Assert.That(fitted.SpreadPpm, Is.EqualTo(1.4826).Within(1e-12));
        Assert.That(fitted.Deviation95Ppm, Is.EqualTo(2).Within(1e-12));
        Assert.That(fitted.Count, Is.EqualTo(100));
    }

    [Test]
    public void ErrorsThatAreNotFiniteAreIgnored()
    {
        var fitted = DiaFragmentMassError.Fit([1, 2, 3, double.NaN, double.PositiveInfinity]);

        Assert.That(fitted.MedianPpm, Is.EqualTo(2));
        Assert.That(fitted.Count, Is.EqualTo(3));
    }

    [Test]
    public void NothingMeasuredIsNaN()
    {
        var fitted = DiaFragmentMassError.Fit([double.NaN]);

        Assert.That(fitted.MedianPpm, Is.NaN);
        Assert.That(fitted.SpreadPpm, Is.NaN);
        Assert.That(fitted.Deviation95Ppm, Is.NaN);
        Assert.That(fitted.Count, Is.EqualTo(0));
    }

    // ---- Reading an anchor ----

    private static MslFragmentIon Ion(int number, float intensity, float mz) =>
        new() { ProductType = ProductType.y, FragmentNumber = number, Intensity = intensity, Mz = mz, Charge = 1 };

    private static MslLibraryEntry Entry(params MslFragmentIon[] ions)
    {
        var entry = new MslLibraryEntry { FullSequence = "PEPTIDEK", BaseSequence = "PEPTIDEK", PrecursorMz = 500, ChargeState = 2 };
        entry.MatchedFragmentIons.AddRange(ions);
        return entry;
    }

    private static double Shifted(double mz, double ppm) => mz * (1 + ppm * 1e-6);

    /// <summary>A window of 21 scans, one a minute: each peak is given per scan as (m/z, height), height 0 for absent.</summary>
    private static MsDataScan[] Window(params Func<int, (double Mz, double Height)>[] peaks) =>
        Enumerable.Range(0, 21).Select(scan =>
        {
            var present = peaks.Select(p => p(scan)).Where(p => p.Height > 0).OrderBy(p => p.Mz).ToArray();
            var spectrum = new MzSpectrum(present.Select(p => p.Mz).ToArray(), present.Select(p => p.Height).ToArray(), false);
            return new MsDataScan(spectrum, scan + 1, 2, true, Polarity.Positive, scan, new MzRange(150, 1500), "synthetic", MZAnalyzerType.Orbitrap,
                present.Sum(p => p.Height), 10, null, $"scan={scan + 1}", isolationMZ: 500, isolationWidth: 25, dissociationType: DissociationType.HCD);
        }).ToArray();

    /// <summary>
    /// The three strongest library fragments are read at the apex scan alone, each at the peak the search would pick; a weaker
    /// fragment is not read, however well it is seen, and the scans either side of the apex do not count.
    /// </summary>
    [Test]
    public void AnAnchorsStrongestFragmentsAreReadAtItsApex()
    {
        var y3 = Ion(3, 1.0f, 300);
        var y4 = Ion(4, 0.9f, 400);
        var y5 = Ion(5, 0.8f, 500);
        var weak = Ion(6, 0.1f, 600);
        var window = Window(
            s => (Shifted(300, s == 10 ? 5 : 12), 1000),
            s => (Shifted(400, s == 10 ? 5 : -12), 1000),
            s => (Shifted(500, -3), 1000),
            s => (Shifted(600, 15), 1000));

        var errors = DiaFragmentMassError.ReadAnchor(Entry(weak, y5, y4, y3), new RtMinutes(10), [window], 20, true);

        Assert.That(errors.Order(), Is.EqualTo(new[] { -3.0, 5, 5 }).Within(1e-3));
    }

    /// <summary>A strong fragment with no peak in tolerance gives nothing, and no weaker fragment takes its place.</summary>
    [Test]
    public void AStrongFragmentNotSeenGivesNothing()
    {
        var window = Window(s => (Shifted(300, 2), 1000), s => (Shifted(400, 25), 1000), s => (Shifted(500, -1), 1000), s => (Shifted(600, 1), 1000));

        var errors = DiaFragmentMassError.ReadAnchor(Entry(Ion(3, 1.0f, 300), Ion(4, 0.9f, 400), Ion(5, 0.8f, 500), Ion(6, 0.1f, 600)), new RtMinutes(10), [window], 20, true);

        Assert.That(errors.Order(), Is.EqualTo(new[] { -1.0, 2 }).Within(1e-3), "y4 is 25 ppm off");
    }

    /// <summary>Of two windows holding the precursor, the one whose scan comes nearest the apex is read, as for the intensity table.</summary>
    [Test]
    public void TheWindowNearestTheApexIsRead()
    {
        var near = Window(s => (Shifted(300, 4), 1000));
        var far = Window(s => (Shifted(300, -8), 1000)).Select(s => new MsDataScan(s.MassSpectrum, s.OneBasedScanNumber + 100, 2, true, Polarity.Positive, s.RetentionTime + 0.5,
            new MzRange(150, 1500), "synthetic", MZAnalyzerType.Orbitrap, s.TotalIonCurrent, 10, null, $"scan={s.OneBasedScanNumber + 100}", isolationMZ: 500, isolationWidth: 25,
            dissociationType: DissociationType.HCD)).ToArray();

        var errors = DiaFragmentMassError.ReadAnchor(Entry(Ion(3, 1.0f, 300)), new RtMinutes(10), [far, near], 20, true);

        Assert.That(errors, Is.EqualTo(new[] { 4.0 }).Within(1e-3));
    }

    [Test]
    public void AnAnchorWithNoWindowReadsNothing()
    {
        Assert.That(DiaFragmentMassError.ReadAnchor(Entry(Ion(3, 1.0f, 300)), new RtMinutes(10), [], 20, true), Is.Empty);
    }

    // Calibration measuring a run: DiaIrtSelfCalibrationTests.CalibrationMeasuresTheRunsFragmentMassError
}
