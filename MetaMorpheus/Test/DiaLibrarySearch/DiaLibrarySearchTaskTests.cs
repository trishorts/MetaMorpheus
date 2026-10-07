#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using TaskLayer;

namespace Test.DiaLibrarySearch;

/// <summary>
/// M7 slice 1 (design/M7.md): the DIA library search as a MetaMorpheus task. It reads and writes TOML like every task,
/// takes the library as a database (`-d library.msl`) with no FASTA, loads DIA scans with every peak, and reports the
/// headline count PXReprise reads from results.txt.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaLibrarySearchTaskTests
{
    private string _directory = "";

    [OneTimeSetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DiaLibrarySearchTaskTests");
        Directory.CreateDirectory(_directory);
    }

    [OneTimeTearDown]
    public void TearDown()
    {
        Directory.Delete(_directory, true);
    }

    /// <summary>Writes a synthetic run as mzML and its library as .msl, as a user would hand them to the task.</summary>
    private (string Mzml, string Library, SyntheticDiaRun Run) WriteInputs(string name, double noisePeaksPerScan = 60,
        Func<Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry, double>? abundance = null)
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, withMs1: true, noisePeaksPerScan: noisePeaksPerScan,
            abundance: abundance);
        string folder = Path.Combine(_directory, name);
        Directory.CreateDirectory(folder);
        string mzml = Path.Combine(folder, name + ".mzML");
        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(AsWritten(run.Scans), new SourceFile("no nativeID format", "mzML format", null, null, null)), mzml, false);
        return (mzml, run.WriteLibrary(folder), run);
    }

    /// <summary>
    /// The scans as a converter writes a DIA run: each MS2 names the window centre as its selected ion and the cycle's MS1
    /// as its precursor scan (the mzML writer requires a selected ion).
    /// </summary>
    private static MsDataScan[] AsWritten(MsDataScan[] scans)
    {
        int? lastMs1 = null;
        return scans.Select(scan =>
        {
            if (scan.MsnOrder == 1)
            {
                lastMs1 = scan.OneBasedScanNumber;
                return scan;
            }
            return new MsDataScan(scan.MassSpectrum, scan.OneBasedScanNumber, scan.MsnOrder, scan.IsCentroid, scan.Polarity, scan.RetentionTime,
                scan.ScanWindowRange, scan.ScanFilter, scan.MzAnalyzer, scan.TotalIonCurrent, scan.InjectionTime, scan.NoiseData, scan.NativeId,
                selectedIonMz: scan.IsolationMz, isolationMZ: scan.IsolationMz, isolationWidth: scan.IsolationWidth,
                dissociationType: scan.DissociationType, oneBasedPrecursorScanNumber: lastMs1);
        }).ToArray();
    }

    /// <summary>The task's settings survive a TOML write and read, and the default task is a DIA library search.</summary>
    [Test]
    public void TheTaskRoundTripsThroughToml()
    {
        var task = new DiaLibrarySearchTask(new DiaLibrarySearchTaskParameters { FragmentTolerancePpm = 15, Ms1TolerancePpm = 4, ClassifierSeed = 3, QValueThreshold = 0.05 });
        string path = Path.Combine(_directory, "DiaLibrarySearchTask.toml");
        Toml.WriteFile(task, path, MetaMorpheusTask.tomlConfig);

        var read = Toml.ReadFile<DiaLibrarySearchTask>(path, MetaMorpheusTask.tomlConfig);

        Assert.That(read.TaskType, Is.EqualTo(MyTask.DiaLibrarySearch));
        Assert.That(read.Parameters.FragmentTolerancePpm, Is.EqualTo(15));
        Assert.That(read.Parameters.Ms1TolerancePpm, Is.EqualTo(4));
        Assert.That(read.Parameters.ClassifierSeed, Is.EqualTo(3));
        Assert.That(read.Parameters.QValueThreshold, Is.EqualTo(0.05));
        Assert.That(new DiaLibrarySearchTask().Parameters.QValueThreshold, Is.EqualTo(0.01));
        Assert.That(new DiaLibrarySearchTask().CommonParameters, Is.Not.Null);
    }

    /// <summary>
    /// The engine's parameters come from the task's, with everything a user cannot set left at the engine's (benchmarked)
    /// defaults.
    /// </summary>
    [Test]
    public void TheEngineParametersFollowTheTasks()
    {
        var engine = new DiaLibrarySearchTaskParameters { FragmentTolerancePpm = 15, Ms1TolerancePpm = 4, ClassifierSeed = 3 }.ToEngineParameters();
        var defaults = new EngineLayer.DiaLibrarySearch.DiaLibrarySearchParameters();

        Assert.That(engine.FragmentTolerancePpm, Is.EqualTo(15));
        Assert.That(engine.Ms1TolerancePpm, Is.EqualTo(4));
        Assert.That(engine.ClassifierSeed, Is.EqualTo(3));
        Assert.That(engine with { FragmentTolerancePpm = defaults.FragmentTolerancePpm, Ms1TolerancePpm = defaults.Ms1TolerancePpm, ClassifierSeed = 0 },
            Is.EqualTo(defaults));
        Assert.That(new DiaLibrarySearchTaskParameters().ToEngineParameters(), Is.EqualTo(defaults));
    }

    /// <summary>
    /// No peak trimming or filtering on DIA input (user, 2026-09-30): the task keeps every MS1 and MS2 peak, although
    /// CommonParameters' default trims MS2 to the top 200 peaks per window, which would bite on this run.
    /// </summary>
    [Test]
    public void DiaScansAreLoadedWithEveryPeak()
    {
        var (mzml, _, _) = WriteInputs("untrimmed", noisePeaksPerScan: 3000);
        var all = MsDataFileReader.GetDataFile(mzml).LoadAllStaticData().GetAllScansList();

        var loaded = DiaLibrarySearchTask.LoadScans(mzml);

        Assert.That(loaded.Select(s => s.MassSpectrum.Size), Is.EqualTo(all.Select(s => s.MassSpectrum.Size)));
        var trimmed = new MyFileManager(true).LoadFile(mzml, new CommonParameters()).GetAllScansList();
        Assert.That(trimmed.Sum(s => s.MassSpectrum.Size), Is.LessThan(all.Sum(s => s.MassSpectrum.Size)), "the default filter would have trimmed");
    }

    /// <summary>
    /// The task searches a run against a library with no FASTA, end to end, and results.txt carries the headline PXReprise
    /// reads: the target precursors at q &lt;= 0.01, per run and in total. Most planted targets are found and few others.
    /// </summary>
    [Test]
    public void TheTaskSearchesARunAgainstALibraryEndToEnd()
    {
        var (count, results, run, matches) = RunTheTask("endtoend");

        Assert.That(count, Is.EqualTo(DiaLibrarySearchTask.CountPassingTargets(matches, 0.01)), "the headline is the engine's count");
        int planted = run.Library.Count(e => !e.IsDecoy && run.PlantedSequences.Contains(e.FullSequence));
        Assert.That(count, Is.InRange(planted * 8 / 10, planted * 11 / 10));
        Assert.That(results, Does.Contain("endtoend.mzML: " + count + " target precursors with q-value <= 0.01"));
    }

    /// <summary>
    /// DiaIrtCalibration.tsv (D12): each run's calibration anchors, so a bad fit can be seen. One row per anchor: the run,
    /// its apex in minutes, its library iRT, the calibrated iRT at that apex and the residual; as many rows as anchors.
    /// </summary>
    [Test]
    public void TheTaskWritesEachRunsCalibration()
    {
        var (mzml, library, _) = WriteInputs("calibration");
        string output = Path.Combine(_directory, "calibration", "output", "Task1");
        Directory.CreateDirectory(output);
        new DiaLibrarySearchTask().RunTask(output, [new DbForTask(library, false)], [mzml], "Task1");

        var lines = File.ReadAllLines(Path.Combine(output, "DiaIrtCalibration.tsv"));
        Assert.That(lines[0].Split('\t'), Is.EqualTo(new[] { "File Name", "AnchorRtMin", "LibraryIrt", "CalibratedIrt", "ResidualIrt" }));

        using var msl = Readers.SpectralLibrary.MslLibrary.Load(library);
        var calibration = EngineLayer.DiaLibrarySearch.DiaIrtSelfCalibration.Calibrate(DiaLibrarySearchTask.LoadScans(mzml), msl,
            new DiaLibrarySearchTaskParameters().ToEngineParameters(), new CommonParameters());
        Assert.That(lines.Length - 1, Is.EqualTo(calibration.Anchors.Count).And.GreaterThan(0));
        for (int i = 0; i < calibration.Anchors.Count; i++)
        {
            var cells = lines[i + 1].Split('\t');
            var (rt, irt) = calibration.Anchors[i];
            double calibrated = calibration.Model.ToIrt(rt).Value;
            Assert.That(cells[0], Is.EqualTo("calibration.mzML"));
            Assert.That(double.Parse(cells[1], System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(rt.Value));
            Assert.That(double.Parse(cells[2], System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(irt.Value));
            Assert.That(double.Parse(cells[3], System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(calibrated));
            Assert.That(double.Parse(cells[4], System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(irt.Value - calibrated).Within(1e-9));
        }
    }

    /// <summary>
    /// AllDiaPrecursors.tsv (M7 slice 2) is written through mzLib's DiaPrecursorFile, in the schema agreed with dataRepo: one
    /// row per target precursor at the threshold, so its rows are the headline; the run key verbatim; each row's apex scan,
    /// retention in minutes and iRT, and quantity; and, with no contaminant database, "contaminants: not assessed".
    /// </summary>
    [Test]
    public void TheTaskWritesThePrecursorTable()
    {
        var (count, _, run, matches) = RunTheTask("table");
        string path = Path.Combine(_directory, "table", "output", "Task1", "AllDiaPrecursors.tsv");

        var table = new DiaPrecursorFile(path);
        table.LoadResults();

        Assert.That(table.ContaminantsAssessed, Is.False);
        Assert.That(table.Results.Count, Is.EqualTo(count));
        var bySequence = matches.Where(m => !m.IsDecoy).ToDictionary(m => (m.FullSequence, m.Charge));
        foreach (var row in table.Results)
        {
            var match = bySequence[(row.FullSequence, row.PrecursorCharge)];
            Assert.That(row.FileName, Is.EqualTo("table.mzML"));
            Assert.That(row.Label, Is.EqualTo("T"));
            Assert.That(row.QValuePrecursorRun, Is.EqualTo(match.QValue).And.LessThanOrEqualTo(0.01));
            Assert.That(row.QValuePrecursorGlobal, Is.EqualTo(row.QValuePrecursorRun), "one run: global is the run's");
            Assert.That(row.ApexScanNumber, Is.EqualTo(match.ApexScanNumber).And.GreaterThan(0));
            Assert.That(row.ApexRtMin, Is.EqualTo(match.ApexRt.Value));
            Assert.That(row.ApexIrt, Is.EqualTo(match.ApexIrt.Value));
            Assert.That(row.LibraryIrt, Is.EqualTo(match.LibraryIrt.Value));
            Assert.That(row.PrecursorMz, Is.EqualTo(match.PrecursorMz));
            Assert.That(row.Score, Is.EqualTo(match.Score));
            Assert.That(row.PrecursorQuantity, Is.EqualTo(double.IsNaN(match.Quantity) ? null : match.Quantity));
            var entry = run.Library.Single(e => e.FullSequence == row.FullSequence && e.ChargeState == row.PrecursorCharge);
            Assert.That(row.BaseSequence, Is.EqualTo(entry.BaseSequence));
            Assert.That(row.ProteinAccession, Is.EqualTo(entry.ProteinAccession ?? ""), "no accession is an empty cell");
        }
    }

    /// <summary>
    /// The headline counts targets at the q-value threshold, not every target scored, and never a decoy (whose q is NaN).
    /// </summary>
    [Test]
    public void TheHeadlineCountsOnlyTargetsPassingTheThreshold()
    {
        static EngineLayer.DiaLibrarySearch.DiaPrecursorMatch Match(bool decoy, double q) =>
            new(0, "PEPTIDE", 2, 400, decoy, new MzLibUtil.Irt(0), new MzLibUtil.RtMinutes(1), new MzLibUtil.Irt(0), 1, [], q);

        var matches = new[] { Match(false, 0.001), Match(false, 0.01), Match(false, 0.02), Match(true, double.NaN), Match(true, 0.0) };

        Assert.That(DiaLibrarySearchTask.CountPassingTargets(matches, 0.01), Is.EqualTo(2));
        Assert.That(DiaLibrarySearchTask.CountPassingTargets(matches, 0.05), Is.EqualTo(3));
        Assert.That(DiaLibrarySearchTask.PassingTargets(matches, 0.01), Is.EqualTo(new[] { matches[0], matches[1] }), "the table's rows");
    }

    /// <summary>
    /// Runs the task on a synthetic run, and the same pipeline by hand on the same files. Returns the headline count,
    /// results.txt, the run, and the hand-run engine's matches.
    /// </summary>
    private (int Count, string Results, SyntheticDiaRun Run, List<EngineLayer.DiaLibrarySearch.DiaPrecursorMatch> Matches) RunTheTask(string name,
        Func<Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry, double>? abundance = null)
    {
        var (mzml, library, run) = WriteInputs(name, abundance: abundance);
        string output = Path.Combine(_directory, name, "output", "Task1");
        Directory.CreateDirectory(output);

        new DiaLibrarySearchTask().RunTask(output, [new DbForTask(library, false)], [mzml], "Task1");

        string results = File.ReadAllText(Path.Combine(output, "results.txt"));
        var headline = Regex.Match(results, @"All target precursors with q-value <= 0\.01: (\d+)");
        Assert.That(headline.Success, results);

        var scans = DiaLibrarySearchTask.LoadScans(mzml);
        using var msl = Readers.SpectralLibrary.MslLibrary.Load(library);
        var parameters = new DiaLibrarySearchTaskParameters().ToEngineParameters();
        var calibration = EngineLayer.DiaLibrarySearch.DiaIrtSelfCalibration.Calibrate(scans, msl, parameters, new CommonParameters());
        var matches = ((EngineLayer.DiaLibrarySearch.DiaLibrarySearchResults)new EngineLayer.DiaLibrarySearch.DiaLibrarySearchEngine(scans, msl,
            calibration.Model, calibration.ApplyTo(parameters), new CommonParameters(), [], []).Run()).Matches;
        return (int.Parse(headline.Groups[1].Value), results, run, matches);
    }

    /// <summary>
    /// The runner refuses a search with no protein database, but a DIA library search needs only its library: a .msl alone
    /// must reach the task.
    /// </summary>
    [Test]
    public void TheRunnerAcceptsALibraryAloneForTheDiaTask()
    {
        var (mzml, library, _) = WriteInputs("runner");
        string output = Path.Combine(_directory, "runner", "output");
        var warnings = new List<string>();
        EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
        EverythingRunnerEngine.WarnHandler += onWarn;
        try
        {
            new EverythingRunnerEngine([("Task1", new DiaLibrarySearchTask())], [mzml], [new DbForTask(library, false)], output).Run();
        }
        finally
        {
            EverythingRunnerEngine.WarnHandler -= onWarn;
        }

        Assert.That(warnings, Has.None.Contains("No protein database"));
        Assert.That(File.Exists(Path.Combine(output, "Task1", "results.txt")));
    }
}
