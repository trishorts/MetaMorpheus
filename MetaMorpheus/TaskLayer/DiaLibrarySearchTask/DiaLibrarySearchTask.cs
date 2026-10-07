#nullable enable
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using Readers;
using Readers.SpectralLibrary;

namespace TaskLayer;

/// <summary>
/// Searches DIA runs against one spectral library (.msl), given as a database. Each run calibrates its own retention time
/// onto the library's iRT and its MS1 offset, then is searched and scored; results.txt reports the target precursors at
/// the q-value threshold, per run and in total. A protein database is not needed.
/// </summary>
public class DiaLibrarySearchTask : MetaMorpheusTask
{
    public DiaLibrarySearchTaskParameters Parameters { get; set; }

    public DiaLibrarySearchTask(DiaLibrarySearchTaskParameters parameters) : base(MyTask.DiaLibrarySearch)
    {
        CommonParameters = new CommonParameters();
        Parameters = parameters;
    }

    /// <summary>Used when reading TOML and when `CMD -g` writes the default task.</summary>
    public DiaLibrarySearchTask() : this(new DiaLibrarySearchTaskParameters())
    {
    }

    /// <summary>
    /// A run's scans with every peak: never through <see cref="MyFileManager"/>'s filter, so CommonParameters' MS2
    /// trimming (on by default) can never reach DIA scans (user, 2026-09-30).
    /// </summary>
    public static MsDataScan[] LoadScans(string path) =>
        MsDataFileReader.GetDataFile(path).LoadAllStaticData().GetAllScansList().ToArray();

    /// <summary>
    /// A match as a row of AllDiaPrecursors.tsv (mzLib's <see cref="DiaPrecursorFromTsv"/>, the schema agreed with dataRepo).
    /// The run key is the file name as given; the quantity is empty when not measured.
    /// </summary>
    internal static DiaPrecursorFromTsv PrecursorRow(DiaPrecursorMatch match, string fileName, MslLibrary library, double globalQValue)
    {
        var entry = library.GetEntry(match.PrecursorIndex);
        return new DiaPrecursorFromTsv
        {
            FileName = fileName,
            FullSequence = match.FullSequence,
            BaseSequence = entry?.BaseSequence ?? "",
            PrecursorCharge = match.Charge,
            PrecursorMz = match.PrecursorMz,
            ProteinAccession = entry?.ProteinAccession ?? "",
            Label = match.IsDecoy ? "D" : "T",
            Score = match.Score,
            QValuePrecursorRun = match.QValue,
            QValuePrecursorGlobal = globalQValue,
            LibraryIrt = match.LibraryIrt.Value,
            ApexIrt = match.ApexIrt.Value,
            ApexRtMin = match.ApexRt.Value,
            ApexScanNumber = match.ApexScanNumber,
            PrecursorQuantity = double.IsNaN(match.Quantity) ? null : match.Quantity,
        };
    }

    /// <summary>The target precursors at or below the q-value threshold, which the headline counts and the table lists; never a decoy.</summary>
    public static IEnumerable<DiaPrecursorMatch> PassingTargets(IEnumerable<DiaPrecursorMatch> matches, double qValueThreshold) =>
        matches.Where(m => !m.IsDecoy && m.QValue <= qValueThreshold);

    /// <summary>
    /// A run's passing targets that are also reported (PLAN D4): those whose global q-value passes too.
    /// </summary>
    public static List<DiaPrecursorMatch> Reported(IEnumerable<DiaPrecursorMatch> runPassing, IReadOnlyDictionary<int, double> globalQValues,
        double qValueThreshold) =>
        runPassing.Where(m => globalQValues[m.PrecursorIndex] <= qValueThreshold).ToList();

    /// <summary>How many <see cref="PassingTargets"/> there are.</summary>
    public static int CountPassingTargets(IEnumerable<DiaPrecursorMatch> matches, double qValueThreshold) =>
        PassingTargets(matches, qValueThreshold).Count();

    protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId,
        FileSpecificParameters[] fileSettingsList)
    {
        MyTaskResults = new MyTaskResults(this);
        var libraries = dbFilenameList.Where(db => db.IsSpectralLibrary).ToList();
        if (libraries.Count != 1 || !string.Equals(Path.GetExtension(libraries[0].FilePath), ".msl", StringComparison.OrdinalIgnoreCase))
            throw new MetaMorpheusException($"A DIA library search needs exactly one .msl library; {libraries.Count} spectral libraries were given.");

        string libraryPath = libraries[0].FilePath;
        Status("Loading library...", taskId);
        // Above 2 GB the whole-file load cannot hold the library in one array (MSL-Q6); index-only reads entries on demand
        using var library = new FileInfo(libraryPath).Length > int.MaxValue ? MslLibrary.LoadIndexOnly(libraryPath) : MslLibrary.Load(libraryPath);
        var parameters = Parameters.ToEngineParameters();
        string threshold = Parameters.QValueThreshold.ToString(CultureInfo.InvariantCulture);

        // Each run's targets passing its own q-value, kept until every run is in and the global q-values are known
        var passing = new List<(string FileName, List<DiaPrecursorMatch> Matches)>();
        var bestScores = new DiaPrecursorFdr.BestScores();
        var calibrationLines = new List<string> { "File Name\tAnchorRtMin\tLibraryIrt\tCalibratedIrt\tResidualIrt" };
        foreach (string rawFile in currentRawFileList)
        {
            if (GlobalVariables.StopLoops)
                break;
            string fileName = Path.GetFileName(rawFile);
            var ids = new List<string> { taskId, "Individual Spectra Files", rawFile };
            StartingDataFile(rawFile, ids);

            Status("Loading spectra file...", ids);
            var scans = LoadScans(rawFile);
            Status("Calibrating retention time...", ids);
            var calibration = DiaIrtSelfCalibration.Calibrate(scans, library, parameters, CommonParameters);
            foreach (var (rt, irt) in calibration.Anchors)
            {
                double calibrated = calibration.Model.ToIrt(rt).Value;
                calibrationLines.Add(string.Join('\t', fileName, rt.Value.ToString("R", CultureInfo.InvariantCulture), irt.Value.ToString("R", CultureInfo.InvariantCulture),
                    calibrated.ToString("R", CultureInfo.InvariantCulture), (irt.Value - calibrated).ToString("R", CultureInfo.InvariantCulture)));
            }
            Status("Searching...", ids);
            var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(scans, library, calibration.Model, calibration.ApplyTo(parameters),
                CommonParameters, FileSpecificParameters, ids).Run();

            bestScores.Add(results.Matches);
            passing.Add((fileName, PassingTargets(results.Matches, Parameters.QValueThreshold).ToList()));
            FinishedDataFile(rawFile, ids);
        }

        // PLAN D4: a precursor is reported in a run when both its run q-value and its global q-value (best score across the
        // runs, then target-decoy competition) pass. The headline is the distinct precursors reported; with one run the
        // global q-value is the run's, so every run-passing target is reported.
        var global = bestScores.GlobalQValues();
        var rows = new List<DiaPrecursorFromTsv>();
        var reported = new HashSet<int>();
        var perRun = new List<string>();
        foreach (var (fileName, matches) in passing)
        {
            var kept = Reported(matches, global, Parameters.QValueThreshold);
            rows.AddRange(kept.Select(m => PrecursorRow(m, fileName, library, global[m.PrecursorIndex])));
            reported.UnionWith(kept.Select(m => m.PrecursorIndex));
            perRun.Add($"{fileName}: {kept.Count} target precursors with q-value <= {threshold}");
        }
        int total = reported.Count;

        string tablePath = Path.Combine(OutputFolder, "AllDiaPrecursors.tsv");
        new DiaPrecursorFile(tablePath, rows, contaminantsAssessed: dbFilenameList.Any(db => db.IsContaminant)).WriteResults(tablePath);
        FinishedWritingFile(tablePath, new List<string> { taskId });

        // Each run's calibration anchors (D12), so a bad retention-time fit can be seen
        string calibrationPath = Path.Combine(OutputFolder, "DiaIrtCalibration.tsv");
        File.WriteAllLines(calibrationPath, calibrationLines);
        FinishedWritingFile(calibrationPath, new List<string> { taskId });

        MyTaskResults.AddTaskSummaryText($"All target precursors with q-value <= {threshold}: {total}");
        foreach (string line in perRun)
            MyTaskResults.AddTaskSummaryText(line);
        return MyTaskResults;
    }
}
