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
    internal static DiaPrecursorFromTsv PrecursorRow(DiaPrecursorMatch match, string fileName, MslLibrary library, bool singleRun)
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
            QValuePrecursorGlobal = singleRun ? match.QValue : double.NaN,
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

        int total = 0;
        var perRun = new List<string>();
        var rows = new List<DiaPrecursorFromTsv>();
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
            Status("Searching...", ids);
            var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(scans, library, calibration.Model, calibration.ApplyTo(parameters),
                CommonParameters, FileSpecificParameters, ids).Run();

            int count = CountPassingTargets(results.Matches, Parameters.QValueThreshold);
            total += count;
            perRun.Add($"{fileName}: {count} target precursors with q-value <= {threshold}");
            rows.AddRange(PassingTargets(results.Matches, Parameters.QValueThreshold)
                .Select(m => PrecursorRow(m, fileName, library, singleRun: currentRawFileList.Count == 1)));
            FinishedDataFile(rawFile, ids);
        }

        // One row per target precursor at the threshold, so the table's rows are the headline. Global q-values across runs
        // come with several-run searches (design/M7.md slice 3); until then they are written only for a single run.
        string tablePath = Path.Combine(OutputFolder, "AllDiaPrecursors.tsv");
        new DiaPrecursorFile(tablePath, rows, contaminantsAssessed: dbFilenameList.Any(db => db.IsContaminant)).WriteResults(tablePath);
        FinishedWritingFile(tablePath, new List<string> { taskId });

        MyTaskResults.AddTaskSummaryText($"All target precursors with q-value <= {threshold}: {total}");
        foreach (string line in perRun)
            MyTaskResults.AddTaskSummaryText(line);
        return MyTaskResults;
    }
}
