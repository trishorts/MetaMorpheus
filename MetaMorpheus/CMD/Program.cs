using CommandLine;
using CommandLine.Text;
using EngineLayer;
using IO.ThermoRawFileReader;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using EngineLayer.DatabaseLoading;
using Omics.Modifications;
using TaskLayer;
// Aliased rather than imported: `using Readers` collides with IO.ThermoRawFileReader over
// ThermoRawFileReaderLicence, which this file uses, and CS0104 is not worth a whole namespace.
using SdrfQuantAudit = Readers.SdrfQuantAudit;
using SdrfQuantAuditor = Readers.SdrfQuantAuditor;
using SdrfSearchDesign = Readers.SdrfSearchDesign;
using SdrfSearchDesignOptions = Readers.SdrfSearchDesignOptions;
using MzLibIsobaricMassTag = Omics.Modifications.IsobaricMassTag;

namespace MetaMorpheusCommandLine
{
    public static class Program
    {
        private static bool InProgress;
        private static CommandLineSettings CommandLineSettings;

        private static System.CodeDom.Compiler.IndentedTextWriter MyWriter = new System.CodeDom.Compiler.IndentedTextWriter(Console.Out, "\t");

        public static int Main(string[] args)
        {
            // an error code of 0 is returned if the program ran successfully.
            // otherwise, an error code of >0 is returned.
            // this makes it easier to determine via scripts when the program fails.
            int errorCode = 0;

            var parser = new Parser(with => with.HelpWriter = null);
            var parserResult = parser.ParseArguments<CommandLineSettings>(args);

            parserResult
              .WithParsed<CommandLineSettings>(options => errorCode = Run(options))
              .WithNotParsed(errs => errorCode = DisplayHelp(parserResult, errs));

            return errorCode;
        }

        public static int DisplayHelp<T>(ParserResult<T> result, IEnumerable<Error> errs)
        {
            Console.WriteLine("Welcome to MetaMorpheus");
            Console.WriteLine(GlobalVariables.MetaMorpheusVersion);

            int errorCode = 0;

            var helpText = HelpText.AutoBuild(result, h =>
            {
                h.AdditionalNewLineAfterOption = false;
                h.Copyright = "";
                return HelpText.DefaultParsingErrorsHandler(result, h);
            }, e => e);

            helpText.MaximumDisplayWidth = 300;

            helpText.AddPostOptionsLine("Example usage (Windows): ");
            helpText.AddPostOptionsLine("CMD.exe -d C:\\ExampleDatabase.fasta -s C:\\ExampleSpectra.mzML -t C:\\ExampleTask.toml");
            helpText.AddPostOptionsLine(Environment.NewLine);

            helpText.AddPostOptionsLine("Example usage (Linux): ");
            helpText.AddPostOptionsLine("dotnet CMD.dll -d home/mydata/ExampleDatabase.fasta -s home/mydata/ExampleSpectra.mzML -t home/mydata/ExampleTask.toml");
            helpText.AddPostOptionsLine(Environment.NewLine);

            Console.WriteLine(helpText);

            if (errs.Any(x => x.Tag != ErrorType.HelpRequestedError))
            {
                errorCode = 1;
            }

            return errorCode;
        }

        private static int Run(CommandLineSettings settings)
        {
            int errorCode = 0;

            if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                Console.WriteLine("Welcome to MetaMorpheus");
            }

            if (settings.CustomDataDirectory != null)
            {
                GlobalVariables.UserSpecifiedDataDir = settings.CustomDataDirectory;
            }

            // subscribed before startup so that what startup notices, such as a custom protease that
            // collided with a built-in, reaches the console; WarnHandler reads the verbosity from here.
            // Removed first so that calling Run more than once in a process does not print it twice.
            CommandLineSettings = settings;
            GlobalVariables.WarnHandler -= WarnHandler;
            GlobalVariables.WarnHandler += WarnHandler;
            GlobalVariables.SetUpGlobalVariables();

            if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                Console.WriteLine(GlobalVariables.MetaMorpheusVersion);
            }

            try
            {
                settings.ValidateCommandLineSettings();
                CommandLineSettings = settings;
            }
            catch (Exception e)
            {
                if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    Console.WriteLine("MetaMorpheus encountered the following error:" + Environment.NewLine + e.Message);
                }
                errorCode = 2;

                return errorCode;
            }

            if (settings.GenerateDefaultTomls)
            {
                if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    Console.WriteLine("Generating default tomls at location: " + settings.OutputFolder);
                }
                CommandLineSettings.GenerateDefaultTaskTomls(settings.OutputFolder);

                return errorCode;
            }

            if (settings.AuditSdrf != null)
            {
                return AuditSdrf(settings, Console.Out);
            }

            if (settings.SdrfDesign != null)
            {
                return WriteDesignFromSdrf(settings, Console.Out);
            }

            // --acceptThermoLicence records the agreement without asking for it, for the situations that
            // cannot answer the prompt further down: a container, a scheduled cluster job, a CI runner,
            // anything with stdin closed. Those are also the situations most likely to need it, because
            // outside a Windows installation SetUpDataDirectory() leaves DataDir at the application
            // folder rather than a per-user one, so the agreement is recorded per extracted copy and
            // every new download, conda environment and container layer starts unagreed.
            if (AgreeToThermoLicence(settings, GlobalVariables.DataDir, Console.Out))
            {
                return errorCode;
            }

            // set up microvignette
            if (settings.RunMicroVignette)
            {
                // set up the spectra file
                settings.Spectra.Clear();
                settings.Spectra.Add(Path.Combine(GlobalVariables.DataDir, @"Data", "SmallCalibratible_Yeast.mzML"));

                // set up the database
                settings.Databases.Clear();
                settings.Databases.Add(Path.Combine(GlobalVariables.DataDir, @"Data", "SmallYeast.fasta"));

                // set up the tasks (calibration, GPTMD, search)
                settings.Tasks.Clear();
                CommandLineSettings.GenerateDefaultTaskTomls(settings.OutputFolder);
                settings.Tasks.Add(Path.Combine(settings.OutputFolder, "CalibrationTask.toml"));
                settings.Tasks.Add(Path.Combine(settings.OutputFolder, "GptmdTask.toml"));
                settings.Tasks.Add(Path.Combine(settings.OutputFolder, "SearchTask.toml"));
            }

            MetaMorpheusEngine.WarnHandler += WarnHandler;
            MetaMorpheusEngine.OutProgressHandler += MyEngine_outProgressHandler;
            MetaMorpheusEngine.StartingSingleEngineHander += MyEngine_startingSingleEngineHander;
            MetaMorpheusEngine.FinishedSingleEngineHandler += MyEngine_finishedSingleEngineHandler;

            MetaMorpheusTask.WarnHandler += WarnHandler;
            MetaMorpheusTask.LogHandler += LogHandler;
            MetaMorpheusTask.StartingSingleTaskHander += MyTaskEngine_startingSingleTaskHander;
            MetaMorpheusTask.FinishedSingleTaskHandler += MyTaskEngine_finishedSingleTaskHandler;
            MetaMorpheusTask.FinishedWritingFileHandler += MyTaskEngine_finishedWritingFileHandler;

            bool containsRawFiles = settings.Spectra.Select(v => Path.GetExtension(v).ToLowerInvariant()).Any(v => v == ".raw");
            if (containsRawFiles && !GlobalVariables.GlobalSettings.UserHasAgreedToThermoRawFileReaderLicence)
            {
                // write the Thermo RawFileReader licence agreement
                Console.WriteLine(ThermoRawFileReaderLicence.ThermoLicenceText);
                Console.WriteLine("\nIn order to search Thermo .raw files, you must agree to the above terms. Do you agree to the above terms? y/n\n");
                string res = Console.ReadLine().ToLowerInvariant();
                if (res == "y")
                {
                    var newGlobalSettings = new GlobalSettings
                    {
                        UserHasAgreedToThermoRawFileReaderLicence = true,
                        WriteExcelCompatibleTSVs = GlobalVariables.GlobalSettings.WriteExcelCompatibleTSVs
                    };

                    Toml.WriteFile<GlobalSettings>(newGlobalSettings, Path.Combine(GlobalVariables.DataDir, @"settings.toml"));
                    GlobalVariables.GlobalSettings = newGlobalSettings;
                }
                else
                {
                    Console.WriteLine("Thermo licence has been declined. Exiting MetaMorpheus. You can still search .mzML and .mgf files without agreeing to the Thermo licence.");
                    errorCode = 3;
                    return errorCode;
                }
            }

            foreach (var db in settings.Databases)
            {
                if (Path.GetExtension(db).Equals(".xml"))
                {
                    GlobalVariables.AddMods(UsefulProteomicsDatabases.ProteinDbLoader.GetPtmListFromProteinXml(db).OfType<Modification>(), true);

                    // print any error messages reading the mods to the console
                    foreach (var error in GlobalVariables.ErrorsReadingMods)
                    {
                        if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                        {
                            Console.WriteLine(error);
                        }
                    }

                    GlobalVariables.ErrorsReadingMods.Clear();
                }
            }

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>();

            var tasks = settings.Tasks.ToList();
            for (int i = 0; i < tasks.Count; i++)
            {
                var filePath = tasks[i];

                var toml = Toml.ReadFile(filePath, MetaMorpheusTask.tomlConfig);

                switch (toml.Get<string>("TaskType"))
                {
                    case "Search":
                        var searchTask = Toml.ReadFile<SearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "SearchTask", searchTask));
                        break;

                    case "Calibrate":
                        var calibrationTask = Toml.ReadFile<CalibrationTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "CalibrationTask", calibrationTask));
                        break;

                    case "Gptmd":
                        var GptmdTask = Toml.ReadFile<GptmdTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "GptmdTask", GptmdTask));
                        break;

                    case "XLSearch":
                        var XlTask = Toml.ReadFile<XLSearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "XLSearchTask", XlTask));
                        break;

                    case "GlycoSearch":
                        var GlycoTask = MetaMorpheusTask.ReadTaskTomlWithLowResFallback<GlycoSearchTask>(filePath);
                        taskList.Add(("Task" + (i + 1) + "GlycoSearchTask", GlycoTask));
                        break;

                    case "Average":
                        var AveragingTask = Toml.ReadFile<SpectralAveragingTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "AveragingTask", AveragingTask));
                        break;

                    case "Truncation":
                        var truncationTask = Toml.ReadFile<TruncationSearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "TruncationSearchTask", truncationTask));
                        break;

                    default:
                        if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                        {
                            Console.WriteLine(toml.Get<string>("TaskType") + " is not a known task type! Skipping.");
                        }
                        break;
                }
            }

            List<string> startingRawFilenameList = settings.Spectra.Select(b => Path.GetFullPath(b)).ToList();
            List<DbForTask> startingXmlDbFilenameList = settings.Databases.Select(b => new DbForTask(Path.GetFullPath(b), IsContaminant(b))).ToList();

            // check that an experimental design is defined if normalization is enabled
            var searchTasks = taskList
                .Where(p => p.Item2.TaskType == MyTask.Search)
                .Select(p => (SearchTask)p.Item2);

            bool reportToConsole = settings.Verbosity == CommandLineSettings.VerbosityType.minimal
                || settings.Verbosity == CommandLineSettings.VerbosityType.normal;

            // An SDRF beside the spectra files makes the design, before the design check reads it.
            int sdrfExitCode = ResolveSdrfDesign(
                Directory.GetParent(startingRawFilenameList.First()).FullName,
                startingRawFilenameList,
                taskList.Where(p => p.Item2.TaskType == MyTask.Search).Select(p => (p.Item1, (SearchTask)p.Item2)).ToList(),
                settings,
                searchTasks.Any(p => p.SearchParameters.Normalize),
                reportToConsole);

            if (sdrfExitCode != 0)
            {
                return sdrfExitCode;
            }

            int designExitCode = ResolveExperimentalDesign(
                Directory.GetParent(startingRawFilenameList.First()).FullName,
                startingRawFilenameList,
                searchTasks.Any(p => p.SearchParameters.Normalize),
                settings.Verbosity == CommandLineSettings.VerbosityType.minimal
                    || settings.Verbosity == CommandLineSettings.VerbosityType.normal);

            if (designExitCode != 0)
            {
                return designExitCode;
            }

            // a task whose digestion gives seed peptides it cannot use would finish with far fewer, wrong identifications
            int digestionExitCode = RefuseSeedDigestion(taskList, Console.Out,
                settings.Verbosity == CommandLineSettings.VerbosityType.minimal
                    || settings.Verbosity == CommandLineSettings.VerbosityType.normal);
            if (digestionExitCode != 0)
            {
                return digestionExitCode;
            }

            EverythingRunnerEngine a = new EverythingRunnerEngine(taskList, startingRawFilenameList, startingXmlDbFilenameList, settings.OutputFolder);

            try
            {
                a.Run();
            }
            catch (Exception e)
            {
                while (e.InnerException != null)
                {
                    e = e.InnerException;
                }

                var message = "Run failed, Exception: " + e.Message;

                if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    Console.WriteLine(message);
                }
                errorCode = 4;
            }

            return errorCode;
        }

        /// <summary>
        /// Refuses the run when a task's digestion settings ask for seed peptides the task cannot use (see
        /// <see cref="MetaMorpheusTask.GetSeedDigestionRefusal"/>). Checked for every task before any runs, as the GUI's Run
        /// button does, and returns the process exit code: 0 to carry on, 6 to stop.
        /// </summary>
        /// <param name="print">Whether to write the reason to <paramref name="output"/>; the run is refused either way.</param>
        public static int RefuseSeedDigestion(List<(string, MetaMorpheusTask)> taskList, TextWriter output, bool print)
        {
            foreach (var (taskName, task) in taskList)
            {
                string refusal = task.GetSeedDigestionRefusal(taskName);
                if (refusal != null)
                {
                    if (print)
                    {
                        output.WriteLine(refusal);
                    }
                    return 6;
                }
            }
            return 0;
        }

        /// <summary>
        /// Works out which experimental design this run will use, reports whatever is wrong with it,
        /// and returns the process exit code: 0 to carry on, 5 to stop.
        ///
        /// Static, and parameterised over its console I/O, so the decision can be tested without
        /// driving a whole command line run. These branches decide whether a search starts at all,
        /// and before this they were reachable only by launching the CLI against real spectra.
        /// </summary>
        /// <param name="designDirectory">The folder beside the spectra files, where a design lives.</param>
        /// <param name="startingRawFilenameList">Full paths of the spectra files this run will search.</param>
        /// <param name="normalizationRequested">True when any search task asks for normalization, which is what makes a design mandatory rather than optional.</param>
        /// <param name="reportToConsole">True at minimal or normal verbosity. At "none" there is nobody to ask, so a recoverable problem is recovered from silently.</param>
        /// <param name="write">Console writer; defaults to <see cref="Console.WriteLine(string)"/>.</param>
        /// <param name="readLine">Console reader; defaults to <see cref="Console.ReadLine"/>.</param>
        public static int ResolveExperimentalDesign(
            string designDirectory,
            List<string> startingRawFilenameList,
            bool normalizationRequested,
            bool reportToConsole,
            Action<string> write = null,
            Func<string> readLine = null)
        {
            write ??= Console.WriteLine;
            readLine ??= Console.ReadLine;

            string pathToExperDesign = Path.Combine(designDirectory, GlobalVariables.ExperimentalDesignFileName);
            string pathToTmtDesign = Path.Combine(designDirectory, GlobalVariables.TmtExperimentalDesignFileName);

            bool hasClassicDesign = File.Exists(pathToExperDesign);
            bool hasTmtDesign = File.Exists(pathToTmtDesign);

            // Only one may be present. Two files carrying replicate structure can drift apart, and
            // nothing here could say which of them the user meant.
            if (hasClassicDesign && hasTmtDesign)
            {
                if (reportToConsole)
                {
                    write("Both ExperimentalDesign and TmtDesign.txt were found. Only one design file type may be present.");
                }

                return 5;
            }

            if (!hasClassicDesign && !hasTmtDesign)
            {
                // A design is optional until something needs it. Normalization is that something.
                if (normalizationRequested)
                {
                    if (reportToConsole)
                    {
                        write("No experimental design file present. Normalization requires a design (ExperimentalDesign.tsv or TmtDesign.txt).");
                    }

                    return 5;
                }

                return 0;
            }

            string designPath = hasClassicDesign ? pathToExperDesign : pathToTmtDesign;
            List<string> errors;
            var warnings = new List<string>();

            if (hasClassicDesign)
            {
                var design = ExperimentalDesign.ReadExperimentalDesign(designPath, startingRawFilenameList, out errors);
                if (!errors.Any())
                {
                    warnings = ExperimentalDesign.GetWarningsInExperimentalDesign(design);
                }
            }
            else
            {
                var tmtFiles = TmtExperimentalDesign.Read(designPath, startingRawFilenameList, out errors);
                if (!errors.Any())
                {
                    warnings = TmtExperimentalDesign.GetWarningsInDesign(tmtFiles);
                }
            }

            if (!errors.Any())
            {
                if (reportToConsole)
                {
                    write("Read " + Path.GetFileName(designPath) + " successfully");
                    foreach (string warning in warnings)
                    {
                        write(warning);
                    }
                }

                return 0;
            }

            // Normalizing against a design that could not be read would silently produce different
            // numbers, so this is the one case that cannot be recovered from.
            if (normalizationRequested)
            {
                if (reportToConsole)
                {
                    foreach (string error in errors)
                    {
                        write(error);
                    }
                }

                return 5;
            }

            // Otherwise the run continues without reading the design, leaving the file alone. It used
            // to be deleted -- silently at verbosity none, and on a y/n whose question was about
            // continuing rather than about deletion. Deleting an input the user cannot be asked about
            // is a lot to infer from a verbosity flag, and #2256 is a user who lost their design that
            // way. The GUI side stops deleting in #2780; the CLI and the GUI have to agree about
            // whether a design that fails to parse is the user's file or ours to discard.
            if (!reportToConsole)
            {
                return 0;
            }

            write((hasClassicDesign ? "An experimental design file" : "A TMT design file")
                + " was found, but an error occurred reading it. Continue without an experimental design? y/n");
            write("First error: " + errors.First());

            string answer = readLine();
            if (answer?.ToLowerInvariant() == "y" || answer?.ToLowerInvariant() == "yes")
            {
                return 0;
            }

            return 5;
        }

        /// <summary>The ending every SDRF file name carries (the SDRF-Proteomics specification's).</summary>
        public const string SdrfFileEnding = ".sdrf.tsv";

        /// <summary>
        /// The automatic design (QuantProject gap AUTODESIGN, ruled 2026-11-07 and 2026-11-11): when an SDRF is beside the
        /// spectra files, the run's experimental design is made from it, and the search's own SDRF is filled from it. With
        /// no SDRF there nothing changes. Runs before <see cref="ResolveExperimentalDesign"/>, which then reads the design
        /// written here exactly as it reads one the user wrote.
        ///
        /// <list type="bullet">
        /// <item>The SDRF is handed to every search task as its <see cref="SearchParameters.InputSdrfPath"/>, with
        /// <see cref="SearchParameters.WriteSdrf"/> on, so the search's own SDRF carries every sample fact the input has. A
        /// task that names a different SDRF of its own keeps it, and a line says so.</item>
        /// <item>The design is read with mzLib's SdrfSearchDesign, for exactly the spectra files given, as --sdrfDesign reads
        /// it. Its condition comes from --sdrfCondition, --sdrfSingleCondition and --sdrfUnknownFactorIsCondition when any
        /// is given, else from the search tasks' settings of the same names, which must agree. The kit is the search tasks'
        /// multiplex label; a plex is read only from --sdrfPlexColumn, --sdrfPlexPattern or --sdrfSinglePlex.</item>
        /// <item>The design is written beside the first spectra file with a <see cref="SdrfDesignMarker"/> holding its and
        /// the SDRF's sha256 (QP-S41). A design already there that the marker shows MetaMorpheus wrote, unedited, is
        /// written again for this run's files; any other design there is the user's, is used as is and never touched, and
        /// a line says the SDRF was not read for it.</item>
        /// <item>A design that cannot be made (refused, unreadable, two SDRFs, search tasks that disagree) is handled as
        /// <see cref="ResolveExperimentalDesign"/> handles a design file it cannot read: with normalization asked for, exit
        /// 5 and every reason; otherwise the reasons and a y/n prompt to continue without a design; at verbosity none, the
        /// run continues. Nothing is written then, and a design MetaMorpheus wrote before, which no longer fits, is
        /// removed with its marker. The SDRF is never edited.</item>
        /// </list>
        /// </summary>
        /// <param name="designDirectory">The folder of the first spectra file, where a run reads its design.</param>
        /// <param name="startingRawFilenameList">Full paths of the spectra files this run will search.</param>
        /// <param name="searchTasks">The run's search tasks, by name; the design serves them all.</param>
        /// <param name="settings">The command line, for the --sdrf... switches.</param>
        /// <param name="normalizationRequested">True when any search task normalizes, which makes a design mandatory.</param>
        /// <param name="reportToConsole">True at minimal or normal verbosity.</param>
        /// <param name="write">Console writer; defaults to <see cref="Console.WriteLine(string)"/>.</param>
        /// <param name="readLine">Console reader; defaults to <see cref="Console.ReadLine"/>.</param>
        /// <returns>0 to carry on, 5 to stop.</returns>
        public static int ResolveSdrfDesign(
            string designDirectory,
            List<string> startingRawFilenameList,
            IReadOnlyList<(string Name, SearchTask Task)> searchTasks,
            CommandLineSettings settings,
            bool normalizationRequested,
            bool reportToConsole,
            Action<string> write = null,
            Func<string> readLine = null)
        {
            write ??= Console.WriteLine;
            readLine ??= Console.ReadLine;

            if (searchTasks.Count == 0)
            {
                return 0;
            }

            var sdrfs = Directory.GetFiles(designDirectory)
                .Where(f => f.EndsWith(SdrfFileEnding, StringComparison.OrdinalIgnoreCase))
                .OrderBy(f => f, StringComparer.OrdinalIgnoreCase)
                .ToList();
            if (sdrfs.Count == 0)
            {
                return 0;
            }

            // Which SDRF, when a folder holds several, is not ruled (QuantProject, 2026-11-07): one is kept per group of
            // files. So none is guessed.
            if (sdrfs.Count > 1)
            {
                return ContinueWithoutSdrfDesign(
                    $"{sdrfs.Count} SDRF files are beside the spectra files ({string.Join(", ", sdrfs.Select(Path.GetFileName))}), "
                    + "and MetaMorpheus does not choose between them, so no experimental design was made from an SDRF. Keep "
                    + "beside the spectra files only the SDRF that describes them.",
                    report: null, normalizationRequested, reportToConsole, write, readLine);
            }

            string sdrf = Path.GetFullPath(sdrfs[0]);
            HandSdrfToSearchTasks(sdrf, searchTasks, reportToConsole, write);

            // A search that quantifies nothing reads no design.
            if (!searchTasks.Any(t => t.Task.SearchParameters.DoLabelFreeQuantification || t.Task.SearchParameters.DoMultiplexQuantification))
            {
                return 0;
            }

            string labelFreePath = Path.Combine(designDirectory, GlobalVariables.ExperimentalDesignFileName);
            string isobaricPath = Path.Combine(designDirectory, GlobalVariables.TmtExperimentalDesignFileName);
            string existing = new[] { labelFreePath, isobaricPath }.Where(File.Exists).ToList() switch
            {
                { Count: 0 } => null,
                { Count: 1 } one => one[0],
                // both kinds: ResolveExperimentalDesign stops the run, and neither is ours to choose between
                _ => string.Empty,
            };
            if (existing == string.Empty)
            {
                return 0;
            }

            var marker = SdrfDesignMarker.ReadFrom(designDirectory);
            bool ours = existing != null && marker != null
                && string.Equals(marker.DesignFileName, Path.GetFileName(existing), StringComparison.Ordinal)
                && marker.DescribesUneditedDesignIn(designDirectory);
            if (existing != null && !ours)
            {
                if (reportToConsole)
                {
                    write($"{Path.GetFileName(existing)} is beside the spectra files and was not written by MetaMorpheus from an SDRF, "
                        + $"or was edited since, so it is the design; the SDRF beside it ({sdrf}) is not read for the design.");
                }
                return 0;
            }

            if (!TryReadGrouping(settings, searchTasks, out var grouping, out string disagreement))
            {
                return ContinueWithoutSdrfDesign(disagreement, report: null, normalizationRequested, reportToConsole, write, readLine,
                    removeOurs: ours ? existing : null, designDirectory);
            }

            var kits = searchTasks
                .Select(t => t.Task.SearchParameters.DoMultiplexQuantification ? t.Task.SearchParameters.MultiplexModId ?? string.Empty : null)
                .Distinct()
                .ToList();
            if (kits.Count > 1)
            {
                return ContinueWithoutSdrfDesign(
                    "The search tasks quantify differently (" + string.Join(", ", searchTasks.Select(t => t.Name + ": "
                        + (t.Task.SearchParameters.DoMultiplexQuantification ? t.Task.SearchParameters.MultiplexModId : "label-free")))
                    + "), and one experimental design cannot serve them all, so none was made from the SDRF.",
                    report: null, normalizationRequested, reportToConsole, write, readLine, removeOurs: ours ? existing : null, designDirectory);
            }

            MzLibIsobaricMassTag tag = null;
            if (kits[0] != null && !MzLibIsobaricMassTag.TryGetIsobaricMassTag(kits[0], GlobalVariables.AllModsKnown, out tag))
            {
                return ContinueWithoutSdrfDesign(
                    $"The search's multiplex label '{kits[0]}' is not an isobaric kit mzLib knows, so the SDRF's channels cannot be read.",
                    report: null, normalizationRequested, reportToConsole, write, readLine, removeOurs: ours ? existing : null, designDirectory);
            }

            var options = new SdrfSearchDesignOptions
            {
                Tag = tag,
                PlexColumn = settings.SdrfPlexColumn,
                PlexFileNamePattern = settings.SdrfPlexPattern,
                SinglePlex = settings.SdrfSinglePlex,
                ConditionColumns = grouping.Columns,
                SingleCondition = grouping.SingleCondition,
                UnknownFactorIsOwnCondition = grouping.UnknownIsOwnCondition,
            };
            var design = SdrfSearchDesign.Read(sdrf, startingRawFilenameList, options);
            string report = design.Report() + SwitchesTheReportNames(design, inATaskFileToo: true);

            if (!design.IsValid)
            {
                return ContinueWithoutSdrfDesign(
                    $"An SDRF is beside the spectra files ({sdrf}), but no experimental design could be made from it, for the reasons above.",
                    report, normalizationRequested, reportToConsole, write, readLine, removeOurs: ours ? existing : null, designDirectory);
            }

            string designPath = tag != null ? isobaricPath : labelFreePath;
            if (existing != null && existing != designPath)
            {
                File.Delete(existing); // ours, unedited, and of the other kind: the search's label changed
            }
            design.WriteDesign(designPath);

            var newMarker = new SdrfDesignMarker(Path.GetFileName(designPath), SdrfDesignMarker.Sha256Of(designPath), sdrf,
                SdrfDesignMarker.Sha256Of(sdrf), DescribeCondition(design, grouping));
            newMarker.WriteTo(designDirectory);
            foreach (var (_, task) in searchTasks)
            {
                task.SearchParameters.DesignWrittenFromSdrf = newMarker.Describe();
            }

            if (reportToConsole)
            {
                write(report.TrimEnd());
                if (tag?.TagType == Omics.Modifications.IsobaricMassTagType.TMT16)
                {
                    write(Tmt16Note);
                }
                write($"Wrote {designPath} from the SDRF {sdrf}" + (ours ? ", again for this run's spectra files." : "."));
            }
            return 0;
        }

        /// <summary>The grouping the user states for a design made from an SDRF: the columns, one condition, or neither.</summary>
        private sealed record SdrfGrouping(List<string> Columns, string SingleCondition, bool UnknownIsOwnCondition);

        /// <summary>
        /// The command line's --sdrfCondition, --sdrfSingleCondition and --sdrfUnknownFactorIsCondition when any is given;
        /// otherwise the search tasks' settings of the same names, which must be the same in every task, because one design
        /// serves them all. The grouping used is written into every search task, so the task files the run writes record it.
        /// </summary>
        private static bool TryReadGrouping(CommandLineSettings settings, IReadOnlyList<(string Name, SearchTask Task)> searchTasks,
            out SdrfGrouping grouping, out string disagreement)
        {
            disagreement = null;
            bool fromCommandLine = settings.SdrfConditionColumns?.Any() == true || settings.SdrfSingleCondition != null
                || settings.SdrfUnknownFactorIsCondition;
            if (fromCommandLine)
            {
                grouping = new SdrfGrouping(settings.SdrfConditionColumns?.ToList() ?? new List<string>(), settings.SdrfSingleCondition,
                    settings.SdrfUnknownFactorIsCondition);
            }
            else
            {
                var stated = searchTasks
                    .Select(t => (t.Name, Grouping: new SdrfGrouping(t.Task.SearchParameters.SdrfConditionColumns?.ToList() ?? new List<string>(),
                        t.Task.SearchParameters.SdrfSingleCondition, t.Task.SearchParameters.SdrfUnknownFactorIsCondition)))
                    .ToList();
                string Key(SdrfGrouping g) => string.Join("\u0001", g.Columns) + "\u0002" + g.SingleCondition + "\u0002" + g.UnknownIsOwnCondition;
                if (stated.Select(s => Key(s.Grouping)).Distinct().Count() > 1)
                {
                    grouping = null;
                    disagreement = "The search tasks build the SDRF's condition differently ("
                        + string.Join("; ", stated.Select(s => $"{s.Name}: {DescribeGrouping(s.Grouping)}"))
                        + "), and one experimental design serves them all, so none was made from the SDRF. Make the tasks' "
                        + "SdrfConditionColumns, SdrfSingleCondition and SdrfUnknownFactorIsCondition the same, or give "
                        + "--sdrfCondition, --sdrfSingleCondition or --sdrfUnknownFactorIsCondition.";
                    return false;
                }
                grouping = stated[0].Grouping;
            }

            foreach (var (_, task) in searchTasks)
            {
                task.SearchParameters.SdrfConditionColumns = grouping.Columns.ToList();
                task.SearchParameters.SdrfSingleCondition = grouping.SingleCondition;
                task.SearchParameters.SdrfUnknownFactorIsCondition = grouping.UnknownIsOwnCondition;
            }
            return true;
        }

        private static string DescribeGrouping(SdrfGrouping g) =>
            (g.SingleCondition != null ? $"one condition, '{g.SingleCondition}'"
                : g.Columns.Count > 0 ? "columns " + string.Join(", ", g.Columns) : "every factor value column")
            + (g.UnknownIsOwnCondition ? ", unknown factors their own condition" : string.Empty);

        /// <summary>How the condition of a design was built, for its marker and results.txt: read or stated.</summary>
        private static string DescribeCondition(SdrfSearchDesign design, SdrfGrouping grouping)
        {
            string single = design.LabelFree?.SingleCondition ?? design.Isobaric?.SingleCondition;
            var columns = design.LabelFree?.ConditionColumns ?? design.Isobaric?.ConditionColumns ?? Array.Empty<string>();
            string text = single != null
                ? $"one condition, '{single}', stated by the user"
                : (grouping.Columns.Count > 0 ? "the column(s) the user named, " : "every factor value column, ") + string.Join(" + ", columns);
            return grouping.UnknownIsOwnCondition
                ? text + "; files whose factor nobody stated are a condition of their own, as stated by the user"
                : text;
        }

        /// <summary>
        /// Gives every search task the SDRF beside the spectra files as its input SDRF, and turns on the search's own SDRF:
        /// with an input SDRF there, the output SDRF is written whether or not the task asked for it (QuantProject, ruled
        /// 2026-11-07). A task naming another SDRF of its own keeps it, and a line says so.
        /// </summary>
        private static void HandSdrfToSearchTasks(string sdrf, IReadOnlyList<(string Name, SearchTask Task)> searchTasks,
            bool reportToConsole, Action<string> write)
        {
            foreach (var (name, task) in searchTasks)
            {
                string own = task.SearchParameters.InputSdrfPath;
                if (string.IsNullOrWhiteSpace(own))
                {
                    task.SearchParameters.InputSdrfPath = sdrf;
                }
                else if (!string.Equals(Path.GetFullPath(own.Trim()), sdrf, StringComparison.OrdinalIgnoreCase) && reportToConsole)
                {
                    write($"{name} names its own InputSdrfPath ({Path.GetFullPath(own.Trim())}), not the SDRF beside the spectra "
                        + $"files ({sdrf}); the search's SDRF copies its samples from the task's own.");
                }

                if (!task.SearchParameters.WriteSdrf)
                {
                    task.SearchParameters.WriteSdrf = true;
                    if (reportToConsole)
                    {
                        write($"{name} writes the search's own SDRF, because an SDRF is beside the spectra files.");
                    }
                }
            }
        }

        /// <summary>
        /// A design that could not be made from the SDRF, handled as <see cref="ResolveExperimentalDesign"/> handles a
        /// design file it cannot read (QuantProject's oracle verdict, 2026-11-07): with normalization asked for, every
        /// reason and exit 5; at verbosity none, carry on; otherwise every reason and a y/n prompt to continue without a
        /// design. A design MetaMorpheus wrote from the SDRF before, unedited, is removed with its marker first: it no
        /// longer fits, and the run must not quantify with it.
        /// </summary>
        private static int ContinueWithoutSdrfDesign(string why, string report, bool normalizationRequested, bool reportToConsole,
            Action<string> write, Func<string> readLine, string removeOurs = null, string designDirectory = null)
        {
            if (removeOurs != null)
            {
                File.Delete(removeOurs);
                File.Delete(SdrfDesignMarker.PathIn(designDirectory));
            }

            if (reportToConsole)
            {
                if (report != null)
                {
                    write(report.TrimEnd());
                }
                if (removeOurs != null)
                {
                    write($"{removeOurs} was removed: MetaMorpheus wrote it from the SDRF before, and it no longer fits.");
                }
            }

            if (normalizationRequested)
            {
                if (reportToConsole)
                {
                    write(why + " Normalization needs an experimental design, so the run stops. No design was written.");
                }
                return 5;
            }

            if (!reportToConsole)
            {
                return 0;
            }

            write(why + " Continue without an experimental design? y/n");
            string answer = readLine();
            return answer?.ToLowerInvariant() == "y" || answer?.ToLowerInvariant() == "yes" ? 0 : 5;
        }

        private static void WriteMultiLineIndented(string toWrite)
        {
            string[] tokens = Regex.Split(toWrite, @"\r?\n|\r");
            foreach (var str in tokens)
            {
                MyWriter.WriteLine(str);
            }
        }

        private static bool IsContaminant(string b)
        {
            if (b.ToUpper().Contains("contaminant".ToUpper())
                || b.ToUpper().Contains("CRAP"))
            {
                return true;
            }

            return false;
        }

        private static void MyTaskEngine_startingSingleTaskHander(object sender, SingleTaskEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("Starting task: " + e.DisplayName);
                MyWriter.Indent++;
            }
        }

        private static void MyTaskEngine_finishedWritingFileHandler(object sender, SingleFileEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("Finished writing file: " + e.WrittenFile);
            }
        }

        /// <summary>
        /// Carries out --sdrfDesign: reads the SDRF for exactly the -s files with mzLib's SdrfSearchDesign, and writes
        /// ExperimentalDesign.tsv (label-free) or TmtDesign.txt (with --sdrfTag) beside the first spectra file, which is
        /// where a run given the same -s looks for it (ResolveExperimentalDesign).
        ///
        /// SdrfSearchDesign is the one reading every SDRF route shares, so a pipeline that checks an SDRF before a search
        /// (PXReprise) and this switch cannot disagree about it: it restricts the SDRF to the searched files with
        /// SdrfSearchScope (a row naming X.raw stands for a search of X.mzML, X-calib.mzML or X-averaged.mzML), then reads
        /// the kept rows with the label-free or isobaric reader, whose rules decide the design. Its report says what the
        /// join set aside and what the reader added, noted or refused.
        ///
        /// Refuses rather than writes anything MetaMorpheus would reject, because an invalid design is worse than none: a
        /// run finding one skips quantification with only a warning. The report is printed whatever the verbosity, as for
        /// --auditSdrf: it is what was asked for, and it is the only record of a biological replicate number added where
        /// the SDRF gives none.
        /// </summary>
        /// <returns>0 when the design was written; 4 when the SDRF could not be read; 5 when the design was refused or a
        /// design file is already there.</returns>
        public static int WriteDesignFromSdrf(CommandLineSettings settings, TextWriter output)
        {
            // -o is not honoured: a design anywhere else is one the run never reads. Said rather than
            // dropped silently, because "Wrote <path>" alone would read as -o having been obeyed.
            if (settings.OutputFolder != null)
            {
                output.WriteLine("-o is ignored with --sdrfDesign: the design is written beside the first spectra file, where a run looks for it.");
            }

            string designDirectory = Directory.GetParent(Path.GetFullPath(settings.Spectra.First())).FullName;
            bool isobaric = settings.SdrfTag != null;
            string designPath = Path.Combine(designDirectory,
                isobaric ? GlobalVariables.TmtExperimentalDesignFileName : GlobalVariables.ExperimentalDesignFileName);

            // Never overwrite a design: it may be the user's own, and a run with both kinds present exits 5.
            foreach (string existing in new[]
                     {
                         Path.Combine(designDirectory, GlobalVariables.ExperimentalDesignFileName),
                         Path.Combine(designDirectory, GlobalVariables.TmtExperimentalDesignFileName)
                     })
            {
                if (File.Exists(existing))
                {
                    output.WriteLine("A design file is already there, and it is not overwritten: " + existing);
                    return 5;
                }
            }

            // validated already, so a kit given is one MetaMorpheus knows
            MzLibIsobaricMassTag tag = null;
            if (isobaric)
            {
                MzLibIsobaricMassTag.TryGetIsobaricMassTag(settings.SdrfTag, out tag);
            }

            var design = SdrfSearchDesign.Read(settings.SdrfDesign, settings.Spectra.Select(Path.GetFullPath).ToList(),
                SdrfDesignOptions(settings, tag));

            output.WriteLine(design.Report());

            if (!design.IsValid)
            {
                output.WriteLine("No design file was written.");

                output.Write(SwitchesTheReportNames(design, inATaskFileToo: false));

                // Neither the join nor a reader ran only when the file itself could not be read: the settings were
                // checked before, so the other such refusals (no searched file, a plex source without a kit) cannot
                // arise here. Told apart by what ran, never by the report's words.
                return design.LabelFree is null && design.Isobaric is null && design.Scoping is null ? 4 : 5;
            }

            // MetaMorpheus has no TMT16 multiplex label: TMTpro 16-plex is searched with TMT18, the same
            // reagent. TmtDesign.txt names no kit, and TMT16's channels are TMT18's first sixteen, so the
            // file is one a TMT18 run reads as is. Said, so the search is not set to a label that is not there.
            if (tag?.TagType == Omics.Modifications.IsobaricMassTagType.TMT16)
            {
                output.WriteLine(Tmt16Note);
            }

            design.WriteDesign(designPath);
            output.WriteLine("Wrote " + designPath);
            return 0;
        }

        /// <summary>
        /// MetaMorpheus has no TMT16 multiplex label: TMTpro 16-plex is searched with TMT18, the same reagent. TmtDesign.txt
        /// names no kit, and TMT16's channels are TMT18's first sixteen, so the file is one a TMT18 run reads as is. Said, so
        /// the search is not set to a label that is not there.
        /// </summary>
        private const string Tmt16Note = "TMT16 is searched in MetaMorpheus with its TMT18 multiplex label (TMTpro, the same reagent). "
            + "This design is read as is by such a run; channels 134C and 135N stay empty.";

        /// <summary>
        /// Where a refusal can be answered by stating a grouping, mzLib's report names its option; this says which switch
        /// (and, for a run, which task setting) that is, one line each, or nothing. Only whether a line is printed hangs on
        /// the report's words.
        /// </summary>
        private static string SwitchesTheReportNames(SdrfSearchDesign design, bool inATaskFileToo)
        {
            var lines = new System.Text.StringBuilder();
            foreach (var (option, commandLine, taskSetting) in new[]
                     {
                         (nameof(SdrfSearchDesignOptions.SingleCondition), "--sdrfSingleCondition", nameof(SearchParameters.SdrfSingleCondition)),
                         (nameof(SdrfSearchDesignOptions.UnknownFactorIsOwnCondition), "--sdrfUnknownFactorIsCondition",
                             nameof(SearchParameters.SdrfUnknownFactorIsCondition)),
                     })
            {
                if (design.Refusals.Any(r => r.Contains("(" + option + ")")))
                {
                    lines.AppendLine($"On the command line, {option} is {commandLine}"
                        + (inATaskFileToo ? $"; in a search task's settings, {taskSetting}." : "."));
                }
            }
            return lines.ToString();
        }

        /// <summary>
        /// What the command line tells SdrfSearchDesign that the SDRF alone cannot say: the kit and where each file's plex
        /// is written (isobaric only), and how the condition is built. The condition columns are passed as given, and
        /// mzLib names them: 'genotype' is 'factor value[genotype]', the 'factor value[' prefix in any case is the
        /// specification's, and any other name with a '[' (such as 'characteristics[cell type]') is that column. With
        /// none, every factor value column makes the condition. A single condition and files whose factor nobody stated
        /// are groupings only the user can state (ruled 2026-11-11): never defaults.
        /// </summary>
        private static SdrfSearchDesignOptions SdrfDesignOptions(CommandLineSettings settings, MzLibIsobaricMassTag tag) => new()
        {
            Tag = tag,
            PlexColumn = settings.SdrfPlexColumn,
            PlexFileNamePattern = settings.SdrfPlexPattern,
            SinglePlex = settings.SdrfSinglePlex,
            ConditionColumns = settings.SdrfConditionColumns?.ToList(),
            SingleCondition = settings.SdrfSingleCondition,
            UnknownFactorIsOwnCondition = settings.SdrfUnknownFactorIsCondition,
        };

        /// <summary>
        /// Prints what one SDRF document says about quantification, and runs nothing else.
        ///
        /// The failure this exists to catch is silent. A file whose comment[label] names the reagent kit
        /// ("TMT10plex" on every row) rather than the channel parses cleanly and validates cleanly, and
        /// then yields one channel for a ten-channel experiment. So does a file with two rows per
        /// (data file, channel), or one that declares no isobaric modification at all. None of that is a
        /// structural defect, so the SDRF validator is right not to report it -- it is a defect in what
        /// the annotation MEANS, and the only way a user finds out today is a wrong number downstream.
        ///
        /// All the reading is mzLib's SdrfQuantAuditor; this prints its report. Scoped to every SDRF
        /// rather than to TMT, because "this is not an isobaric experiment" is a useful answer to get.
        /// </summary>
        /// <returns>0 when the report was written, 4 when the document could not be read.</returns>
        public static int AuditSdrf(CommandLineSettings settings, TextWriter output)
        {
            SdrfQuantAudit audit;

            try
            {
                audit = SdrfQuantAuditor.Audit(settings.AuditSdrf, settings.AuditDataDirectory);
            }
            catch (Exception e)
            {
                // The auditor does not throw on a document it can read; getting here means it could not
                // be read at all -- not an SDRF, truncated, locked. Reported with the path, because the
                // usual cause is pointing at the wrong file.
                if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    output.WriteLine("The SDRF file could not be read: " + settings.AuditSdrf
                                     + Environment.NewLine + e.Message);
                }

                return 4;
            }

            // Written whatever the verbosity, including 'none'. -v suppresses MetaMorpheus talking about
            // itself; the report is the thing that was asked for, and a run that produced it silently
            // would have produced nothing at all.
            output.WriteLine(audit.ToReport());

            return 0;
        }

        /// <summary>
        /// Carries out --acceptThermoLicence: prints the licence, records the agreement, and reports
        /// whether this was the flag on its own, in which case the caller has no run to continue into.
        /// Does nothing and reports false when the flag was not given.
        /// The data directory and the output writer are parameters rather than reached for through
        /// Console and GlobalVariables so that this can be exercised directly, Run() being private.
        /// </summary>
        public static bool AgreeToThermoLicence(CommandLineSettings settings, string dataDirectory, TextWriter output)
        {
            if (!settings.AcceptThermoLicence)
            {
                return false;
            }

            if (!GlobalVariables.GlobalSettings.UserHasAgreedToThermoRawFileReaderLicence)
            {
                // Print the licence whatever the verbosity. The flag is an affirmative act by
                // whoever typed it, and this is what puts the terms they agreed to in the log.
                output.WriteLine(ThermoRawFileReaderLicence.ThermoLicenceText);
                output.WriteLine("\nThe --acceptThermoLicence flag was given, which agrees to the above terms.");

                RecordThermoLicenceAgreement(dataDirectory);
            }

            if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                output.WriteLine("Agreed to the Thermo RawFileReader licence. Recorded in "
                    + Path.Combine(dataDirectory, @"settings.toml"));
            }

            // On its own the flag is a setup step, so there is no run to continue into.
            return settings.Tasks.Count + settings.Databases.Count + settings.Spectra.Count == 0;
        }

        /// <summary>
        /// Records agreement to the Thermo RawFileReader licence, in memory and in settings.toml,
        /// carrying the one other setting that file holds so recording the agreement does not reset it.
        /// </summary>
        private static void RecordThermoLicenceAgreement(string dataDirectory)
        {
            var newGlobalSettings = new GlobalSettings
            {
                UserHasAgreedToThermoRawFileReaderLicence = true,
                WriteExcelCompatibleTSVs = GlobalVariables.GlobalSettings.WriteExcelCompatibleTSVs
            };

            Toml.WriteFile<GlobalSettings>(newGlobalSettings, Path.Combine(dataDirectory, @"settings.toml"));
            GlobalVariables.GlobalSettings = newGlobalSettings;
        }

        private static void MyTaskEngine_finishedSingleTaskHandler(object sender, SingleTaskEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal || CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.minimal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;
            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal || CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.minimal)
            {
                MyWriter.Indent--;
                WriteMultiLineIndented("Finished task: " + e.DisplayName);
            }
        }

        private static void MyEngine_startingSingleEngineHander(object sender, SingleEngineEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("Starting engine: " + e.MyEngine.GetType().Name + " " + e.MyEngine.GetId());
                MyWriter.Indent++;
            }
        }

        private static void MyEngine_finishedSingleEngineHandler(object sender, SingleEngineFinishedEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("Engine results: " + e);
                MyWriter.Indent--;
                WriteMultiLineIndented("Finished engine: " + e.MyResults.MyEngine.GetType().Name + " " + e.MyResults.MyEngine.GetId());
            }
        }

        private static void MyEngine_outProgressHandler(object sender, ProgressEventArgs e)
        {
            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                MyWriter.Write(e.NewProgress + " ");
            }
            InProgress = true;
        }

        private static void WarnHandler(object sender, StringEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.minimal || CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.minimal || CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("WARN: " + e.S);
            }
        }

        private static void LogHandler(object sender, StringEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("Log: " + e.S);
            }
        }
    }
}
