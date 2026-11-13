using EngineLayer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using SdrfSearchDesign = Readers.SdrfSearchDesign;
using SdrfSearchDesignOptions = Readers.SdrfSearchDesignOptions;
using MzLibIsobaricMassTag = Omics.Modifications.IsobaricMassTag;

namespace TaskLayer
{
    /// <summary>
    /// What the user stated about an SDRF's design, outside the SDRF: the condition (columns, one condition, or files
    /// whose factor nobody stated as their own condition) and where each file's plex is written. The command line fills
    /// it from its --sdrf... switches; empty means nothing was stated there, and the search tasks' settings are read.
    /// </summary>
    public sealed record SdrfDesignStatements
    {
        public IReadOnlyList<string> ConditionColumns { get; init; } = Array.Empty<string>();
        public string SingleCondition { get; init; }
        public bool UnknownFactorIsOwnCondition { get; init; }
        public string PlexColumn { get; init; }
        public string PlexPattern { get; init; }
        public string SinglePlex { get; init; }

        internal bool StatesGrouping => ConditionColumns.Any() || SingleCondition != null || UnknownFactorIsOwnCondition;
    }

    /// <summary>
    /// How a front end names the ways to state a grouping, for the lines <see cref="SdrfAutomaticDesign.Resolve"/> writes:
    /// the command line names its switches, the GUI its own controls.
    /// </summary>
    public sealed record SdrfDesignWording
    {
        /// <summary>Lines added after mzLib's report, saying how to answer a refusal it names (may return empty).</summary>
        public Func<SdrfSearchDesign, string> ReportHint { get; init; } = _ => string.Empty;

        /// <summary>The end of the sentence refusing search tasks that build the condition differently.</summary>
        public string HowToAgreeOnGrouping { get; init; } =
            "Make the tasks' SdrfConditionColumns, SdrfSingleCondition and SdrfUnknownFactorIsCondition the same.";
    }

    /// <summary>
    /// The search-task settings <see cref="SdrfAutomaticDesign.Resolve"/> changes, taken before a run so the GUI can put
    /// them back when it ends (QuantProject, ruled 2026-11-13): the SDRF describes the spectra files, not the task, so a
    /// task loaded in the GUI must not keep pointing at one run's SDRF, or keep Write SDRF on, for the next run's files.
    /// The run's own task files record the settings it used.
    /// </summary>
    public sealed class SdrfTaskSettingsSnapshot
    {
        private readonly List<(SearchParameters Parameters, string InputSdrfPath, bool WriteSdrf, List<string> ConditionColumns,
            string SingleCondition, bool UnknownFactorIsCondition, string DesignWrittenFromSdrf)> _saved;

        private SdrfTaskSettingsSnapshot(List<(SearchParameters, string, bool, List<string>, string, bool, string)> saved) =>
            _saved = saved;

        /// <summary>The SDRF settings of each of <paramref name="searchTasks"/> as they are now.</summary>
        public static SdrfTaskSettingsSnapshot Take(IEnumerable<SearchTask> searchTasks) =>
            new(searchTasks.Select(t => t.SearchParameters)
                .Select(p => (p, p.InputSdrfPath, p.WriteSdrf, p.SdrfConditionColumns?.ToList(), p.SdrfSingleCondition,
                    p.SdrfUnknownFactorIsCondition, p.DesignWrittenFromSdrf))
                .ToList());

        /// <summary>Puts every setting back as it was when the snapshot was taken.</summary>
        public void Restore()
        {
            foreach (var s in _saved)
            {
                s.Parameters.InputSdrfPath = s.InputSdrfPath;
                s.Parameters.WriteSdrf = s.WriteSdrf;
                s.Parameters.SdrfConditionColumns = s.ConditionColumns?.ToList();
                s.Parameters.SdrfSingleCondition = s.SingleCondition;
                s.Parameters.SdrfUnknownFactorIsCondition = s.UnknownFactorIsCondition;
                s.Parameters.DesignWrittenFromSdrf = s.DesignWrittenFromSdrf;
            }
        }
    }

    /// <summary>What <see cref="SdrfAutomaticDesign.Resolve"/> did, and what a front end has to tell the user.</summary>
    public sealed class SdrfDesignOutcome
    {
        internal SdrfDesignOutcome(List<string> messages, string refusal, string sdrfPath, string designPath)
        {
            Messages = messages;
            Refusal = refusal;
            SdrfPath = sdrfPath;
            DesignPath = designPath;
        }

        /// <summary>The lines to show, in order (nothing when no SDRF is beside the spectra files).</summary>
        public IReadOnlyList<string> Messages { get; }

        /// <summary>
        /// Why no design could be made from the SDRF, when one should have been; null otherwise. A front end handles it as
        /// a design file it cannot read: stop when normalizing, otherwise ask whether to continue without a design.
        /// </summary>
        public string Refusal { get; }

        /// <summary>The SDRF beside the spectra files, as a full path; null when there is none (or several).</summary>
        public string SdrfPath { get; }

        /// <summary>The design written from the SDRF; null when none was.</summary>
        public string DesignPath { get; }
    }

    /// <summary>
    /// The automatic design (QuantProject gap AUTODESIGN, ruled 2026-11-07, 2026-11-11 and 2026-11-13): when an SDRF is
    /// beside the spectra files (in the first spectra file's folder, named by <see cref="IsSdrfFileName"/>), the run's
    /// experimental design is made from it, and the search's own SDRF is filled from it. With no SDRF there nothing
    /// changes. Shared by the command line and the GUI, which only show the <see cref="SdrfDesignOutcome"/>; it runs before
    /// the run's design check, which then reads the design written here exactly as it reads one the user wrote.
    ///
    /// <list type="bullet">
    /// <item>The SDRF is handed to every search task as its <see cref="SearchParameters.InputSdrfPath"/>, with
    /// <see cref="SearchParameters.WriteSdrf"/> on, so the search's own SDRF carries every sample fact the input has. A
    /// task that names a different SDRF of its own keeps it, and a line says so.</item>
    /// <item>The design is read with mzLib's SdrfSearchDesign, for exactly the spectra files given. Its condition comes
    /// from <see cref="SdrfDesignStatements"/> when they state one, else from the search tasks' settings of the same
    /// names, which must agree. The kit is the search tasks' multiplex label; the plex is read from a stated column,
    /// pattern or single plex, else from the SDRF's <see cref="DefaultPlexColumn"/>.</item>
    /// <item>The design is written beside the first spectra file with a <see cref="SdrfDesignMarker"/> holding its and
    /// the SDRF's sha256 (QP-S41). A design already there that the marker shows MetaMorpheus wrote, unedited, is written
    /// again for this run's files; any other design there is the user's, is used as is and never touched, and a line says
    /// the SDRF was not read for it.</item>
    /// <item>A design that cannot be made (refused, unreadable, two SDRFs, search tasks that disagree) is returned as a
    /// <see cref="SdrfDesignOutcome.Refusal"/>. Nothing is written then, and a design MetaMorpheus wrote before, which
    /// no longer fits, is removed with its marker. The SDRF is never edited.</item>
    /// </list>
    /// </summary>
    public static class SdrfAutomaticDesign
    {
        /// <summary>The ending every SDRF file name carries (the SDRF-Proteomics specification's).</summary>
        public const string SdrfFileEnding = ".sdrf.tsv";

        /// <summary>
        /// A file named just this is an SDRF too (QuantProject, ruled 2026-11-13). It lacks the dot
        /// <see cref="SdrfFileEnding"/> starts with, and it is the name a person most often gives one, so missing it would
        /// make no design and say nothing.
        /// </summary>
        public const string SdrfBareFileName = "sdrf.tsv";

        /// <summary>
        /// TMT16 has no multiplex label of its own in MetaMorpheus: its channels are TMT18's first sixteen, so the file is
        /// one a TMT18 run reads as is. Said, so the search is not set to a label that is not there.
        /// </summary>
        public const string Tmt16Note = "TMT16 is searched in MetaMorpheus with its TMT18 multiplex label (TMTpro, the same reagent). "
            + "This design is read as is by such a run; channels 134C and 135N stay empty.";

        /// <summary>
        /// The column a TMT SDRF's plex is read from when no plex source is stated (QuantProject, ruled 2026-11-13): the
        /// column PXReprise writes in every TMT SDRF it drafts. A fixed column name is the SDRF stating the plex, not
        /// MetaMorpheus guessing it; a stated column, file-name pattern or single plex wins over it.
        /// </summary>
        public const string DefaultPlexColumn = "comment[plex]";

        /// <summary>
        /// The SDRF's own spelling of <see cref="DefaultPlexColumn"/> (case ignored), or null when it has no such column or
        /// cannot be read (the SDRF reader then says why).
        /// </summary>
        public static string DefaultPlexColumnIn(string sdrfPath)
        {
            string header;
            try
            {
                header = File.ReadLines(sdrfPath).FirstOrDefault() ?? string.Empty;
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException or ArgumentException or NotSupportedException)
            {
                return null;
            }
            return header.TrimStart('\uFEFF').Split('\t').Select(c => c.Trim())
                .FirstOrDefault(c => string.Equals(c, DefaultPlexColumn, StringComparison.OrdinalIgnoreCase));
        }

        /// <summary>
        /// Where a TMT design's plex is read: the stated column, pattern or single plex when one is stated, otherwise the
        /// SDRF's <see cref="DefaultPlexColumn"/> when it has one, otherwise nothing (mzLib then refuses, naming what is
        /// missing; <see cref="DefaultPlexNote"/> adds that comment[plex] was looked for).
        /// </summary>
        public static (string Column, string Pattern, string SinglePlex) PlexSource(string sdrfPath, string column, string pattern,
            string singlePlex) =>
            column == null && pattern == null && singlePlex == null
                ? (DefaultPlexColumnIn(sdrfPath), null, null)
                : (column, pattern, singlePlex);

        /// <summary>
        /// A line for a report mzLib refused because no plex source was declared: the SDRF has no comment[plex], the column
        /// read when none is named. Empty otherwise.
        /// </summary>
        public static string DefaultPlexNote(SdrfSearchDesign design) =>
            design.Refusals.Any(r => r.StartsWith("No plex source was declared", StringComparison.Ordinal))
                ? $"The SDRF has no {DefaultPlexColumn} column, which MetaMorpheus reads for the plex when none is named."
                    + Environment.NewLine
                : string.Empty;

        /// <summary>
        /// True for a file the automatic design reads as an SDRF: a name ending <see cref="SdrfFileEnding"/>, or exactly
        /// <see cref="SdrfBareFileName"/>. Case is ignored.
        /// </summary>
        public static bool IsSdrfFileName(string path)
        {
            string name = Path.GetFileName(path);
            return name.EndsWith(SdrfFileEnding, StringComparison.OrdinalIgnoreCase)
                || string.Equals(name, SdrfBareFileName, StringComparison.OrdinalIgnoreCase);
        }

        /// <summary>
        /// True for a file the GUI recognises when it is dropped, though it adds nothing to the run itself: an SDRF, an
        /// experimental design file, or the marker. Each is read from the spectra files' folder, never from where it was
        /// dropped (QuantProject, ruled 2026-11-13).
        /// </summary>
        public static bool IsDesignOrSdrfFile(string path)
        {
            string name = Path.GetFileName(path);
            return IsSdrfFileName(path)
                || string.Equals(name, GlobalVariables.ExperimentalDesignFileName, StringComparison.OrdinalIgnoreCase)
                || string.Equals(name, GlobalVariables.TmtExperimentalDesignFileName, StringComparison.OrdinalIgnoreCase)
                || string.Equals(name, SdrfDesignMarker.FileName, StringComparison.OrdinalIgnoreCase);
        }

        /// <summary>
        /// The line to show for a dropped SDRF or design file that is not in <paramref name="designDirectory"/>, the folder of
        /// the first spectra file, and so is not used; null when it is there.
        /// </summary>
        /// <param name="path">The dropped file.</param>
        /// <param name="designDirectory">The first checked spectra file's folder; null when no spectra file is loaded.</param>
        public static string NoteForDroppedFile(string path, string designDirectory)
        {
            string folder = Path.GetDirectoryName(Path.GetFullPath(path));
            if (designDirectory != null && string.Equals(folder, Path.TrimEndingDirectorySeparator(Path.GetFullPath(designDirectory)),
                    StringComparison.OrdinalIgnoreCase))
            {
                return null;
            }
            return $"{Path.GetFileName(path)} is not used: SDRF and experimental design files are read from the folder of the first "
                + "spectra file" + (designDirectory != null ? $" ({designDirectory}), and it is in {folder}." : ", and no spectra file is loaded.");
        }

        /// <summary>
        /// Makes the run's design from the SDRF beside the spectra files, when there is one (see the class summary), and says
        /// what it did. Changes the search tasks (their input SDRF, Write SDRF and grouping settings) and the design files
        /// beside the spectra; never the SDRF.
        /// </summary>
        /// <param name="designDirectory">The folder of the first spectra file, where a run reads its design.</param>
        /// <param name="spectraFilePaths">Full paths of the spectra files this run will search.</param>
        /// <param name="searchTasks">The run's search tasks, by name; the design serves them all.</param>
        /// <param name="stated">What the user stated outside the SDRF and the tasks (the command line's switches).</param>
        /// <param name="wording">How the front end names the ways to state a grouping.</param>
        public static SdrfDesignOutcome Resolve(string designDirectory, List<string> spectraFilePaths,
            IReadOnlyList<(string Name, SearchTask Task)> searchTasks, SdrfDesignStatements stated, SdrfDesignWording wording)
        {
            stated ??= new SdrfDesignStatements();
            wording ??= new SdrfDesignWording();
            var messages = new List<string>();

            if (searchTasks.Count == 0)
            {
                return new SdrfDesignOutcome(messages, null, null, null);
            }

            if (!TryFindSdrf(designDirectory, out string sdrf, out string several))
            {
                return new SdrfDesignOutcome(messages, several, null, null);
            }

            HandSdrfToSearchTasks(sdrf, searchTasks, messages);

            // A search that quantifies nothing reads no design.
            if (!searchTasks.Any(t => t.Task.SearchParameters.DoLabelFreeQuantification || t.Task.SearchParameters.DoMultiplexQuantification))
            {
                return new SdrfDesignOutcome(messages, null, sdrf, null);
            }

            if (!TryFindOurDesign(designDirectory, sdrf, messages, out string existing, out bool ours))
            {
                return new SdrfDesignOutcome(messages, null, sdrf, null);
            }
            string removeOurs = ours ? existing : null;

            if (!TryReadGrouping(stated, searchTasks, wording, out var grouping, out string disagreement))
            {
                return Refuse(messages, disagreement, null, sdrf, removeOurs, designDirectory);
            }

            var kits = searchTasks
                .Select(t => t.Task.SearchParameters.DoMultiplexQuantification ? t.Task.SearchParameters.MultiplexModId ?? string.Empty : null)
                .Distinct()
                .ToList();
            if (kits.Count > 1)
            {
                return Refuse(messages,
                    "The search tasks quantify differently (" + string.Join(", ", searchTasks.Select(t => t.Name + ": "
                        + (t.Task.SearchParameters.DoMultiplexQuantification ? t.Task.SearchParameters.MultiplexModId : "label-free")))
                    + "), and one experimental design cannot serve them all, so none was made from the SDRF.",
                    null, sdrf, removeOurs, designDirectory);
            }

            // The kit is the search task's multiplex label. A design made with a kit guessed from the SDRF when the files
            // were added is made again with it here; an SDRF naming a channel this kit does not have is refused by mzLib.
            MzLibIsobaricMassTag tag = null;
            if (kits[0] != null && !MzLibIsobaricMassTag.TryGetIsobaricMassTag(kits[0], GlobalVariables.AllModsKnown, out tag))
            {
                return Refuse(messages,
                    $"The search's multiplex label '{kits[0]}' is not an isobaric kit mzLib knows, so the SDRF's channels cannot be read.",
                    null, sdrf, removeOurs, designDirectory);
            }

            var outcome = MakeDesign(designDirectory, spectraFilePaths, sdrf, tag, grouping, stated, wording, existing, ours, messages);
            if (outcome.DesignPath != null)
            {
                string record = SdrfDesignMarker.ReadFrom(designDirectory)?.Describe();
                foreach (var (_, task) in searchTasks)
                {
                    task.SearchParameters.DesignWrittenFromSdrf = record;
                }
            }
            return outcome;
        }

        /// <summary>
        /// The GUI's automatic design when spectra files are added (QuantProject, ruled 2026-11-13), before any search task
        /// need be loaded: the same SDRF, design and marker rules as <see cref="Resolve"/>, with the kit guessed from the
        /// SDRF's channel labels (<see cref="GuessIsobaricKit"/>) and the condition from every factor value column, since no
        /// task states one yet. Nothing is handed to a task; <see cref="Resolve"/> does that at Run, with the tasks' own kit
        /// and grouping, making the design again when they differ. A design that cannot be made is returned as
        /// <see cref="SdrfDesignOutcome.Refusal"/>, for the GUI to note; a kit that cannot be guessed only defers the design to
        /// Run.
        /// </summary>
        /// <param name="designDirectory">The folder of the first checked spectra file, where a run reads its design.</param>
        /// <param name="spectraFilePaths">Full paths of the checked spectra files.</param>
        /// <param name="wording">How the front end names the ways to state a grouping.</param>
        public static SdrfDesignOutcome ResolveOnAddingFiles(string designDirectory, List<string> spectraFilePaths,
            SdrfDesignWording wording)
        {
            wording ??= new SdrfDesignWording();
            var messages = new List<string>();

            if (spectraFilePaths.Count == 0)
            {
                return new SdrfDesignOutcome(messages, null, null, null);
            }
            if (!TryFindSdrf(designDirectory, out string sdrf, out string several))
            {
                return new SdrfDesignOutcome(messages, several, null, null);
            }

            if (!TryFindOurDesign(designDirectory, sdrf, messages, out string existing, out bool ours))
            {
                return new SdrfDesignOutcome(messages, null, sdrf, null);
            }

            if (!GuessIsobaricKit(sdrf, out MzLibIsobaricMassTag tag, out string undecided))
            {
                messages.Add(undecided + " The experimental design is made from the SDRF at Run, with the search task's label.");
                return new SdrfDesignOutcome(messages, null, sdrf, null);
            }

            var everyFactorColumn = new SdrfGrouping(new List<string>(), null, false);
            return MakeDesign(designDirectory, spectraFilePaths, sdrf, tag, everyFactorColumn, new SdrfDesignStatements(), wording,
                existing, ours, messages);
        }

        /// <summary>
        /// The kit an SDRF's channel labels name, guessed (QuantProject, ruled 2026-11-13, for the GUI before a search task
        /// is known): every comment[label] cell 'label free sample' is label-free (true, <paramref name="tag"/> null);
        /// otherwise the kit with the fewest channels that holds every label (TMT10's TMT131 and TMT11's TMT131N/131C are told
        /// apart by mzLib's label table, <c>SdrfIsobaricLabels</c>). False, with the reason, when the SDRF has no label
        /// column, mixes label-free and labelled rows, or names labels no single kit holds.
        /// </summary>
        public static bool GuessIsobaricKit(string sdrfPath, out MzLibIsobaricMassTag tag, out string undecided)
        {
            tag = null;
            undecided = null;
            List<string> labels;
            try
            {
                var lines = File.ReadLines(sdrfPath).Where(l => l.Length > 0).ToList();
                var header = lines.Count == 0 ? Array.Empty<string>() : lines[0].TrimStart('\uFEFF').Split('\t');
                int column = Array.FindIndex(header, h => string.Equals(h.Trim(), LabelColumn, StringComparison.OrdinalIgnoreCase));
                if (column < 0)
                {
                    undecided = $"The SDRF {sdrfPath} has no {LabelColumn} column, so its kit cannot be told.";
                    return false;
                }
                labels = lines.Skip(1).Select(l => l.Split('\t')).Where(cells => cells.Length > column)
                    .Select(cells => cells[column].Trim()).Where(c => c.Length > 0)
                    .Distinct(StringComparer.OrdinalIgnoreCase).ToList();
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                undecided = $"The SDRF {sdrfPath} cannot be read: {e.Message}";
                return false;
            }

            var labelFree = labels.Where(l => l.IndexOf(LabelFreeSample, StringComparison.OrdinalIgnoreCase) >= 0).ToList();
            if (labels.Count > 0 && labelFree.Count == labels.Count)
            {
                return true;
            }
            if (labelFree.Count > 0 || labels.Count == 0)
            {
                undecided = $"The SDRF {sdrfPath} " + (labels.Count == 0 ? "names no labels" : "mixes label-free and labelled rows")
                    + ", so its kit cannot be told.";
                return false;
            }

            var kit = Enum.GetValues<Omics.Modifications.IsobaricMassTagType>()
                .OrderBy(k => MzLibIsobaricMassTag.GetReporterIonLabels(k).Count)
                .Cast<Omics.Modifications.IsobaricMassTagType?>()
                .FirstOrDefault(k => labels.All(l => Readers.SdrfIsobaricLabels.TryGetChannel(k.Value, l, out _)));
            if (kit == null || !MzLibIsobaricMassTag.TryGetIsobaricMassTag(kit.Value, out tag)) // mzLib's kit, as --sdrfDesign takes it
            {
                undecided = $"No isobaric kit holds every label the SDRF {sdrfPath} names ({string.Join(", ", labels)}), so its kit cannot be told.";
                tag = null;
                return false;
            }
            return true;
        }

        /// <summary>The SDRF column of each row's label, and the label of a label-free sample.</summary>
        private const string LabelColumn = "comment[label]";
        private const string LabelFreeSample = "label free sample";

        /// <summary>
        /// The one SDRF beside the spectra files. False with no SDRF (<paramref name="several"/> null) or with several
        /// (<paramref name="several"/> says so): which SDRF, when a folder holds several, is not ruled (QuantProject,
        /// 2026-11-07), so none is guessed.
        /// </summary>
        private static bool TryFindSdrf(string designDirectory, out string sdrf, out string several)
        {
            sdrf = null;
            several = null;
            var sdrfs = Directory.GetFiles(designDirectory)
                .Where(IsSdrfFileName)
                .OrderBy(f => f, StringComparer.OrdinalIgnoreCase)
                .ToList();
            if (sdrfs.Count == 1)
            {
                sdrf = Path.GetFullPath(sdrfs[0]);
                return true;
            }
            if (sdrfs.Count > 1)
            {
                several = $"{sdrfs.Count} SDRF files are beside the spectra files ({string.Join(", ", sdrfs.Select(Path.GetFileName))}), "
                    + "and MetaMorpheus does not choose between them, so no experimental design was made from an SDRF. Keep "
                    + "beside the spectra files only the SDRF that describes them.";
            }
            return false;
        }

        /// <summary>
        /// Whether the design beside the spectra files may be made from the SDRF: true with none there, or with one the
        /// marker shows MetaMorpheus wrote, unedited (<paramref name="ours"/>). False when both kinds are there (the design
        /// check stops the run, and neither is ours to choose between) or the design is the user's (a line says the SDRF
        /// is not read for it).
        /// </summary>
        private static bool TryFindOurDesign(string designDirectory, string sdrf, List<string> messages, out string existing, out bool ours)
        {
            string labelFreePath = Path.Combine(designDirectory, GlobalVariables.ExperimentalDesignFileName);
            string isobaricPath = Path.Combine(designDirectory, GlobalVariables.TmtExperimentalDesignFileName);
            var present = new[] { labelFreePath, isobaricPath }.Where(File.Exists).ToList();
            existing = present.Count == 1 ? present[0] : null;
            ours = false;
            if (present.Count > 1)
            {
                return false;
            }

            var marker = SdrfDesignMarker.ReadFrom(designDirectory);
            ours = existing != null && marker != null
                && string.Equals(marker.DesignFileName, Path.GetFileName(existing), StringComparison.Ordinal)
                && marker.DescribesUneditedDesignIn(designDirectory);
            if (existing != null && !ours)
            {
                messages.Add($"{Path.GetFileName(existing)} is beside the spectra files and was not written by MetaMorpheus from an SDRF, "
                    + $"or was edited since, so it is the design; the SDRF beside it ({sdrf}) is not read for the design.");
                return false;
            }
            return true;
        }

        /// <summary>
        /// Reads the design from the SDRF and writes it with its marker, or refuses. A design of ours whose content would not
        /// change is left as it is (QuantProject, ruled 2026-11-13: re-made only when it no longer fits), and its marker is
        /// brought up to date.
        /// </summary>
        private static SdrfDesignOutcome MakeDesign(string designDirectory, List<string> spectraFilePaths, string sdrf,
            MzLibIsobaricMassTag tag, SdrfGrouping grouping, SdrfDesignStatements stated, SdrfDesignWording wording,
            string existing, bool ours, List<string> messages)
        {
            string removeOurs = ours ? existing : null;
            var plex = tag != null
                ? PlexSource(sdrf, stated.PlexColumn, stated.PlexPattern, stated.SinglePlex)
                : (stated.PlexColumn, stated.PlexPattern, stated.SinglePlex);
            var options = new SdrfSearchDesignOptions
            {
                Tag = tag,
                PlexColumn = plex.Item1,
                PlexFileNamePattern = plex.Item2,
                SinglePlex = plex.Item3,
                ConditionColumns = grouping.Columns,
                SingleCondition = grouping.SingleCondition,
                UnknownFactorIsOwnCondition = grouping.UnknownIsOwnCondition,
            };
            var design = SdrfSearchDesign.Read(sdrf, spectraFilePaths, options);
            string report = design.Report() + DefaultPlexNote(design) + wording.ReportHint(design);

            if (!design.IsValid)
            {
                return Refuse(messages,
                    $"An SDRF is beside the spectra files ({sdrf}), but no experimental design could be made from it, for the reasons above.",
                    report, sdrf, removeOurs, designDirectory);
            }

            string designPath = Path.Combine(designDirectory,
                tag != null ? GlobalVariables.TmtExperimentalDesignFileName : GlobalVariables.ExperimentalDesignFileName);
            string newPath = designPath + ".new";
            bool unchanged;
            try
            {
                design.WriteDesign(newPath);
                unchanged = ours && existing == designPath
                    && SdrfDesignMarker.ContentSha256Of(newPath) == SdrfDesignMarker.ContentSha256Of(designPath);
                if (!unchanged)
                {
                    if (existing != null && existing != designPath)
                    {
                        File.Delete(existing); // ours, unedited, and of the other kind: the search's label changed
                    }
                    File.Move(newPath, designPath, overwrite: true);
                }
            }
            finally
            {
                if (File.Exists(newPath))
                {
                    File.Delete(newPath);
                }
            }

            var newMarker = new SdrfDesignMarker(Path.GetFileName(designPath), SdrfDesignMarker.ContentSha256Of(designPath), sdrf,
                SdrfDesignMarker.Sha256Of(sdrf), DescribeCondition(design, grouping));
            newMarker.WriteTo(designDirectory);

            messages.Add(report.TrimEnd());
            if (tag?.TagType == Omics.Modifications.IsobaricMassTagType.TMT16)
            {
                messages.Add(Tmt16Note);
            }
            messages.Add(unchanged
                ? $"{designPath}, made by MetaMorpheus from the SDRF {sdrf}, already fits this run's spectra files."
                : $"Wrote {designPath} from the SDRF {sdrf}" + (ours ? ", again for this run's spectra files." : "."));
            return new SdrfDesignOutcome(messages, null, sdrf, designPath);
        }

        /// <summary>
        /// A design that could not be made from the SDRF. A design MetaMorpheus wrote from the SDRF before, unedited, is
        /// removed with its marker first: it no longer fits, and the run must not quantify with it.
        /// </summary>
        private static SdrfDesignOutcome Refuse(List<string> messages, string why, string report, string sdrf, string removeOurs,
            string designDirectory)
        {
            if (removeOurs != null)
            {
                File.Delete(removeOurs);
                File.Delete(SdrfDesignMarker.PathIn(designDirectory));
            }

            if (report != null)
            {
                messages.Add(report.TrimEnd());
            }
            if (removeOurs != null)
            {
                messages.Add($"{removeOurs} was removed: MetaMorpheus wrote it from the SDRF before, and it no longer fits.");
            }
            return new SdrfDesignOutcome(messages, why, sdrf, null);
        }

        /// <summary>The grouping the user states for a design made from an SDRF: the columns, one condition, or neither.</summary>
        private sealed record SdrfGrouping(List<string> Columns, string SingleCondition, bool UnknownIsOwnCondition);

        /// <summary>
        /// The stated grouping when there is one; otherwise the search tasks' settings of the same names, which must be the
        /// same in every task, because one design serves them all. The grouping used is written into every search task, so
        /// the task files the run writes record it.
        /// </summary>
        private static bool TryReadGrouping(SdrfDesignStatements stated, IReadOnlyList<(string Name, SearchTask Task)> searchTasks,
            SdrfDesignWording wording, out SdrfGrouping grouping, out string disagreement)
        {
            disagreement = null;
            if (stated.StatesGrouping)
            {
                grouping = new SdrfGrouping(stated.ConditionColumns.ToList(), stated.SingleCondition, stated.UnknownFactorIsOwnCondition);
            }
            else
            {
                var fromTasks = searchTasks
                    .Select(t => (t.Name, Grouping: new SdrfGrouping(t.Task.SearchParameters.SdrfConditionColumns?.ToList() ?? new List<string>(),
                        t.Task.SearchParameters.SdrfSingleCondition, t.Task.SearchParameters.SdrfUnknownFactorIsCondition)))
                    .ToList();
                string Key(SdrfGrouping g) => string.Join("\u0001", g.Columns) + "\u0002" + g.SingleCondition + "\u0002" + g.UnknownIsOwnCondition;
                if (fromTasks.Select(s => Key(s.Grouping)).Distinct().Count() > 1)
                {
                    grouping = null;
                    disagreement = "The search tasks build the SDRF's condition differently ("
                        + string.Join("; ", fromTasks.Select(s => $"{s.Name}: {DescribeGrouping(s.Grouping)}"))
                        + "), and one experimental design serves them all, so none was made from the SDRF. "
                        + wording.HowToAgreeOnGrouping;
                    return false;
                }
                grouping = fromTasks[0].Grouping;
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
            List<string> messages)
        {
            foreach (var (name, task) in searchTasks)
            {
                string own = task.SearchParameters.InputSdrfPath;
                if (string.IsNullOrWhiteSpace(own))
                {
                    task.SearchParameters.InputSdrfPath = sdrf;
                }
                else if (!string.Equals(Path.GetFullPath(own.Trim()), sdrf, StringComparison.OrdinalIgnoreCase))
                {
                    messages.Add($"{name} names its own InputSdrfPath ({Path.GetFullPath(own.Trim())}), not the SDRF beside the spectra "
                        + $"files ({sdrf}); the search's SDRF copies its samples from the task's own.");
                }

                if (!task.SearchParameters.WriteSdrf)
                {
                    task.SearchParameters.WriteSdrf = true;
                    messages.Add($"{name} writes the search's own SDRF, because an SDRF is beside the spectra files.");
                }
            }
        }
    }
}
