using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using IsobaricMassTagType = Omics.Modifications.IsobaricMassTagType;

namespace Test
{
    /// <summary>
    /// The GUI's half of the automatic design (QuantProject, ruled 2026-11-13), as far as it lives outside GUI.csproj, which
    /// this project cannot reference: the design made when spectra files are added (kit guessed from the SDRF), made again at
    /// Run only when it no longer fits, the task settings put back after a run, dropped SDRF and design files, and the
    /// label-free design window's Save keeping a Sample Type column and keeping the design ours.
    /// </summary>
    [TestFixture]
    internal static class GuiSdrfAutomaticDesignTests
    {
        private static string NewFolder(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "GuiSdrfAutomaticDesign", name);
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            return folder;
        }

        private static List<string> Spectra(string folder, params string[] names) =>
            names.Select(n =>
            {
                string path = Path.Combine(folder, n);
                File.WriteAllText(path, string.Empty);
                return path;
            }).ToList();

        private static string WriteSdrf(string folder, string[] header, params string[][] rows)
        {
            string path = Path.Combine(folder, "deposit.sdrf.tsv");
            File.WriteAllLines(path, new[] { string.Join("\t", header) }.Concat(rows.Select(r => string.Join("\t", r))));
            return path;
        }

        private static readonly string[] LabelFreeHeader =
            { "source name", "characteristics[biological replicate]", "comment[label]", "comment[data file]", "factor value[genotype]" };

        private static string WriteLabelFreeSdrf(string folder) => WriteSdrf(folder, LabelFreeHeader,
            new[] { "WT 1", "1", "label free sample", "a1.raw", "WT" },
            new[] { "KO 1", "1", "label free sample", "b1.raw", "KO" });

        private static readonly string[] TmtHeader =
            { "source name", "characteristics[biological replicate]", "comment[label]", "comment[data file]", "comment[plex]", "factor value[genotype]" };

        /// <summary>One plex of the given labels, alternating WT and KO, each its own biological replicate within its group.</summary>
        private static string WriteTmtSdrf(string folder, params string[] labels) => WriteSdrf(folder, TmtHeader,
            labels.Select((l, i) => new[] { "S" + i, (i / 2 + 1).ToString(), l, "p1.raw", "P1", i % 2 == 0 ? "WT" : "KO" }).ToArray());

        private static string LabelFreeDesign(string folder) => Path.Combine(folder, GlobalVariables.ExperimentalDesignFileName);
        private static string TmtDesign(string folder) => Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);

        private static SearchTask Tmt(string kit)
        {
            var task = new SearchTask();
            task.SearchParameters.DoMultiplexQuantification = true;
            task.SearchParameters.MultiplexModId = kit;
            return task;
        }

        // ---------------------------------------------------------------- the kit, guessed from the SDRF's labels

        [TestCase(new[] { "TMT126", "TMT127N", "TMT131" }, IsobaricMassTagType.TMT10, TestName = "TMT10's own 131 spelling is TMT10")]
        [TestCase(new[] { "TMT126", "TMT127N", "TMT131C" }, IsobaricMassTagType.TMT11, TestName = "131C makes it TMT11")]
        [TestCase(new[] { "TMT126", "TMT134N" }, IsobaricMassTagType.TMT16, TestName = "134N fits TMT16, the smaller TMTpro kit")]
        [TestCase(new[] { "TMT126", "TMT135N" }, IsobaricMassTagType.TMT18, TestName = "135N makes it TMT18")]
        public static void TheKitIsTheSmallestThatHoldsEveryLabel(string[] labels, IsobaricMassTagType expected)
        {
            string folder = NewFolder("Kit" + expected);
            string sdrf = WriteTmtSdrf(folder, labels);

            Assert.That(SdrfAutomaticDesign.GuessIsobaricKit(sdrf, out var tag, out string undecided), Is.True, undecided);
            Assert.That(tag.TagType, Is.EqualTo(expected));
        }

        [Test]
        public static void ALabelFreeSdrfHasNoKit()
        {
            string folder = NewFolder("KitLabelFree");
            Assert.That(SdrfAutomaticDesign.GuessIsobaricKit(WriteLabelFreeSdrf(folder), out var tag, out _), Is.True);
            Assert.That(tag, Is.Null);
        }

        [TestCase("label free sample", "TMT126", "mixes label-free and labelled rows")]
        [TestCase("TMT126", "XYZ", "No isobaric kit holds every label")]
        public static void AKitThatCannotBeToldIsNotGuessed(string first, string second, string why)
        {
            string folder = NewFolder("KitUndecided" + second);
            string sdrf = WriteTmtSdrf(folder, first, second);

            Assert.That(SdrfAutomaticDesign.GuessIsobaricKit(sdrf, out var tag, out string undecided), Is.False);
            Assert.That(tag, Is.Null);
            Assert.That(undecided, Does.Contain(why));
        }

        // ---------------------------------------------------------------- the design, made when the files are added

        /// <summary>The user's example: the files are added, and the design window opens filled in, before Run.</summary>
        [Test]
        public static void AddingFilesBesideALabelFreeSdrfWritesTheDesignAndItsMarker()
        {
            string folder = NewFolder("AddLabelFree");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteLabelFreeSdrf(folder);

            var outcome = SdrfAutomaticDesign.ResolveOnAddingFiles(folder, spectra, null);

            Assert.That(outcome.Refusal, Is.Null, string.Join("\n", outcome.Messages));
            var design = ExperimentalDesign.ReadExperimentalDesign(LabelFreeDesign(folder), spectra, out var errors);
            Assert.That(errors, Is.Empty);
            Assert.That(design.Select(f => f.Condition), Is.EquivalentTo(new[] { "WT", "KO" }));
            Assert.That(SdrfDesignMarker.ReadFrom(folder).DescribesUneditedDesignIn(folder), Is.True);
            Assert.That(outcome.Messages.Last(), Does.StartWith("Wrote " + LabelFreeDesign(folder)));
        }

        [Test]
        public static void AddingFilesBesideATmtSdrfWritesTmtDesignWithTheGuessedKitAndCommentPlex()
        {
            string folder = NewFolder("AddTmt");
            var spectra = Spectra(folder, "p1.mzML");
            WriteTmtSdrf(folder, "TMT126", "TMT127N", "TMT127C", "TMT128N");

            var outcome = SdrfAutomaticDesign.ResolveOnAddingFiles(folder, spectra, null);

            Assert.That(outcome.Refusal, Is.Null, string.Join("\n", outcome.Messages));
            var files = TmtExperimentalDesign.Read(TmtDesign(folder), spectra, out var errors);
            Assert.That(errors, Is.Empty);
            Assert.That(files.Single().Plex, Is.EqualTo("P1"));
            Assert.That(File.Exists(LabelFreeDesign(folder)), Is.False);
        }

        [Test]
        public static void AddingFilesBesideAUsersDesignLeavesItAlone()
        {
            string folder = NewFolder("AddUsers");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteLabelFreeSdrf(folder);
            File.WriteAllLines(LabelFreeDesign(folder), new[]
            {
                "FileName\tCondition\tBiorep\tFraction\tTechrep", "a1.mzML\tmine\t1\t1\t1", "b1.mzML\tmine\t2\t1\t1",
            });
            string users = File.ReadAllText(LabelFreeDesign(folder));

            var outcome = SdrfAutomaticDesign.ResolveOnAddingFiles(folder, spectra, null);

            Assert.That(File.ReadAllText(LabelFreeDesign(folder)), Is.EqualTo(users));
            Assert.That(File.Exists(SdrfDesignMarker.PathIn(folder)), Is.False);
            Assert.That(outcome.Messages.Single(), Does.Contain("so it is the design"));
        }

        [Test]
        public static void AKitThatCannotBeGuessedLeavesTheDesignToRun()
        {
            string folder = NewFolder("AddUndecided");
            var spectra = Spectra(folder, "p1.mzML");
            WriteTmtSdrf(folder, "TMT126", "XYZ");

            var outcome = SdrfAutomaticDesign.ResolveOnAddingFiles(folder, spectra, null);

            Assert.That(outcome.Refusal, Is.Null);
            Assert.That(outcome.Messages.Single(), Does.EndWith("The experimental design is made from the SDRF at Run, with the search task's label."));
            Assert.That(File.Exists(TmtDesign(folder)) || File.Exists(LabelFreeDesign(folder)), Is.False);
        }

        [Test]
        public static void AnSdrfThatGivesNoDesignIsAReasonNotADesign()
        {
            string folder = NewFolder("AddRefused");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, LabelFreeHeader.Take(4).ToArray(), new[] { "WT 1", "1", "label free sample", "a1.raw" });

            var outcome = SdrfAutomaticDesign.ResolveOnAddingFiles(folder, spectra, null);

            Assert.That(outcome.Refusal, Does.Contain("no experimental design could be made from it"));
            Assert.That(File.Exists(LabelFreeDesign(folder)), Is.False);
        }

        // ---------------------------------------------------------------- at Run: the task's kit decides

        /// <summary>
        /// Ruled 2026-11-13: the SDRF names only ten channels, so TMT10 was guessed when the files were added; the search is
        /// TMT11. Every channel the SDRF names is in TMT11, so the run goes ahead, its design made with the task's kit.
        /// </summary>
        [Test]
        public static void AtRunATaskKitHoldingEveryNamedChannelGoesAhead()
        {
            string folder = NewFolder("RunTmt11");
            var spectra = Spectra(folder, "p1.mzML");
            WriteTmtSdrf(folder, "TMT126", "TMT127N", "TMT131N");
            Assume.That(SdrfAutomaticDesign.GuessIsobaricKit(Path.Combine(folder, "deposit.sdrf.tsv"), out var guessed, out _)
                && guessed.TagType == IsobaricMassTagType.TMT10, "ten channels' labels guess TMT10");
            SdrfAutomaticDesign.ResolveOnAddingFiles(folder, spectra, null);
            Assume.That(File.Exists(TmtDesign(folder)));

            var outcome = SdrfAutomaticDesign.Resolve(folder, spectra, new List<(string, SearchTask)> { ("Task1", Tmt("TMT11")) }, null, null);

            Assert.That(outcome.Refusal, Is.Null, string.Join("\n", outcome.Messages));
            TmtExperimentalDesign.Read(TmtDesign(folder), spectra, out var errors);
            Assert.That(errors, Is.Empty);
        }

        /// <summary>Ruled 2026-11-13: refused at Run when the SDRF names a channel the task's kit does not have.</summary>
        [Test]
        public static void AtRunAChannelTheTaskKitLacksIsRefusedAndOurDesignRemoved()
        {
            string folder = NewFolder("RunTmt6");
            var spectra = Spectra(folder, "p1.mzML");
            WriteTmtSdrf(folder, "TMT126", "TMT127N", "TMT127C");
            SdrfAutomaticDesign.ResolveOnAddingFiles(folder, spectra, null);
            Assume.That(File.Exists(TmtDesign(folder)));

            var outcome = SdrfAutomaticDesign.Resolve(folder, spectra, new List<(string, SearchTask)> { ("Task1", Tmt("TMT6")) }, null, null);

            Assert.That(outcome.Refusal, Is.Not.Null);
            Assert.That(File.Exists(TmtDesign(folder)), Is.False);
            Assert.That(outcome.Messages, Has.Some.EndsWith("was removed: MetaMorpheus wrote it from the SDRF before, and it no longer fits."));
        }

        /// <summary>Ruled 2026-11-13: at Run, our design is made again only when it no longer fits; otherwise it is left as is.</summary>
        [Test]
        public static void AtRunOurDesignIsLeftAsItIsWhenItStillFits()
        {
            string folder = NewFolder("RunFits");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteLabelFreeSdrf(folder);
            SdrfAutomaticDesign.ResolveOnAddingFiles(folder, spectra, null);
            var written = File.GetLastWriteTimeUtc(LabelFreeDesign(folder));

            var outcome = SdrfAutomaticDesign.Resolve(folder, spectra, new List<(string, SearchTask)> { ("Task1", new SearchTask()) }, null, null);

            Assert.That(outcome.Messages.Last(), Does.EndWith("already fits this run's spectra files."));
            Assert.That(File.GetLastWriteTimeUtc(LabelFreeDesign(folder)), Is.EqualTo(written));
            Assert.That(Directory.GetFiles(folder, "*.new"), Is.Empty);
        }

        // ---------------------------------------------------------------- task settings, for one run only

        [Test]
        public static void TheTaskSettingsARunChangesComeBack()
        {
            string folder = NewFolder("Snapshot");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteLabelFreeSdrf(folder);
            var task = new SearchTask();
            task.SearchParameters.SdrfConditionColumns = new List<string> { "genotype" };
            var snapshot = SdrfTaskSettingsSnapshot.Take(new[] { task });

            SdrfAutomaticDesign.Resolve(folder, spectra, new List<(string, SearchTask)> { ("Task1", task) }, null, null);
            Assume.That(task.SearchParameters.InputSdrfPath, Is.Not.Null);
            Assume.That(task.SearchParameters.WriteSdrf, Is.True);
            Assume.That(task.SearchParameters.DesignWrittenFromSdrf, Is.Not.Null);

            snapshot.Restore();

            Assert.That(task.SearchParameters.InputSdrfPath, Is.Null);
            Assert.That(task.SearchParameters.WriteSdrf, Is.False);
            Assert.That(task.SearchParameters.DesignWrittenFromSdrf, Is.Null);
            Assert.That(task.SearchParameters.SdrfConditionColumns, Is.EqualTo(new[] { "genotype" }));
        }

        // ---------------------------------------------------------------- dropped files

        [TestCase("deposit.sdrf.tsv", true)]
        [TestCase("sdrf.tsv", true)]
        [TestCase("ExperimentalDesign.tsv", true)]
        [TestCase("TmtDesign.txt", true)]
        [TestCase("DesignFromSdrf.txt", true)]
        [TestCase("results.tsv", false)]
        public static void DroppedSdrfAndDesignFilesAreRecognised(string name, bool recognised)
        {
            Assert.That(SdrfAutomaticDesign.IsDesignOrSdrfFile(Path.Combine("C:", "data", name)), Is.EqualTo(recognised));
        }

        [Test]
        public static void ADroppedSdrfIsUsedOnlyBesideTheSpectra()
        {
            string spectraFolder = NewFolder("DropSpectra");
            string elsewhere = NewFolder("DropElsewhere");

            Assert.That(SdrfAutomaticDesign.NoteForDroppedFile(Path.Combine(spectraFolder, "x.sdrf.tsv"), spectraFolder), Is.Null);
            Assert.That(SdrfAutomaticDesign.NoteForDroppedFile(Path.Combine(elsewhere, "x.sdrf.tsv"), spectraFolder),
                Does.StartWith("x.sdrf.tsv is not used: SDRF and experimental design files are read from the folder of the first spectra file")
                    .And.EndWith("and it is in " + elsewhere + "."));
            Assert.That(SdrfAutomaticDesign.NoteForDroppedFile(Path.Combine(elsewhere, "x.sdrf.tsv"), null), Does.EndWith("no spectra file is loaded."));
        }

        // ---------------------------------------------------------------- the label-free design window's Save

        /// <summary>
        /// The window saves through ExperimentalDesign.WriteExperimentalDesignToFile, which has no Sample Type in its rows: a
        /// column the file has is kept by file name, so a pooled reference stays marked.
        /// </summary>
        [Test]
        public static void SavingKeepsTheSampleTypeColumn()
        {
            string folder = NewFolder("SampleType");
            var spectra = Spectra(folder, "a1.mzML", "pool.mzML");
            File.WriteAllLines(LabelFreeDesign(folder), new[]
            {
                "FileName\tCondition\tBiorep\tFraction\tTechrep\tSample Type", "a1.mzML\tWT\t1\t1\t1\tstudy sample", "pool.mzML\tWT\t2\t1\t1\treference",
            });
            var design = ExperimentalDesign.ReadExperimentalDesign(LabelFreeDesign(folder), spectra, out var errors);
            Assume.That(errors, Is.Empty);

            ExperimentalDesign.WriteExperimentalDesignToFile(design);

            Assert.That(ExperimentalDesign.ReadSampleTypes(LabelFreeDesign(folder)),
                Is.EquivalentTo(new Dictionary<string, string> { ["a1.mzML"] = "study sample", ["pool.mzML"] = "reference" }));
        }

        [Test]
        public static void SavingWithoutASampleTypeColumnAddsNone()
        {
            string folder = NewFolder("NoSampleType");
            var spectra = Spectra(folder, "a1.mzML");

            ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo> { new(spectra[0], "WT", 0, 0, 0) });

            Assert.That(File.ReadAllLines(LabelFreeDesign(folder))[0], Is.EqualTo("FileName\tCondition\tBiorep\tFraction\tTechrep"));
        }

        /// <summary>
        /// Ruled 2026-11-13: only a real edit makes a design the user's. Our design, opened in the window and saved unchanged
        /// (read, then written back in the spectra files' order, with Windows line endings), is still ours.
        /// </summary>
        [Test]
        public static void OurDesignSavedUnchangedByTheWindowIsStillOurs()
        {
            string folder = NewFolder("WindowSave");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteLabelFreeSdrf(folder);
            SdrfAutomaticDesign.ResolveOnAddingFiles(folder, spectra, null);
            string before = SdrfDesignMarker.Sha256Of(LabelFreeDesign(folder));

            var design = ExperimentalDesign.ReadExperimentalDesign(LabelFreeDesign(folder), spectra, out _);
            ExperimentalDesign.WriteExperimentalDesignToFile(spectra.Select(s => design.Single(d => d.FullFilePathWithExtension == s)).ToList());
            Assume.That(SdrfDesignMarker.Sha256Of(LabelFreeDesign(folder)), Is.Not.EqualTo(before), "the bytes changed");

            Assert.That(SdrfDesignMarker.ReadFrom(folder).DescribesUneditedDesignIn(folder), Is.True);
        }
    }
}
