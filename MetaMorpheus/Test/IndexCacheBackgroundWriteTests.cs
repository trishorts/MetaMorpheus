using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Indexing;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// GenerateIndexes can leave the index cache writing on another thread while the search goes ahead with the index in memory.
    /// The cache must come out exactly as a write on the calling thread makes it. Nothing may look for the folder before the write
    /// is finished. A write that fails must still fail the task, at the point where it is waited for.
    /// </summary>
    [TestFixture]
    public static class IndexCacheBackgroundWriteTests
    {
        /// <summary>Minimal task so the index-cache methods can be exercised without running a search.</summary>
        private class IndexProbeTask : MetaMorpheusTask
        {
            public IndexProbeTask() : base(MyTask.Search) { CommonParameters = new CommonParameters(); }

            protected override MyTaskResults RunSpecific(string outputFolder, List<DbForTask> dbFilenameList,
                List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
                => throw new NotSupportedException("probe only");
        }

        /// <summary>A database folder holding a copy of the small yeast fasta, so each call caches into its own folder.</summary>
        private static (string Folder, List<DbForTask> Databases, List<Protein> Proteins) NewDatabaseFolder()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "IndexCacheBackground-" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(folder);
            string fasta = Path.Combine(folder, "smalldb.fasta");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta"), fasta);
            var proteins = ProteinDbLoader.LoadProteinFasta(fasta, true, DecoyType.Reverse, false, out _);
            return (folder, new List<DbForTask> { new DbForTask(fasta, false) }, proteins);
        }

        private static IndexingEngine MakeEngine(List<Protein> proteins, List<DbForTask> databases, bool generatePrecursorIndex)
        {
            var parameters = new CommonParameters(maxThreadsToUsePerFile: 4,
                digestionParams: new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 7));
            var variable = GlobalVariables.AllModsKnown.Where(m => m.ModificationType == "Common Variable" && m.IdWithMotif == "Oxidation on M").ToList();
            var fixedMods = GlobalVariables.AllModsKnown.Where(m => m.ModificationType == "Common Fixed" && m.IdWithMotif == "Carbamidomethyl on C").ToList();
            return new IndexingEngine(proteins, variable, fixedMods, null, null, null, 0, DecoyType.Reverse, parameters,
                new List<(string, CommonParameters)> { ("", parameters) }, 30000, generatePrecursorIndex,
                databases.Select(p => new FileInfo(p.FilePath)).ToList(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
        }

        private static string OnlyCacheFolder(string databaseFolder)
        {
            string[] folders = Directory.GetDirectories(Path.Combine(databaseFolder, "DatabaseIndex"));
            Assert.That(folders, Has.Length.EqualTo(1), "exactly one cache folder");
            return folders[0];
        }

        /// <summary>
        /// The same database indexed twice, once writing on the calling thread and once in the background: every cache file is
        /// byte for byte the same, for a fragment-only index and for one with a precursor index.
        /// </summary>
        [Test]
        public static void ABackgroundWriteLeavesTheSameCacheFilesAsAWriteOnTheCallingThread([Values(false, true)] bool generatePrecursorIndex)
        {
            var (syncFolder, syncDatabases, syncProteins) = NewDatabaseFolder();
            var (backgroundFolder, backgroundDatabases, backgroundProteins) = NewDatabaseFolder();
            try
            {
                List<IBioPolymerWithSetMods> syncPeptides = null, backgroundPeptides = null;
                FragmentIndex syncFragments = null, backgroundFragments = null;
                List<int>[] syncPrecursors = null, backgroundPrecursors = null;

                new IndexProbeTask().GenerateIndexes(MakeEngine(syncProteins, syncDatabases, generatePrecursorIndex), syncDatabases,
                    ref syncPeptides, ref syncFragments, ref syncPrecursors, syncProteins, "sync");

                var task = new IndexProbeTask();
                task.GenerateIndexes(MakeEngine(backgroundProteins, backgroundDatabases, generatePrecursorIndex), backgroundDatabases,
                    ref backgroundPeptides, ref backgroundFragments, ref backgroundPrecursors, backgroundProteins, "background", writeCacheInBackground: true);

                // the index is returned before its cache is written, and it is the one being written
                Assert.That(backgroundPeptides.Select(p => p.FullSequence), Is.EqualTo(syncPeptides.Select(p => p.FullSequence)));
                task.FinishWritingIndexCache();

                string syncCache = OnlyCacheFolder(syncFolder);
                string backgroundCache = OnlyCacheFolder(backgroundFolder);
                var names = new List<string> { MetaMorpheusTask.PeptideIndexFileName, MetaMorpheusTask.FragmentIndexFileName, "DigestionParameters.toml" };
                if (generatePrecursorIndex)
                {
                    names.Add(MetaMorpheusTask.PrecursorIndexFileName);
                }

                Assert.That(Directory.GetFiles(backgroundCache).Select(Path.GetFileName).OrderBy(n => n),
                    Is.EqualTo(Directory.GetFiles(syncCache).Select(Path.GetFileName).OrderBy(n => n)), "the same files, and no partial file left");
                foreach (string name in names)
                {
                    Assert.That(File.ReadAllBytes(Path.Combine(backgroundCache, name)), Is.EqualTo(File.ReadAllBytes(Path.Combine(syncCache, name))), name);
                }
            }
            finally
            {
                Directory.Delete(syncFolder, true);
                Directory.Delete(backgroundFolder, true);
            }
        }

        /// <summary>
        /// The next GenerateIndexes waits for the write before it looks for the folder. So it finds the cache the write left,
        /// reads it, and does not build and write the same index a second time.
        /// </summary>
        [Test]
        public static void TheNextLookupWaitsForTheWriteAndReadsWhatItWrote()
        {
            var (folder, databases, proteins) = NewDatabaseFolder();
            try
            {
                var task = new IndexProbeTask();
                List<IBioPolymerWithSetMods> built = null, read = null;
                FragmentIndex builtFragments = null, readFragments = null;
                List<int>[] builtPrecursors = null, readPrecursors = null;

                task.GenerateIndexes(MakeEngine(proteins, databases, false), databases, ref built, ref builtFragments, ref builtPrecursors, proteins, "first",
                    writeCacheInBackground: true);
                task.GenerateIndexes(MakeEngine(proteins, databases, false), databases, ref read, ref readFragments, ref readPrecursors, proteins, "second",
                    writeCacheInBackground: true);
                task.FinishWritingIndexCache();

                OnlyCacheFolder(folder);
                Assert.That(read, Is.Not.SameAs(built), "the second call must read the cache, not build again");
                Assert.That(read.Select(p => p.FullSequence), Is.EqualTo(built.Select(p => p.FullSequence)));
                Assert.That(readFragments.PeptideIds, Is.EqualTo(builtFragments.PeptideIds));
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// A write that fails in the background fails the task where the write is waited for: the next GenerateIndexes, before
        /// it looks for or builds anything, or the explicit wait. It fails with the exception itself, not one wrapped around it, and
        /// only once. The real writer still cleans up its partial file.
        /// </summary>
        [Test]
        public static void AFailedBackgroundWriteSurfacesAtTheWait()
        {
            var (folder, databases, proteins) = NewDatabaseFolder();
            try
            {
                var task = new IndexProbeTask();
                string target = Path.Combine(folder, MetaMorpheusTask.FragmentIndexFileName);
                var writeThroughTemporaryFile = typeof(MetaMorpheusTask).GetMethod("WriteThroughTemporaryFile", BindingFlags.NonPublic | BindingFlags.Static);

                // the real writer, failing part way as a full disk would
                task.StartWritingIndexCache(() =>
                {
                    try
                    {
                        writeThroughTemporaryFile.Invoke(null, new object[] { target, (Action<FileStream>)(file =>
                        {
                            file.Write(new byte[64]);
                            throw new IOException("disk full");
                        }) });
                    }
                    catch (TargetInvocationException e)
                    {
                        System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(e.InnerException).Throw();
                    }
                });

                List<IBioPolymerWithSetMods> peptides = null;
                FragmentIndex fragments = null;
                List<int>[] precursors = null;
                var thrown = Assert.Throws<IOException>(() => task.GenerateIndexes(MakeEngine(proteins, databases, false), databases,
                    ref peptides, ref fragments, ref precursors, proteins, "next"));
                Assert.That(thrown.Message, Is.EqualTo("disk full"));
                Assert.That(Directory.Exists(Path.Combine(folder, "DatabaseIndex")), Is.False, "the failure surfaces before anything is looked for or built");
                Assert.That(File.Exists(target) || File.Exists(target + ".partial"), Is.False, "the failed write leaves no file");

                Assert.DoesNotThrow(() => task.FinishWritingIndexCache(), "a failure is reported once");

                task.StartWritingIndexCache(() => throw new InvalidOperationException("second failure"));
                Assert.That(() => task.FinishWritingIndexCache(), Throws.InstanceOf<InvalidOperationException>().With.Message.EqualTo("second failure"));
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }
    }
}
