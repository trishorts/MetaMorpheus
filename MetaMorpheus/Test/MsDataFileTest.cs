using Easy.Common.Extensions;
using EngineLayer;
using MassSpectrometry;
using NetSerializer;
using Nett;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MathNet.Numerics;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class MsDataFileTest
    {
        [OneTimeSetUp]
        public static void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
        }

        [Test]
        public static void TestLoadAndRunMgf()
        {
            //The purpose of this test is to ensure that mgfs can be run without crashing.
            //Whenever a new feature is added that may require things an mgf does not have,
            //there should be a check that prevents mgfs from using that feature.
            string mgfName = @"TestData\ok.mgf";
            string xmlName = @"TestData\okk.xml";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestLoadAndRunMgf");

            SearchTask task1 = new()
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DoLabelFreeQuantification = true
                }
            };
            List<(string, MetaMorpheusTask)> taskList = new()
            {
                ("task1", task1),
            };
            //run!

            var engine = new EverythingRunnerEngine(taskList, new List<string> { mgfName }, new List<DbForTask> { new DbForTask(xmlName, false) }, outputFolder);
            engine.Run();
            //Just don't crash! There should also be at least one psm at 1% FDR, but can't check for that.
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestCompressionDecompression()
        {
            string testInputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"CompressionTest");
            DirectoryInfo testDirectory = new(testInputFolder);
            MyFileManager.CompressDirectory(testDirectory);

            foreach (FileInfo file in testDirectory.GetFiles())
            {
                Assert.AreEqual(".gz", file.Extension);
            }

            MyFileManager.DecompressDirectory(testDirectory);

            foreach (FileInfo file in testDirectory.GetFiles())
            {
                Assert.AreNotEqual(".gz", file.Extension);
            }
        }

        [Test]
        public static void TestMs2ScanWithSpecificMass()
        {
            Ms2ScanWithSpecificMass scanB = new(
                new MsDataScan(
                    new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);

            var closestExperimentalMassB = scanB.GetClosestExperimentalIsotopicEnvelope(10);

            Assert.IsNull(closestExperimentalMassB);
        }

        [Test]
        public static void TestMs2ScanWithSpecificMassHiResVersusLowRes()
        {
            //load hi res data
            MyFileManager myFileManager = new MyFileManager(true);
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\sliced_b6.mzML");
            MsDataFile myFile = myFileManager.LoadFile(origDataFile, new CommonParameters());

            var highMs2Scans = SearchTask.GetMs2Scans(myFile,origDataFile,new CommonParameters());

            Assert.AreEqual(447, highMs2Scans.Count());

            List<double> highMs2PrecursorMasses = highMs2Scans.Select(b => b.PrecursorMass.Round(3)).ToList();

            //switch to LowCID
            MsDataFile lowFile = myFileManager.LoadFile(origDataFile, new CommonParameters(dissociationType:DissociationType.LowCID));

            var lowMs2Scans = SearchTask.GetMs2Scans(lowFile, origDataFile, new CommonParameters(dissociationType: DissociationType.LowCID));

            Assert.AreEqual(447, lowMs2Scans.Count());

            List<double> lowMs2PrecursorMasses = lowMs2Scans.Select(b => b.PrecursorMass.Round(3)).ToList();

            CollectionAssert.AreEquivalent(highMs2PrecursorMasses, lowMs2PrecursorMasses);
        }


        [Test]
        public static void TestMs2ScanWithSpecificMassHiResVersusLowResCompIons()
        {
            //load hi res data
            MyFileManager myFileManager = new MyFileManager(true);
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\sliced_b6.mzML");
            MsDataFile myFile = myFileManager.LoadFile(origDataFile, new CommonParameters(addCompIons: true));

            var highMs2Scans = SearchTask.GetMs2Scans(myFile, origDataFile, new CommonParameters(addCompIons: true));

            Assert.AreEqual(447, highMs2Scans.Count());

            List<double> highMs2PrecursorMasses = highMs2Scans.Select(b => b.PrecursorMass.Round(3)).ToList();

            //switch to LowCID
            MsDataFile lowFile = myFileManager.LoadFile(origDataFile, new CommonParameters(dissociationType: DissociationType.LowCID, addCompIons: true));

            var lowMs2Scans = SearchTask.GetMs2Scans(lowFile, origDataFile, new CommonParameters(dissociationType: DissociationType.LowCID, addCompIons: true));

            Assert.AreEqual(447, lowMs2Scans.Count());

            List<double> lowMs2PrecursorMasses = lowMs2Scans.Select(b => b.PrecursorMass.Round(3)).ToList();

            CollectionAssert.AreEquivalent(highMs2PrecursorMasses, lowMs2PrecursorMasses);
        }
    }
}