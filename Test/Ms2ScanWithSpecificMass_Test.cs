using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class Ms2ScanWithSpecificMass_Test
    {
        [Test]
        public static void MyTest()
        {
            string spectrumFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Ms2ScanTestData\PEPTIDE_MS1_Z1.txt");
            MzSpectrum mzs = ReadMzAndIntensity(spectrumFilePath);
            MzRange range = mzs.Range;
            double newRetentionTime = 1;
            int newOneBasedScanNumber = 1;
            bool isCentroid = true;
            double totalIonCurrent = 1;
            double injectionTime = 1;
            double[,] noiseData = new double[,] { { 0,0} };
            double precursorMonoisotopicPeakMz = 1;
            int precursorCharge = 1;

            MsDataScan msd = new MsDataScan(mzs, newOneBasedScanNumber, 2, isCentroid, Polarity.Positive, newRetentionTime, range, "", MZAnalyzerType.Orbitrap, totalIonCurrent, injectionTime, noiseData, "", null, null, null, null, null, DissociationType.HCD, null, null);
            Ms2ScanWithSpecificMass mwsm = new Ms2ScanWithSpecificMass(msd, precursorMonoisotopicPeakMz, precursorCharge, spectrumFilePath, new CommonParameters(), null);
            IsotopicEnvelope closestExperimentalMass = mwsm.GetClosestExperimentalIsotopicEnvelope(range.Minimum);




            Assert.AreEqual(1, 2);
        }

        public static MzSpectrum ReadMzAndIntensity(string file)
        {
            string[] mzIntensityInput = File.ReadAllLines(file);
            double[] mzValues = new double[mzIntensityInput.Length];
            double[] intensities = new double[mzIntensityInput.Length];
            for (int i = 0; i < mzIntensityInput.Length; i++)
            {
                string[] pair = mzIntensityInput[i].Split('\t');
                mzValues[i] = Convert.ToDouble(pair[0]);
                intensities[i] = Convert.ToDouble(pair[1]);
            }
            return new MzSpectrum(mzValues, intensities, true);
        }
    }
}
