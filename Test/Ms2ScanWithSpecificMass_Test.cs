using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
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
        /// <summary>
        /// Generated an isotope distribution for a charge=1 peptide with sequence "PEPTIDE" using protein prospector.
        /// There were seven m/z peaks with greater than 0 intensity. These I stored in the included text file.
        /// I computed the IsotopicEnvelope from this data providing the most abundant mass as the closest.
        /// I then compared the monoisotopic mass to that computed from the PeptideWithSetMods function.
        /// </summary>
        [Test]
        public static void TestOneIsotopicEnvelope()
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

            PeptideWithSetModifications pep = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>(), 0, null, null, 1, 7, 0, CleavageSpecificity.None, null);
            Assert.That(pep.MonoisotopicMass, Is.EqualTo(closestExperimentalMass.MonoisotopicMass).Within(0.001));
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
