using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using StatisticalModels;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// The PEP column reports a posterior error probability read off the target-decoy distribution of the classifier's
    /// PEP (<see cref="FdrInfo.EstimatedPep"/>), at PSM level and peptide level separately. The classifier's value
    /// (<see cref="FdrInfo.PEP"/>) stays what every ranking reads, so the estimate must be a function of it alone within
    /// each population, never fall as it rises, and be left out where PEP was not trained.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class PepEstimateTests
    {
        private const string DataFileName = "TaGe_SA_HeLa_04_subset_longestSeq.mzML";

        private static PeptideWithSetModifications Peptide(string sequence, bool isDecoy = false) =>
            new(sequence, [], p: new Protein(sequence, (isDecoy ? "DECOY_" : "") + sequence, isDecoy: isDecoy));

        private static SpectralMatch Match(int scanNumber, double classifierPep, PeptideWithSetModifications peptide)
        {
            var scan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), scanNumber, 2, true,
                Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=" + scanNumber,
                double.NaN, null, null, double.NaN, null, DissociationType.HCD, 1, null);
            var ms2 = new Ms2ScanWithSpecificMass(scan, peptide.MonoisotopicMass.ToMz(1), 1, "A.mzML", new CommonParameters());
            var match = new PeptideSpectralMatch(peptide, 0, 10, scanNumber, ms2, new CommonParameters(), new List<MatchedFragmentIon>());
            match.PsmFdrInfo = new FdrInfo { PEP = classifierPep };
            match.PeptideFdrInfo = new FdrInfo { PEP = classifierPep };
            return match;
        }

        /// <summary>A decoy by its hypotheses, as the FDR code counts it (SpectralMatch.IsDecoy is set only on resolving).</summary>
        private static bool IsDecoyMatch(SpectralMatch m) => m.BestMatchingBioPolymersWithSetMods.Any(h => h.IsDecoy);

        /// <summary>Twelve matches: classifier PEPs rising, decoys mixed in toward the bad end, one tie.</summary>
        private static List<SpectralMatch> TwelveMatches()
        {
            double[] pep = { 0.001, 0.002, 0.004, 0.01, 0.02, 0.05, 0.05, 0.2, 0.4, 0.6, 0.8, 0.95 };
            bool[] decoy = { false, false, false, false, true, false, false, false, true, false, true, true };
            return Enumerable.Range(0, pep.Length).Select(i => Match(i + 1, pep[i], Peptide("PEPTIDE" + (char)('A' + i) + "K", decoy[i]))).ToList();
        }

        [Test]
        public static void TheEstimateIsTargetDecoyPepOfTheClassifierPepAndIsWrittenAtOneLevelOnly()
        {
            var matches = TwelveMatches();
            FdrAnalysisEngine.EstimatePep(matches, peptideLevel: false);

            var expected = TargetDecoyPep.Estimate(matches.Select(m => -m.PsmFdrInfo.PEP).ToArray(), matches.Select(IsDecoyMatch).ToArray()).Peps;
            for (int i = 0; i < matches.Count; i++)
            {
                Assert.That(matches[i].PsmFdrInfo.EstimatedPep, Is.EqualTo(expected[i]), $"match {i}");
                Assert.That(matches[i].PeptideFdrInfo.EstimatedPep, Is.Null, "the peptide level is a different population and is left alone");
                Assert.That(matches[i].PsmFdrInfo.PEP, Is.EqualTo(new[] { 0.001, 0.002, 0.004, 0.01, 0.02, 0.05, 0.05, 0.2, 0.4, 0.6, 0.8, 0.95 }[i]),
                    "the classifier's value, which every ranking reads, is not touched");
            }
            Assert.That(matches[0].PsmFdrInfo.EstimatedPep, Is.EqualTo(0), "above every decoy");
            Assert.That(matches[5].PsmFdrInfo.EstimatedPep, Is.EqualTo(matches[6].PsmFdrInfo.EstimatedPep), "a tied classifier PEP gets one estimate");
            for (int i = 1; i < matches.Count; i++)
            {
                Assert.That(matches[i].PsmFdrInfo.EstimatedPep, Is.GreaterThanOrEqualTo(matches[i - 1].PsmFdrInfo.EstimatedPep),
                    "the estimate never falls as the classifier PEP rises, so no PEP-ordered result can change");
            }
        }

        [Test]
        public static void AMatchHalfTargetHalfDecoyCountsAsHalfADecoy()
        {
            var matches = TwelveMatches();
            matches[1].AddOrReplace(Peptide("DECOYPEPTIDEK", isDecoy: true), 10, 0, true, new List<MatchedFragmentIon>());
            Assert.That(matches[1].BestMatchingBioPolymersWithSetMods.Count(), Is.EqualTo(2));

            FdrAnalysisEngine.EstimatePep(matches, peptideLevel: true);

            double[] fraction = matches.Select(m => m.BestMatchingBioPolymersWithSetMods.Count(h => h.IsDecoy) / (double)m.BestMatchingBioPolymersWithSetMods.Count()).ToArray();
            Assert.That(fraction[1], Is.EqualTo(0.5));
            var expected = TargetDecoyPep.Estimate(matches.Select(m => -m.FdrInfo.PEP).ToArray(), fraction).Peps;
            for (int i = 0; i < matches.Count; i++)
            {
                Assert.That(matches[i].PeptideFdrInfo.EstimatedPep, Is.EqualTo(expected[i]), $"match {i}");
                Assert.That(matches[i].PsmFdrInfo.EstimatedPep, Is.Null);
            }
            Assert.That(matches[1].PeptideFdrInfo.EstimatedPep, Is.GreaterThan(0), "the half decoy is the first decoy weight in the ranking");
        }

        [Test]
        public static void NothingIsEstimatedWhenPepWasNotTrainedOrThereAreNoDecoys()
        {
            // PEP not trained: every classifier PEP is the default 0.
            var untrained = TwelveMatches();
            untrained.ForEach(m => m.PsmFdrInfo.PEP = 0);
            FdrAnalysisEngine.EstimatePep(untrained, peptideLevel: false);
            Assert.That(untrained.All(m => m.PsmFdrInfo.EstimatedPep == null), Is.True, "one classifier value says nothing about any match");

            var targetsOnly = TwelveMatches().Where(m => !IsDecoyMatch(m)).ToList();
            FdrAnalysisEngine.EstimatePep(targetsOnly, peptideLevel: false);
            Assert.That(targetsOnly.All(m => m.PsmFdrInfo.EstimatedPep == null), Is.True, "no decoys, no estimate");

            // An earlier pass's estimate is kept, not erased, when a later pass cannot estimate.
            var kept = TwelveMatches();
            FdrAnalysisEngine.EstimatePep(kept, peptideLevel: false);
            double? before = kept[3].PsmFdrInfo.EstimatedPep;
            var noDecoysNow = kept.Where(m => !IsDecoyMatch(m)).ToList();
            FdrAnalysisEngine.EstimatePep(noDecoysNow, peptideLevel: false);
            Assert.That(kept[3].PsmFdrInfo.EstimatedPep, Is.EqualTo(before));
        }

        [Test]
        public static void ThePepColumnReportsTheEstimateAndFallsBackToTheClassifier()
        {
            var match = TwelveMatches()[3];
            match.ResolveAllAmbiguities();
            match.PsmFdrInfo.EstimatedPep = 0.25;
            var columns = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchScoreData(columns, match, writePeptideLevelFdr: false);
            Assert.That(columns[SpectrumMatchFromTsvHeader.PEP], Is.EqualTo(0.25.ToString()));

            columns.Clear();
            PsmTsvWriter.AddMatchScoreData(columns, match, writePeptideLevelFdr: true);
            Assert.That(columns[SpectrumMatchFromTsvHeader.PEP], Is.EqualTo(0.01.ToString()), "no peptide-level estimate: the classifier's value, as before");
        }

        [Test]
        public static void OnARealSearchEachLevelIsEstimatedOnItsOwnPopulation()
        {
            var commonParameters = new CommonParameters(digestionParams: new DigestionParams());
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)> { (DataFileName, commonParameters) };
            var dataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData", DataFileName);
            var msDataFile = new MyFileManager(true).LoadFile(dataFile, commonParameters);
            List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta"),
                true, DecoyType.Reverse, false, out _, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex,
                ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex, -1);
            var scans = MetaMorpheusTask.GetMs2Scans(msDataFile, dataFile, commonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] psmArray = new PeptideSpectralMatch[scans.Length];
            new ClassicSearchEngine(psmArray, scans, new List<Modification>(), new List<Modification>(), null, null, null,
                proteins, new SinglePpmAroundZeroSearchMode(5), commonParameters, fsp, null, new List<string>(), false).Run();
            var psms = psmArray.Where(p => p != null).ToList();
            new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>()).Run();

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPepEstimate");
            Directory.CreateDirectory(outputFolder);
            new PepAnalysisEngine(psms, "standard", fsp, outputFolder) { UsePeptideLevelQValueForTraining = true }.ComputePEPValuesForAllPSMs();
            Assert.That(psms.Select(p => p.PsmFdrInfo.PEP).Distinct().Count(), Is.GreaterThan(1), "PEP trained");

            // PEPs exist and are not retrained: the PEP q-value passes run at both levels, and the estimate with them.
            new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>(), doPEP: false).Run();

            static double DecoyShare(SpectralMatch m) =>
                m.BestMatchingBioPolymersWithSetMods.Count(h => h.IsDecoy) / (double)m.BestMatchingBioPolymersWithSetMods.Count();

            var psmExpected = TargetDecoyPep.Estimate(psms.Select(p => -p.FdrInfo.PEP).ToArray(), psms.Select(DecoyShare).ToArray()).Peps;
            for (int i = 0; i < psms.Count; i++)
                Assert.That(psms[i].PsmFdrInfo.EstimatedPep, Is.EqualTo(psmExpected[i]), $"PSM {i}");

            // The peptide population: the lowest-PEP match of each full sequence, as the PEP q-value pass takes it.
            var peptides = psms.OrderBy(p => p.FdrInfo.PEP).ThenByDescending(p => p).GroupBy(p => p.FullSequence).Select(g => g.First()).ToList();
            var peptideExpected = TargetDecoyPep.Estimate(peptides.Select(p => -p.FdrInfo.PEP).ToArray(), peptides.Select(DecoyShare).ToArray()).Peps;
            for (int i = 0; i < peptides.Count; i++)
                Assert.That(peptides[i].PeptideFdrInfo.EstimatedPep, Is.EqualTo(peptideExpected[i]), $"peptide {i}");

            Assert.That(peptides.Zip(peptideExpected).Any(t => t.Second != psmExpected[psms.IndexOf(t.First)]), Is.True,
                "the two populations differ, so the two levels do too");
        }

        [Test]
        public static void TheBestPeptidePepColumnReportsTheEstimateWhenThereIsOne()
        {
            var protein = new Protein("PEPTIDEK", "P1");
            var group = new ProteinGroup(new HashSet<Omics.IBioPolymer> { protein }, new HashSet<Omics.IBioPolymerWithSetMods>(), new HashSet<Omics.IBioPolymerWithSetMods>())
            {
                BestPeptidePEP = 0.01
            };
            Assert.That(group.ToString().Split('\t').Last(), Is.EqualTo(0.01.ToString()), "no estimate: the classifier's value");
            group.BestPeptideEstimatedPEP = 0.003;
            Assert.That(group.ToString().Split('\t').Last(), Is.EqualTo(0.003.ToString()));
        }
    }
}
