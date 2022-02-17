﻿using EngineLayer;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer
{
    public class SingleCellMBRTask : MetaMorpheusTask
    {
        public SingleCellMBRTask() : base(MyTask.SingleCellProteomics)
        {
            CommonParameters = new CommonParameters();
            SearchParameters = new SearchParameters();
        }

        public SearchParameters SearchParameters { get; set; }

        public static MassDiffAcceptor GetMassDiffAcceptor(Tolerance precursorMassTolerance, MassDiffAcceptorType massDiffAcceptorType, string customMdac)
        {
            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact:
                    if (precursorMassTolerance is PpmTolerance)
                    {
                        return new SinglePpmAroundZeroSearchMode(precursorMassTolerance.Value);
                    }
                    else
                    {
                        return new SingleAbsoluteAroundZeroSearchMode(precursorMassTolerance.Value);
                    }

                case MassDiffAcceptorType.OneMM:
                    return new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, precursorMassTolerance);

                case MassDiffAcceptorType.TwoMM:
                    return new DotMassDiffAcceptor("2mm", new List<double> { 0, 1.0029, 2.0052 }, precursorMassTolerance);

                case MassDiffAcceptorType.ThreeMM:
                    return new DotMassDiffAcceptor("3mm", new List<double> { 0, 1.0029, 2.0052, 3.0077 }, precursorMassTolerance);

                case MassDiffAcceptorType.ModOpen:
                    return new IntervalMassDiffAcceptor("-187andUp", new List<DoubleRange> { new DoubleRange(-187, double.PositiveInfinity) });

                case MassDiffAcceptorType.Open:
                    return new OpenSearchMode();

                case MassDiffAcceptorType.Custom:
                    return ParseSearchMode(customMdac);

                case MassDiffAcceptorType.PlusOrMinusThreeMM:
                    return new DotMassDiffAcceptor(
                        "PlusOrMinus3Da",
                        new List<double>
                        {
                            -3 * Chemistry.Constants.C13MinusC12,
                            -2 * Chemistry.Constants.C13MinusC12,
                            -1 * Chemistry.Constants.C13MinusC12,
                            0,
                            1 * Chemistry.Constants.C13MinusC12,
                            2 * Chemistry.Constants.C13MinusC12,
                            3 * Chemistry.Constants.C13MinusC12
                        },
                        precursorMassTolerance);

                default:
                    throw new MetaMorpheusException("Unknown MassDiffAcceptorType");
            }
        }

        private static MassDiffAcceptor ParseSearchMode(string text)
        {
            MassDiffAcceptor massDiffAcceptor = null;

            try
            {
                var split = text.Split(' ');

                switch (split[1])
                {
                    case "dot":
                        double[] massShifts = split[4].Split(',').Select(p => double.Parse(p, CultureInfo.InvariantCulture)).ToArray();
                        string newString = split[2].Replace("�", "");
                        double toleranceValue = double.Parse(newString, CultureInfo.InvariantCulture);
                        if (split[3].ToUpperInvariant().Equals("PPM"))
                        {
                            massDiffAcceptor = new DotMassDiffAcceptor(split[0], massShifts, new PpmTolerance(toleranceValue));
                        }
                        else if (split[3].ToUpperInvariant().Equals("DA"))
                        {
                            massDiffAcceptor = new DotMassDiffAcceptor(split[0], massShifts, new AbsoluteTolerance(toleranceValue));
                        }

                        break;

                    case "interval":
                        IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(';'), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split(',')[0],
                            CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split(',')[1], CultureInfo.InvariantCulture)));
                        massDiffAcceptor = new IntervalMassDiffAcceptor(split[0], doubleRanges);
                        break;

                    case "OpenSearch":
                        massDiffAcceptor = new OpenSearchMode();
                        break;

                    case "daltonsAroundZero":
                        massDiffAcceptor = new SingleAbsoluteAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                        break;

                    case "ppmAroundZero":
                        massDiffAcceptor = new SinglePpmAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                        break;

                    default:
                        throw new MetaMorpheusException("Unrecognized search mode type: " + split[1]);
                }
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException("Could not parse search mode string: " + e.Message);
            }

            return massDiffAcceptor;
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList, List<PeptideSpectralMatch> allPeptides, FlashLfqResults flashLfqResults)
        {

            // write prose settings
            ProseCreatedWhileRunning.Append("Starting single-cell proteomics MBR analysis.");

            if (SearchParameters.DoQuantification)
            {
                // disable quantification if a .mgf is being used
                if (currentRawFileList.Any(x => Path.GetExtension(x).Equals(".mgf", StringComparison.OrdinalIgnoreCase)))
                {
                    SearchParameters.DoQuantification = false;
                }

            }

            LoadModifications(taskId, out var variableModifications, out var fixedModifications, out var localizeableModificationTypes);
            ProseCreatedWhileRunning.Append("Modifications loaded.");

            // load proteins
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, SearchParameters.SearchTarget, SearchParameters.DecoyType, localizeableModificationTypes, CommonParameters);
            SanitizeProteinDatabase(proteinList, SearchParameters.TCAmbiguity);
            ProseCreatedWhileRunning.Append("Protein list loaded.");


            // load spectral libraries
            var spectralLibrary = LoadSpectralLibraries(taskId, dbFilenameList);
            ProseCreatedWhileRunning.Append("Spectral library loaded.");

            // start the single cell MBR task
            MyTaskResults = new MyTaskResults(this);


            MyFileManager myFileManager = new MyFileManager(SearchParameters.DisposeOfFileWhenDone);

            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();

            Status("Searching files...", taskId);
            Status("Searching files...", new List<string> { taskId, "Individual Spectra Files" });

            Dictionary<string, int[]> numMs2SpectraPerFile = new Dictionary<string, int[]>();
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                if (GlobalVariables.StopLoops) { break; }

                var origDataFile = currentRawFileList[spectraFileIndex];

                // mark the file as in-progress
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                MassDiffAcceptor massDiffAcceptor = GetMassDiffAcceptor(combinedParams.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams);
                Status("Getting ms2 scans...", thisId);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy(b => b.PrecursorMass).ToArray();
                numMs2SpectraPerFile.Add(Path.GetFileNameWithoutExtension(origDataFile), new int[] { myMsDataFile.GetAllScansList().Count(p => p.MsnOrder == 2), arrayOfMs2ScansSortedByMass.Length });
                myFileManager.DoneWithFile(origDataFile);

                PeptideSpectralMatch[] fileSpecificPsms = new PeptideSpectralMatch[arrayOfMs2ScansSortedByMass.Length];

                // modern search
                if (SearchParameters.SearchType == SearchType.Modern)
                {
                    for (int currentPartition = 0; currentPartition < combinedParams.TotalPartitions; currentPartition++)
                    {
                        List<PeptideWithSetModifications> peptideIndex = null;
                        List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count / combinedParams.TotalPartitions,
                            ((currentPartition + 1) * proteinList.Count / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count / combinedParams.TotalPartitions));

                        Status("Getting fragment dictionary...", new List<string> { taskId });
                        var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, SearchParameters.SilacLabels,
                            SearchParameters.StartTurnoverLabel, SearchParameters.EndTurnoverLabel, currentPartition, SearchParameters.DecoyType, combinedParams, FileSpecificParameters,
                            SearchParameters.MaxFragmentSize, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), SearchParameters.TCAmbiguity, new List<string> { taskId });
                        List<int>[] fragmentIndex = null;
                        List<int>[] precursorIndex = null;

                        lock (indexLock)
                        {
                            GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);
                        }

                        Status("Searching files...", taskId);

                        new ModernSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, currentPartition,
                            combinedParams, this.FileSpecificParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                        ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + combinedParams.TotalPartitions + "!", thisId));
                        if (GlobalVariables.StopLoops) { break; }
                    }
                }
                // nonspecific search
                else if (SearchParameters.SearchType == SearchType.NonSpecific)
                {
                    PeptideSpectralMatch[][] fileSpecificPsmsSeparatedByFdrCategory = new PeptideSpectralMatch[numFdrCategories][]; //generate an array of all possible locals
                    for (int i = 0; i < numFdrCategories; i++) //only add if we're using for FDR, else ignore it as null.
                    {
                        fileSpecificPsmsSeparatedByFdrCategory[i] = new PeptideSpectralMatch[arrayOfMs2ScansSortedByMass.Length];
                    }

                    //create params for N, C, or both if semi
                    List<CommonParameters> paramsToUse = new List<CommonParameters> { combinedParams };
                    if (combinedParams.DigestionParams.SearchModeType == CleavageSpecificity.Semi) //if semi, we need to do both N and C to hit everything
                    {
                        paramsToUse.Clear();
                        List<FragmentationTerminus> terminiToUse = new List<FragmentationTerminus> { FragmentationTerminus.N, FragmentationTerminus.C };
                        foreach (FragmentationTerminus terminus in terminiToUse) //set both termini
                        {
                            paramsToUse.Add(combinedParams.CloneWithNewTerminus(terminus));
                        }
                    }

                    //Compress array of deconvoluted ms2 scans to avoid searching the same ms2 multiple times while still identifying coisolated peptides
                    List<int>[] coisolationIndex = new List<int>[] { new List<int>() };
                    if (arrayOfMs2ScansSortedByMass.Length != 0)
                    {
                        int maxScanNumber = arrayOfMs2ScansSortedByMass.Max(x => x.OneBasedScanNumber);
                        coisolationIndex = new List<int>[maxScanNumber + 1];
                        for (int i = 0; i < arrayOfMs2ScansSortedByMass.Length; i++)
                        {
                            int scanNumber = arrayOfMs2ScansSortedByMass[i].OneBasedScanNumber;
                            if (coisolationIndex[scanNumber] == null)
                            {
                                coisolationIndex[scanNumber] = new List<int> { i };
                            }
                            else
                            {
                                coisolationIndex[scanNumber].Add(i);
                            }
                        }
                        coisolationIndex = coisolationIndex.Where(x => x != null).ToArray();
                    }

                    //foreach terminus we're going to look at
                    foreach (CommonParameters paramToUse in paramsToUse)
                    {
                        //foreach database partition
                        for (int currentPartition = 0; currentPartition < paramToUse.TotalPartitions; currentPartition++)
                        {
                            List<PeptideWithSetModifications> peptideIndex = null;

                            List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count / paramToUse.TotalPartitions,
                                ((currentPartition + 1) * proteinList.Count / paramToUse.TotalPartitions) - (currentPartition * proteinList.Count / paramToUse.TotalPartitions));

                            List<int>[] fragmentIndex = null;
                            List<int>[] precursorIndex = null;

                            Status("Getting fragment dictionary...", new List<string> { taskId });
                            var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, SearchParameters.SilacLabels,
                                SearchParameters.StartTurnoverLabel, SearchParameters.EndTurnoverLabel, currentPartition, SearchParameters.DecoyType, paramToUse, FileSpecificParameters,
                                SearchParameters.MaxFragmentSize, true, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), SearchParameters.TCAmbiguity, new List<string> { taskId });
                            lock (indexLock)
                            {
                                GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);
                            }

                            Status("Searching files...", taskId);

                            new NonSpecificEnzymeSearchEngine(fileSpecificPsmsSeparatedByFdrCategory, arrayOfMs2ScansSortedByMass, coisolationIndex, peptideIndex, fragmentIndex,
                                precursorIndex, currentPartition, paramToUse, this.FileSpecificParameters, variableModifications, massDiffAcceptor,
                                SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                            ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + paramToUse.TotalPartitions + "!", thisId));
                            if (GlobalVariables.StopLoops) { break; }
                        }
                    }
                    lock (psmLock)
                    {
                        for (int i = 0; i < allCategorySpecificPsms.Length; i++)
                        {
                            if (allCategorySpecificPsms[i] != null)
                            {
                                allCategorySpecificPsms[i].AddRange(fileSpecificPsmsSeparatedByFdrCategory[i]);
                            }
                        }
                    }
                }
                // classic search
                else
                {
                    Status("Starting search...", thisId);
                    var newClassicSearchEngine = new ClassicSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, SearchParameters.SilacLabels,
                       SearchParameters.StartTurnoverLabel, SearchParameters.EndTurnoverLabel, proteinList, massDiffAcceptor, combinedParams, this.FileSpecificParameters, spectralLibrary, thisId, SearchParameters.WriteSpectralLibrary);
                    newClassicSearchEngine.Run();

                    ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));
                }

                //look for internal fragments
                if (SearchParameters.MinAllowedInternalFragmentLength != 0)
                {
                    MatchInternalFragmentIons(fileSpecificPsms, arrayOfMs2ScansSortedByMass, combinedParams, SearchParameters.MinAllowedInternalFragmentLength);
                }

                // calculate/set spectral angles if there is a spectral library being used
                if (spectralLibrary != null)
                {
                    Status("Calculating spectral library similarity...", thisId);
                }
                SpectralLibrarySearchFunction.CalculateSpectralAngles(spectralLibrary, fileSpecificPsms, arrayOfMs2ScansSortedByMass, combinedParams);

                lock (psmLock)
                {
                    allPeptides.AddRange(fileSpecificPsms);
                }

                completedFiles++;
                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            if (spectralLibrary != null)
            {
                spectralLibrary.CloseConnections();
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            int numNotches = GetNumNotches(SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            //resolve category specific fdrs (for speedy semi and nonspecific
            if (SearchParameters.SearchType == SearchType.NonSpecific)
            {
                allPeptides = NonSpecificEnzymeSearchEngine.ResolveFdrCategorySpecificPsms(allCategorySpecificPsms, numNotches, taskId, CommonParameters, FileSpecificParameters);
            }

            PostSearchAnalysisParameters parameters = new PostSearchAnalysisParameters
            {
                SearchTaskResults = MyTaskResults,
                SearchTaskId = taskId,
                SearchParameters = SearchParameters,
                ProteinList = proteinList,
                AllPsms = allPeptides,
                VariableModifications = variableModifications,
                FixedModifications = fixedModifications,
                ListOfDigestionParams = new HashSet<DigestionParams>(fileSpecificCommonParams.Select(p => p.DigestionParams)),
                CurrentRawFileList = currentRawFileList,
                MyFileManager = myFileManager,
                NumNotches = numNotches,
                OutputFolder = OutputFolder,
                IndividualResultsOutputFolder = Path.Combine(OutputFolder, "Individual File Results"),
                FlashLfqResults = flashLfqResults,
                FileSettingsList = fileSettingsList,
                NumMs2SpectraPerFile = numMs2SpectraPerFile,
                DatabaseFilenameList = dbFilenameList
            };
            PostSearchAnalysisTask postProcessing = new PostSearchAnalysisTask
            {
                Parameters = parameters,
                FileSpecificParameters = this.FileSpecificParameters,
                CommonParameters = CommonParameters
            };
            return postProcessing.Run();
        }
    }
}