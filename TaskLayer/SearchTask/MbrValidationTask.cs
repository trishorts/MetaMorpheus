using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace TaskLayer.SearchTask
{
    class MbrValidationTask : MetaMorpheusTask
    {
        public readonly FlashLfqResults FlashLfqResults;
        public readonly SpectraFileInfo AcceptorFile;
        public readonly Ms2ScanWithSpecificMass[] MS2ScansInWindow;
        public readonly Dictionary<Proteomics.AminoAcidPolymer.Peptide, SpectralSimilarity> BestSimilarityForPeptides;


        public MbrValidationTask(FlashLfqResults flashLfqResults, SpectraFileInfo acceptorFile){

            AcceptorFile = acceptorFile;
            FlashLfqResults = flashLfqResults;

         }

        public override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList,
            List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyFileManager myFileManager = new MyFileManager(SearchParameters.DisposeOfFileWhenDone);
            MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams);
            var mbrPeaks = FlashLfqResults.Peaks[AcceptorFile].Where(s => s.IsMbrPeak);

            foreach (ChromatographicPeak peak in mbrPeaks)
            {
                double apexRT = peak.Apex.IndexedPeak.RetentionTime;
                double apexMz = peak.Apex.IndexedPeak.Mz;
                double peakHalfWidth = 1.0; //Placeholder value to determine retention time window

                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy(b => b.PrecursorMass).ToArray();

            }
        }


        // mark the file as in-progress
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                MassDiffAcceptor massDiffAcceptor = GetMassDiffAcceptor(combinedParams.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams);
                Status("Getting ms2 scans...", thisId);
                
                numMs2SpectraPerFile.Add(Path.GetFileNameWithoutExtension(origDataFile), new int[] { myMsDataFile.GetAllScansList().Count(p => p.MsnOrder == 2), arrayOfMs2ScansSortedByMass.Length });
                myFileManager.DoneWithFile(origDataFile);
    }
}
