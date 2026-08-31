using EngineLayer.GlycoSearch;
using System.Collections.Generic;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class GlycoSearchParameters
    {
        public GlycoSearchParameters()
        {
            OGlycanDatabasefile = "OGlycan.gdb";
            NGlycanDatabasefile = "NGlycan.gdb";
            GlycoSearchType = GlycoSearchType.OGlycanSearch;
            OxoniumIonFilt = true;
            DecoyType = DecoyType.Reverse;
            GlycoSearchTopNum = 50;
            RetainedGsmsPerScan = 25;
            MaximumOGlycanAllowed = 4;
            DoParsimony = true;
            NoOneHitWonders = false;
            ModPeptidesAreDifferent = false;

            //quantification options
            DoQuantification = false;
            DoMbrAnalysis = true;
            QuantifyPpmTol = 5;
            Normalize = false;

            //output options
            WriteIndividualFiles = false;
            WriteDecoys = true;
            WriteContaminants = true;
            WriteSpectrumLibrary = false;
            DisposeOfFileWhenDone = true;
            WritePrunedDataBase = false;

            ModsToWriteSelection = SearchParameters.DefaultModsToWriteSelection();
        }
        public string OGlycanDatabasefile { get; set; }
        public string NGlycanDatabasefile { get; set; }
        public GlycoSearchType GlycoSearchType { get; set; }
        public bool OxoniumIonFilt { get; set; }
        public DecoyType DecoyType { get; set; }
        public int GlycoSearchTopNum { get; set; }

        /// <summary>
        /// How many GlycoSpectralMatch objects to retain per MS2 scan. Was a hardcoded 10.
        ///
        /// 25, not 10, because in a combined target+decoy-glycosite search the two interleave roughly
        /// half and half, so the 10th TARGET match sits around position 20 and a retention of 10 yields
        /// about 5 targets. 25 gives ~10 targets plus headroom. Deeper buys little: the competition
        /// signal is concentrated in the top handful.
        /// </summary>
        public int RetainedGsmsPerScan { get; set; }
        public int MaximumOGlycanAllowed { get; set; }

        public bool DoParsimony { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool ModPeptidesAreDifferent { get; set; }
        
        //quantification options
        public bool DoQuantification { get; set; }
        public bool DoMbrAnalysis { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool Normalize { get; set; }

        //output options
        public bool WriteIndividualFiles { get; set; }
        public bool WriteDecoys { get; set; }
        public bool WriteContaminants { get; set; }
        public bool WriteSpectrumLibrary { get; set; }
        public bool WritePrunedDataBase { get; set; }
        public bool DisposeOfFileWhenDone { get; set; }

        public Dictionary<string, int> ModsToWriteSelection { get; set; }
    }
}