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
            DecoyGlycositeResidues = new string[0];
            DecoyGlycositesAdjacentToRealSites = false;
            DecoyGlycositeCanonicalTarget = "T";
            MaxDecoyGlycositesPerPeptide = 0;
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

        /// <summary>
        /// Residues added as DECOY glycosites -- e.g. ["A"]. A localization that lands on one is wrong
        /// by construction, which is the ground truth a false-localization rate needs.
        ///
        /// EMPTY BY DEFAULT: an ordinary search is completely unaffected.
        ///
        /// This is construction (b) of at least five candidates (design/PLAN_GLYCO_LOCALIZATION.md
        /// section 1.0) and it is a PROBE, not a settled choice. The decoy-amino-acid family has known
        /// unexplained behaviour -- see the localization project's M14/M15.
        /// </summary>
        public string[] DecoyGlycositeResidues { get; set; }

        /// <summary>
        /// Construction (f): place decoy sites at positions IMMEDIATELY FLANKING a real candidate site,
        /// whatever residue is there (S/T excluded, so the site stays wrong by construction).
        ///
        /// Why position rather than residue: M13 measured that residue-chosen decoys never reach
        /// p >= 0.90 in any stratum, because a non-glycosylatable residue rarely sits where the
        /// fragment ions fail to exclude it. Adjacency is the regime where the evidence genuinely
        /// cannot separate two sites, and it is the only regime that can populate the high-confidence
        /// bins a Phred-like score depends on.
        ///
        /// NOTE the rate this produces is CONSERVATIVE, not unbiased: it deliberately samples the
        /// hardest subset. Never quote it beside a residue-decoy rate as the same quantity.
        /// </summary>
        public bool DecoyGlycositesAdjacentToRealSites { get; set; }

        /// <summary>
        /// Which single glycan-target letter a decoy site may draw from. "T" or "S". The restriction to
        /// one letter preserves competition parity; WHICH letter is arbitrary and is exposed so it can
        /// be flipped as a control -- see design/MEASUREMENTS.md M16.
        /// </summary>
        public string DecoyGlycositeCanonicalTarget { get; set; }

        /// <summary>
        /// Cap on decoy sites added per peptide. 0 = uncapped (every flanking position).
        ///
        /// WHY THIS IS A VALIDITY REQUIREMENT, not tidiness. The FLR estimator scales a decoy's
        /// probability mass by that peptide's target:decoy site ratio, which assumes decoy and
        /// wrong-target sites are EXCHANGEABLE. Uncapped adjacency put decoy sites at ~50% of all
        /// candidate sites (M14: 55%, M16: 42%, first construction-(a) run: 48.9%), and at that
        /// density the estimator overflows -- it reported 234 wrong sites among 121 target sites and
        /// clipped the answer. Half the localization graph being decoy is not a perturbation of the
        /// graph, it is a rewrite of it.
        ///
        /// Capping keeps the adjacency that makes decoys informative while leaving the graph close to
        /// what a real search builds.
        /// </summary>
        public int MaxDecoyGlycositesPerPeptide { get; set; }
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