using System.Collections.Generic;
using Omics;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;

namespace EngineLayer.ParallelSearch
{
    /// <summary>
    /// Reads a MetaMorpheus <c>.msl</c> spectral library as a PEPTIDE SOURCE for the parallel search
    /// (Strategy A): each library entry's full sequence is reconstructed into a searchable
    /// <see cref="PeptideWithSetModifications"/>. The search then re-fragments these in double
    /// precision (so the stored float32 library fragments are not used and byte-identity is preserved);
    /// the library simply supplies the precomputed (e.g. 0-missed-cleavage) target+decoy peptide set,
    /// removing per-database digestion.
    ///
    /// Accession caveat: the current .msl does not store the source-protein accession, so each peptide
    /// is given its own synthetic <see cref="Protein"/> parent (decoy flag preserved). PSM- and
    /// peptide-level results are unaffected, but protein parsimony/grouping will treat each peptide as
    /// its own protein. To restore real protein grouping, store the accession in the .msl (the
    /// MslLibraryEntry model has the field) and key the synthetic proteins by it.
    /// </summary>
    public static class MslPeptideReader
    {
        /// <summary>
        /// Reconstructs all peptides from a .msl library. Mods are resolved against
        /// <see cref="GlobalVariables.AllModsKnownDictionary"/> so they match the search exactly.
        /// </summary>
        /// <param name="mslPath">Path to the .msl library.</param>
        /// <param name="databaseName">Used to namespace the synthetic protein accessions.</param>
        public static List<IBioPolymerWithSetMods> ReadPeptides(string mslPath, string databaseName)
        {
            var mods = GlobalVariables.AllModsKnownDictionary;
            var peptides = new List<IBioPolymerWithSetMods>();

            // The .msl peptides were digested with trypsin / 0 missed cleavages. A non-null
            // DigestionParams is required downstream (e.g. parsimony keys on DigestionAgent).
            var digestionParams = new DigestionParams("trypsin", maxMissedCleavages: 0);

            using var library = new SpectralLibrary(new List<string> { mslPath });

            int index = 0;
            foreach (var spectrum in library.GetAllLibrarySpectra())
            {
                string fullSequence = spectrum.Sequence;
                string baseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSequence);

                // One synthetic protein per peptide; decoy flag carried through for FDR.
                var protein = new Protein(baseSequence, $"{databaseName}_{index++}", isDecoy: spectrum.IsDecoy);

                var peptide = new PeptideWithSetModifications(
                    fullSequence, mods, p: protein, digestionParams: digestionParams,
                    oneBasedStartResidueInProtein: 1,
                    oneBasedEndResidueInProtein: baseSequence.Length,
                    missedCleavages: 0);

                peptides.Add(peptide);
            }

            return peptides;
        }
    }
}
