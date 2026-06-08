using System.Collections.Generic;
using Chemistry;
using Omics;
using Omics.Fragmentation;
using Omics.Fragmentation.Peptide;
using Omics.SpectralMatch.MslSpectralLibrary;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;

namespace EngineLayer.ParallelSearch
{
    /// <summary>
    /// Reads a MetaMorpheus <c>.msl</c> spectral library as a PEPTIDE SOURCE for the parallel search
    /// (Strategy A): each library entry's full sequence is reconstructed into a searchable
    /// <see cref="PeptideWithSetModifications"/>, and its STORED fragment ions are returned alongside it
    /// as a ready <see cref="Product"/> list. The search matches these stored fragments directly instead
    /// of re-fragmenting the peptide — this is the actual speedup (digestion AND fragmentation are skipped).
    ///
    /// The .msl stores fragment m/z as float32 (~0.5 ppm vs a double recomputation). Because the product
    /// mass tolerance (±20–30 ppm) is far larger than that error, the matched-peak set is unchanged, so
    /// IDs and scores are preserved while fragmentation cost is eliminated.
    ///
    /// Accession handling: entries carry their source-protein accession (DECOY_… for decoys). One shared
    /// <see cref="Protein"/> is built per accession so protein parsimony/grouping is correct. The shared
    /// protein's "sequence" is just one of its peptides, so per-peptide residue start/end positions
    /// (used for protein-coverage display) are approximate; PSM/peptide IDs, masses, and scores are exact
    /// because a string-constructed peptide derives its base sequence and mass from its own full sequence,
    /// not from the parent protein.
    /// </summary>
    public static class MslPeptideReader
    {
        /// <summary>
        /// Reconstructs all peptides from a .msl library together with their stored (float) fragments.
        /// Mods are resolved against <see cref="GlobalVariables.AllModsKnownDictionary"/> so the
        /// reconstructed peptides match the search exactly.
        /// </summary>
        /// <param name="mslPath">Path to the .msl library.</param>
        /// <param name="databaseName">Namespaces the fallback accession for entries that lack one.</param>
        public static List<(IBioPolymerWithSetMods Peptide, List<Product> Fragments)> ReadPeptides(
            string mslPath, string databaseName)
        {
            var mods = GlobalVariables.AllModsKnownDictionary;
            var result = new List<(IBioPolymerWithSetMods, List<Product>)>();

            // The .msl peptides were digested with trypsin / 0 missed cleavages. A non-null
            // DigestionParams is required downstream (e.g. parsimony keys on DigestionAgent).
            var digestionParams = new DigestionParams("trypsin", maxMissedCleavages: 0);

            // One shared Protein per accession → correct parsimony / protein grouping.
            var proteinsByAccession = new Dictionary<string, Protein>();

            using var library = MslLibrary.Load(mslPath);
            foreach (var entry in library.GetAllEntries(includeDecoys: true))
            {
                string fullSequence = entry.FullSequence;
                string baseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSequence);

                string accession = string.IsNullOrEmpty(entry.ProteinAccession)
                    ? $"{databaseName}_UNKNOWN"
                    : entry.ProteinAccession;

                if (!proteinsByAccession.TryGetValue(accession, out var protein))
                {
                    protein = new Protein(baseSequence, accession, isDecoy: entry.IsDecoy);
                    proteinsByAccession[accession] = protein;
                }

                var peptide = new PeptideWithSetModifications(
                    fullSequence, mods, p: protein, digestionParams: digestionParams,
                    oneBasedStartResidueInProtein: 1,
                    oneBasedEndResidueInProtein: baseSequence.Length,
                    missedCleavages: 0);

                result.Add((peptide, BuildProducts(entry.MatchedFragmentIons)));
            }

            return result;
        }

        /// <summary>
        /// Turns the library's stored fragment ions into <see cref="Product"/> objects the search can
        /// match directly. The stored float m/z is converted to neutral mass via
        /// <see cref="ClassExtensions.ToMass(double, int)"/> — the on-disk reconstruction leaves
        /// <see cref="Product.NeutralMass"/> at 0, so it MUST be set here for the scorer to match.
        /// </summary>
        private static List<Product> BuildProducts(List<MslFragmentIon> ions)
        {
            var products = new List<Product>(ions.Count);
            foreach (var ion in ions)
            {
                FragmentationTerminus terminus =
                    TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus.TryGetValue(
                        ion.ProductType, out var t)
                        ? t
                        : FragmentationTerminus.None;

                products.Add(new Product(
                    ion.ProductType, terminus, ion.Mz.ToMass(ion.Charge), ion.FragmentNumber,
                    ion.ResiduePosition, ion.NeutralLoss, ion.SecondaryProductType,
                    ion.SecondaryFragmentNumber));
            }
            return products;
        }
    }
}
