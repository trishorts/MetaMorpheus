using System.Collections.Generic;
using System.Text;
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
    /// <see cref="Protein"/> is built per accession so protein parsimony/grouping is correct. The .msl does
    /// not store the protein sequence, so each shared protein's sequence is the CONCATENATION of its
    /// peptides' base sequences and each peptide is given its slice's start/end. With 0-missed-cleavage
    /// tryptic peptides (non-overlapping) this concatenation approximates the original protein, which keeps
    /// sequence-coverage well-defined and in-bounds. PSM/peptide IDs, masses, and scores are exact because a
    /// string-constructed peptide derives its base sequence and mass from its own full sequence, not from
    /// the parent protein.
    /// </summary>
    public static class MslPeptideReader
    {
        // One entry's reconstructed data, pending creation of its shared parent protein.
        private readonly struct PendingPeptide
        {
            public readonly string FullSequence;
            public readonly int Start;          // 1-based start of this peptide in the concatenated protein
            public readonly int Length;         // base-sequence length
            public readonly List<Product> Fragments;

            public PendingPeptide(string fullSequence, int start, int length, List<Product> fragments)
            {
                FullSequence = fullSequence;
                Start = start;
                Length = length;
                Fragments = fragments;
            }
        }

        private sealed class AccessionGroup
        {
            public readonly StringBuilder Sequence = new();
            public readonly List<PendingPeptide> Peptides = new();
            public bool IsDecoy;
        }

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

            // The .msl peptides were digested with trypsin / 0 missed cleavages. A non-null
            // DigestionParams is required downstream (e.g. parsimony keys on DigestionAgent).
            var digestionParams = new DigestionParams("trypsin", maxMissedCleavages: 0);

            // Pass 1: read entries, grouping peptides by accession and laying them out end-to-end so each
            // accession's shared protein gets one concatenated sequence with valid per-peptide offsets.
            var groups = new Dictionary<string, AccessionGroup>();
            string fallbackAccession = $"{databaseName}_UNKNOWN";

            using (var library = MslLibrary.Load(mslPath))
            {
                foreach (var entry in library.GetAllEntries(includeDecoys: true))
                {
                    string fullSequence = entry.FullSequence;
                    string baseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSequence);

                    string accession = string.IsNullOrEmpty(entry.ProteinAccession)
                        ? fallbackAccession
                        : entry.ProteinAccession;

                    if (!groups.TryGetValue(accession, out var group))
                    {
                        group = new AccessionGroup { IsDecoy = entry.IsDecoy };
                        groups[accession] = group;
                    }

                    int start = group.Sequence.Length + 1; // 1-based
                    group.Sequence.Append(baseSequence);
                    group.Peptides.Add(new PendingPeptide(
                        fullSequence, start, baseSequence.Length, BuildProducts(entry.MatchedFragmentIons)));
                }
            }

            // Pass 2: materialize one shared Protein per accession, then the peptides over their slices.
            var result = new List<(IBioPolymerWithSetMods, List<Product>)>();
            foreach (var kvp in groups)
            {
                var protein = new Protein(kvp.Value.Sequence.ToString(), kvp.Key, isDecoy: kvp.Value.IsDecoy);

                foreach (var pending in kvp.Value.Peptides)
                {
                    var peptide = new PeptideWithSetModifications(
                        pending.FullSequence, mods, p: protein, digestionParams: digestionParams,
                        oneBasedStartResidueInProtein: pending.Start,
                        oneBasedEndResidueInProtein: pending.Start + pending.Length - 1,
                        missedCleavages: 0);

                    result.Add((peptide, pending.Fragments));
                }
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
