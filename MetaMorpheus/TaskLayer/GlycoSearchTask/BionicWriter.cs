using EngineLayer;
using EngineLayer.GlycoSearch;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace TaskLayer.GlycoSearchTask
{
    /// <summary>
    /// Writes glycopeptide results in Bionic format for downstream analysis.
    /// Format includes: Protein name, UniProtKB accession, peptide position, glycosylation sites,
    /// base sequence, modified sequence with dalton shifts, and peptide-relative glyco sites.
    /// </summary>
    public static class BionicWriter
    {
        public const string Header = "Protein name\tUniProtKB accession #\tIdentified peptide sequence starting position (UniProtKB numbering)\tGlycosylation site(s) (numbering with respect to UniprotKB numbering)\tIdentified peptide base sequence\tIdentified peptide with modifications (as dalton shifts)\tGlycosylation site(s) (numbering with respect to peptide sequence)";

        /// <summary>
        /// Writes glycopeptide results to a tab-separated file in Bionic format.
        /// </summary>
        /// <param name="gsms">List of GlycoSpectralMatches to write</param>
        /// <param name="filePath">Output file path</param>
        public static void WriteGlycoPeptidesToTsv(List<GlycoSpectralMatch> gsms, string filePath)
        {
            if (gsms == null || gsms.Count == 0)
            {
                return;
            }

            // Filter to only glycopeptides (those with glycan modifications)
            var glycoPeptides = gsms.Where(g => g.Routes != null || g.GlycanScore > 0).ToList();

            if (!glycoPeptides.Any())
            {
                return;
            }

            using StreamWriter output = new(filePath);
            output.WriteLine(Header);

            foreach (var gsm in glycoPeptides)
            {
                var line = GetBionicLine(gsm);
                if (!string.IsNullOrEmpty(line))
                {
                    output.WriteLine(line);
                }
            }
        }

        /// <summary>
        /// Generates a single Bionic-format line for a GlycoSpectralMatch.
        /// </summary>
        private static string GetBionicLine(GlycoSpectralMatch gsm)
        {
            var sb = new StringBuilder();

            // Get protein information
            var peptide = gsm.BestMatchingBioPolymersWithSetMods.FirstOrDefault()?.SpecificBioPolymer as PeptideWithSetModifications;
            if (peptide == null)
            {
                return null;
            }

            // 1. Protein name
            string proteinName = peptide.Protein?.FullName ?? "";
            sb.Append(proteinName);
            sb.Append('\t');

            // 2. UniProtKB accession #
            string accession = gsm.Accession ?? peptide.Protein?.Accession ?? "";
            sb.Append(accession);
            sb.Append('\t');

            // 3. Identified peptide sequence starting position (UniProtKB numbering)
            int startPosition = gsm.OneBasedStartResidue ?? peptide.OneBasedStartResidue;
            sb.Append(startPosition);
            sb.Append('\t');

            // 4. Glycosylation site(s) (numbering with respect to UniProtKB numbering)
            var glycoSitesProtein = GetGlycosylationSitesProteinNumbering(gsm, startPosition);
            sb.Append(string.Join(";", glycoSitesProtein.Select(s => s.ToString(CultureInfo.InvariantCulture))));
            sb.Append('\t');

            // 5. Identified peptide base sequence
            sb.Append(gsm.BaseSequence);
            sb.Append('\t');

            // 6. Identified peptide with modifications (as dalton shifts)
            string modifiedSequence = GetSequenceWithDaltonShifts(gsm, peptide);
            sb.Append(modifiedSequence);
            sb.Append('\t');

            // 7. Glycosylation site(s) (numbering with respect to peptide sequence)
            var glycoSitesPeptide = GetGlycosylationSitesPeptideNumbering(gsm);
            sb.Append(string.Join(";", glycoSitesPeptide));

            return sb.ToString();
        }

        /// <summary>
        /// Gets glycosylation sites in protein numbering (UniProtKB).
        /// </summary>
        private static List<double> GetGlycosylationSitesProteinNumbering(GlycoSpectralMatch gsm, int peptideStartPosition)
        {
            var sites = new List<double>();

            if (gsm.Routes != null && gsm.LocalizedGlycan != null)
            {
                // O-glycan or combined search with localization
                foreach (var localizedMod in gsm.LocalizedGlycan.Where(lg => lg.Confident))
                {
                    // SiteIndex is 1-based within peptide, convert to protein position
                    int proteinSite = peptideStartPosition + localizedMod.SiteIndex - 2;
                    sites.Add(proteinSite);
                }

                // If no confident localizations, use all possible sites
                if (!sites.Any() && gsm.Routes.Any())
                {
                    foreach (var route in gsm.Routes.Take(1)) // Use first route as representative
                    {
                        foreach (var modSitePair in route.ModSitePairs)
                        {
                            int proteinSite = peptideStartPosition + modSitePair.SiteIndex - 2;
                            if (!sites.Contains(proteinSite))
                            {
                                sites.Add(proteinSite);
                            }
                        }
                    }
                }
            }
            else if (gsm.NGlycanLocalizations != null && gsm.NGlycanLocalizations.Any())
            {
                // N-glycan search
                foreach (var site in gsm.NGlycanLocalizations)
                {
                    int proteinSite = peptideStartPosition + site - 1;
                    sites.Add(proteinSite);
                }
            }
            else if (gsm.GlycanScore > 0)
            {
                // N-glycan without explicit localization - try to find N-X-S/T motif
                var nGlycoSites = FindNGlycanMotifSites(gsm.BaseSequence);
                foreach (var site in nGlycoSites)
                {
                    int proteinSite = peptideStartPosition + site - 1;
                    sites.Add(proteinSite);
                }
            }

            return sites.OrderBy(s => s).ToList();
        }

        /// <summary>
        /// Gets glycosylation sites in peptide numbering (1-based).
        /// </summary>
        private static List<int> GetGlycosylationSitesPeptideNumbering(GlycoSpectralMatch gsm)
        {
            var sites = new List<int>();

            if (gsm.Routes != null && gsm.LocalizedGlycan != null)
            {
                // O-glycan or combined search with localization
                foreach (var localizedMod in gsm.LocalizedGlycan.Where(lg => lg.Confident))
                {
                    // SiteIndex is 1-based, but output should be 0-based per example (subtract 1)
                    int peptideSite = localizedMod.SiteIndex - 1;
                    sites.Add(peptideSite);
                }

                // If no confident localizations, use all possible sites from first route
                if (!sites.Any() && gsm.Routes.Any())
                {
                    foreach (var modSitePair in gsm.Routes.First().ModSitePairs)
                    {
                        int peptideSite = modSitePair.SiteIndex - 1;
                        if (!sites.Contains(peptideSite))
                        {
                            sites.Add(peptideSite);
                        }
                    }
                }
            }
            else if (gsm.NGlycanLocalizations != null && gsm.NGlycanLocalizations.Any())
            {
                // N-glycan search - these are already peptide-relative
                sites.AddRange(gsm.NGlycanLocalizations);
            }
            else if (gsm.GlycanScore > 0)
            {
                // N-glycan without explicit localization - find N-X-S/T motif
                sites.AddRange(FindNGlycanMotifSites(gsm.BaseSequence));
            }

            return sites.OrderBy(s => s).ToList();
        }

        /// <summary>
        /// Finds N-glycan sequon (N-X-S/T where X != P) positions in peptide sequence.
        /// Returns 1-based positions of the asparagine.
        /// </summary>
        private static List<int> FindNGlycanMotifSites(string baseSequence)
        {
            var sites = new List<int>();

            for (int i = 0; i < baseSequence.Length - 2; i++)
            {
                if (baseSequence[i] == 'N' &&
                    baseSequence[i + 1] != 'P' &&
                    (baseSequence[i + 2] == 'S' || baseSequence[i + 2] == 'T'))
                {
                    sites.Add(i + 1); // 1-based position
                }
            }

            return sites;
        }

        /// <summary>
        /// Generates the peptide sequence with modifications shown as dalton shifts.
        /// Example: GLYCN[1216.4228]ASM[15.9949]PEPTIDEK
        /// </summary>
        private static string GetSequenceWithDaltonShifts(GlycoSpectralMatch gsm, PeptideWithSetModifications peptide)
        {
            var sb = new StringBuilder();
            string baseSequence = gsm.BaseSequence;

            // Build a dictionary of position -> total modification mass
            var positionToMass = new Dictionary<int, double>();

            // Add standard modifications from the peptide
            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                int position = mod.Key - 1; // Convert from 1-based to 0-based
                if (position >= 1 && position <= baseSequence.Length)
                {
                    int seqIndex = position - 1; // 0-based index in base sequence
                    if (!positionToMass.ContainsKey(seqIndex))
                    {
                        positionToMass[seqIndex] = 0;
                    }
                    positionToMass[seqIndex] += mod.Value.MonoisotopicMass ?? 0;
                }
            }

            // Add glycan modifications
            if (gsm.Routes != null && gsm.Routes.Any())
            {
                // O-glycan or combined search
                var route = gsm.Routes.First();
                foreach (var modSitePair in route.ModSitePairs)
                {
                    int seqIndex = modSitePair.SiteIndex - 2; // Convert to 0-based
                    if (seqIndex >= 0 && seqIndex < baseSequence.Length)
                    {
                        double glycanMass = modSitePair.ModId >= 0
                            ? GlycanBox.GlobalOGlycans[modSitePair.ModId].Mass / 1E5
                            : GlycanBox.GlobalNGlycans[modSitePair.ModId].Mass / 1E5;

                        if (!positionToMass.ContainsKey(seqIndex))
                        {
                            positionToMass[seqIndex] = 0;
                        }
                        positionToMass[seqIndex] += glycanMass;
                    }
                }
            }
            else if (gsm.NGlycan != null && gsm.NGlycan.Any())
            {
                // N-glycan search
                var nGlycoSites = gsm.NGlycanLocalizations ?? FindNGlycanMotifSites(baseSequence);
                if (nGlycoSites.Any())
                {
                    int siteIndex = nGlycoSites.First() - 1; // 0-based
                    if (siteIndex >= 0 && siteIndex < baseSequence.Length)
                    {
                        double glycanMass = gsm.NGlycan.First().Mass / 1E5;
                        if (!positionToMass.ContainsKey(siteIndex))
                        {
                            positionToMass[siteIndex] = 0;
                        }
                        positionToMass[siteIndex] += glycanMass;
                    }
                }
            }

            // Build the sequence string with modifications
            for (int i = 0; i < baseSequence.Length; i++)
            {
                sb.Append(baseSequence[i]);

                if (positionToMass.TryGetValue(i, out double mass) && Math.Abs(mass) > 0.0001)
                {
                    sb.Append('[');
                    sb.Append(mass.ToString("F4", CultureInfo.InvariantCulture));
                    sb.Append(']');
                }
            }

            return sb.ToString();
        }

        /// <summary>
        /// Gets the tab-separated header for the Bionic format output file.
        /// </summary>
        public static string GetHeader()
        {
            return Header;
        }
    }
}
