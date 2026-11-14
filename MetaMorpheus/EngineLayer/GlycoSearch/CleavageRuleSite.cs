using Omics.Digestion;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.GlycoSearch
{
    /// <summary>
    /// A protein residue where a cut the protease made says a glycan must be (OpeRATOR at P1', StcE at P2), or must not be
    /// (IMPa at P1), with what the identifications say about it. The rule is checked against the spectra here, never
    /// enforced: occupancy varies between molecules, and a peptide end can be made by something other than the protease.
    /// </summary>
    public class CleavageRuleSite
    {
        public string ProteinAccession { get; init; }
        public int OneBasedResidue { get; init; }
        public char AminoAcid { get; init; }
        public bool IsForbidden { get; init; }

        /// <summary>PSMs one of whose ends, if the protease cut it, implies this.</summary>
        public int ImpliedBy { get; set; }

        /// <summary>Of those, Level 1 PSMs that place the residue as the rule says, and that place it the other way. Only a PSM
        /// containing the residue counts; a cut's residue on the neighbouring peptide is not in it.</summary>
        public int ImpliersAgreeing { get; set; }
        public int ImpliersAgainst { get; set; }

        /// <summary>Other Level 1 PSMs containing the residue: with a glycan on it, and without.</summary>
        public int OthersWithGlycan { get; set; }
        public int OthersWithout { get; set; }

        /// <param name="psms">For each confident target PSM: its peptide, the agent that digested it, and, if it is Level 1, the
        /// two-based site keys its route places glycans on (<see cref="Route.ModSitePairs"/>; null otherwise).</param>
        public static List<CleavageRuleSite> Find(IEnumerable<(DigestionProduct Peptide, DigestionAgent Agent, IEnumerable<int> Level1SiteKeys)> psms)
        {
            var list = psms.Select(p => (p.Peptide, p.Agent,
                Glycosylated: p.Level1SiteKeys?.Select(key => p.Peptide.OneBasedStartResidue + key - 2).ToHashSet())).ToList();
            var sites = new Dictionary<(string Accession, int Residue, bool IsForbidden), (CleavageRuleSite Site, HashSet<int> Impliers)>();
            for (int i = 0; i < list.Count; i++)
            {
                var (peptide, agent, _) = list[i];
                foreach (var (residue, requirement) in peptide.GetCleavageConditions(agent))
                {
                    var key = (peptide.Parent.Accession, residue, requirement.IsForbidden);
                    if (!sites.TryGetValue(key, out var entry))
                    {
                        entry = (new CleavageRuleSite { ProteinAccession = key.Accession, OneBasedResidue = residue,
                            AminoAcid = peptide.Parent.BaseSequence[residue - 1], IsForbidden = requirement.IsForbidden }, new HashSet<int>());
                        sites[key] = entry;
                    }
                    if (entry.Impliers.Add(i))
                    {
                        entry.Site.ImpliedBy++;
                    }
                }
            }

            foreach (var (key, (site, impliers)) in sites)
            {
                for (int i = 0; i < list.Count; i++)
                {
                    var (peptide, _, glycosylated) = list[i];
                    if (glycosylated is null || peptide.Parent.Accession != key.Accession
                        || key.Residue < peptide.OneBasedStartResidue || key.Residue > peptide.OneBasedEndResidue)
                    {
                        continue;
                    }
                    bool hasGlycan = glycosylated.Contains(key.Residue);
                    if (impliers.Contains(i))
                    {
                        if (hasGlycan != site.IsForbidden)
                            site.ImpliersAgreeing++;
                        else
                            site.ImpliersAgainst++;
                    }
                    else if (hasGlycan)
                    {
                        site.OthersWithGlycan++;
                    }
                    else
                    {
                        site.OthersWithout++;
                    }
                }
            }

            return sites.Values.Select(v => v.Site)
                .OrderBy(s => s.ProteinAccession).ThenBy(s => s.OneBasedResidue).ThenBy(s => s.IsForbidden).ToList();
        }
    }
}
