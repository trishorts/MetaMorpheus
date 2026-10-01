#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using EngineLayer.FdrAnalysis;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>Protein groups inferred from one DIA run's peptides, and the spectral matches they were inferred from.</summary>
/// <param name="ProteinGroups">MetaMorpheus's scored protein groups, with its protein q-values, best first.</param>
/// <param name="SpectralMatches">One per mapped peptide, in the order of the peptides given.</param>
/// <param name="UnmappedPeptides">Peptides no protein of their own kind (target or decoy) yields; not used as evidence.</param>
public sealed record DiaProteinInferenceResult(List<ProteinGroup> ProteinGroups, List<SpectralMatch> SpectralMatches, int UnmappedPeptides);

/// <summary>
/// Protein inference for the DIA search, by MetaMorpheus's DDA engines unchanged: <see cref="ProteinParsimonyEngine"/>
/// then <see cref="ProteinScoringAndFdrEngine"/>, as PostSearchAnalysisTask runs them. This class is only the boundary:
/// it presents each DIA peptide as the spectral match those engines expect.
/// </summary>
/// <remarks>
/// One spectral match per peptide, not per precursor: the peptide carries a q-value for targets and decoys alike, which
/// protein FDR needs. Proteins come from digesting the database with the search's digestion parameters, as in DDA, so
/// the library's decoys must be digests of decoy proteins (as the library builder makes them).
/// </remarks>
public static class DiaProteinInference
{
    /// <summary>Every peptide the proteins yield under these digestion parameters, keyed by full sequence.</summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static Dictionary<string, List<PeptideWithSetModifications>> IndexPeptides(IEnumerable<Protein> proteins,
        CommonParameters commonParameters, List<Modification> fixedModifications, List<Modification> variableModifications)
    {
        ArgumentNullException.ThrowIfNull(proteins);
        ArgumentNullException.ThrowIfNull(commonParameters);
        ArgumentNullException.ThrowIfNull(fixedModifications);
        ArgumentNullException.ThrowIfNull(variableModifications);

        var index = new Dictionary<string, List<PeptideWithSetModifications>>();
        foreach (var protein in proteins)
        {
            foreach (var peptide in protein.Digest(commonParameters.DigestionParams, fixedModifications, variableModifications).Cast<PeptideWithSetModifications>())
            {
                if (!index.TryGetValue(peptide.FullSequence, out var list))
                {
                    index[peptide.FullSequence] = list = [];
                }
                list.Add(peptide);
            }
        }
        return index;
    }

    /// <summary>
    /// Infers and scores protein groups from one run's peptides (<see cref="DiaPeptideFdr.Assign"/>), targets and decoys.
    /// Peptides passing <see cref="CommonParameters.QValueThreshold"/> are the evidence, as PSMs are in DDA.
    /// </summary>
    /// <param name="proteins">Target and decoy proteins: the database the library was built from.</param>
    /// <param name="fullFilePath">The run the peptides came from.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">
    /// A PEP q-value threshold is set. The DIA search has no PEP, so the protein engine would filter every peptide out.
    /// </exception>
    public static DiaProteinInferenceResult Infer(IReadOnlyList<DiaPeptideMatch> peptides, IEnumerable<Protein> proteins,
        CommonParameters commonParameters, string fullFilePath, List<Modification> fixedModifications, List<Modification> variableModifications)
    {
        ArgumentNullException.ThrowIfNull(peptides);
        ArgumentNullException.ThrowIfNull(fullFilePath);
        if (commonParameters is not null && commonParameters.PepQValueThreshold < commonParameters.QValueThreshold)
        {
            throw new ArgumentException("DIA protein inference filters on q-value; the DIA search computes no PEP, so a PEP q-value threshold is not supported.", nameof(commonParameters));
        }

        var index = IndexPeptides(proteins, commonParameters!, fixedModifications, variableModifications);
        var matches = BuildSpectralMatches(peptides, index, commonParameters!, fullFilePath, out int unmapped);

        var parsimony = (ProteinParsimonyResults)new ProteinParsimonyEngine(matches, modPeptidesAreDifferent: false,
            commonParameters!, null, []).Run();
        var scored = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(parsimony.ProteinGroups, matches,
            noOneHitWonders: false, treatModPeptidesAsDifferentPeptides: false, mergeIndistinguishableProteinGroups: true,
            commonParameters!, null, []).Run();

        return new DiaProteinInferenceResult(scored.SortedAndScoredProteinGroups, matches, unmapped);
    }

    /// <summary>
    /// One spectral match per peptide that a protein of its own kind yields, attached to every such protein. Scores are
    /// shifted by one constant so the lowest is 1: MetaMorpheus ignores a match scored at or below zero.
    /// </summary>
    private static List<SpectralMatch> BuildSpectralMatches(IReadOnlyList<DiaPeptideMatch> peptides,
        Dictionary<string, List<PeptideWithSetModifications>> index, CommonParameters commonParameters, string fullFilePath, out int unmapped)
    {
        double shift = peptides.Count == 0 ? 0 : Math.Max(0, 1 - peptides.Min(p => p.Score));
        var matches = new List<SpectralMatch>();
        unmapped = 0;
        foreach (var peptide in peptides)
        {
            var yields = index.TryGetValue(peptide.FullSequence, out var all)
                ? all.Where(p => p.Parent.IsDecoy == peptide.IsDecoy).ToList()
                : [];
            if (yields.Count == 0)
            {
                unmapped++;
                continue;
            }

            double score = peptide.Score + shift;
            var scan = SyntheticScan(peptide, yields[0], commonParameters, fullFilePath);
            var match = new PeptideSpectralMatch(yields[0], 0, score, peptide.BestPrecursorIndex, scan, commonParameters, new List<MatchedFragmentIon>());
            foreach (var other in yields.Skip(1))
            {
                match.AddOrReplace(other, score, 0, true, new List<MatchedFragmentIon>());
            }
            match.ResolveAllAmbiguities();
            match.SetFdrValues(0, 0, peptide.QValue, 0, 0, peptide.QValue, double.NaN, pepQValue: 2);
            match.PeptideFdrInfo = new FdrInfo { QValue = peptide.QValue, QValueNotch = peptide.QValue, PEP = double.NaN, PEP_QValue = 2 };
            matches.Add(match);
        }
        return matches;
    }

    /// <summary>
    /// A stand-in scan: a DIA peptide is a chromatographic peak group, not one scan. Its scan number is the one-based
    /// library index of the peptide's best precursor, so each match points back to the precursor it came from.
    /// </summary>
    private static Ms2ScanWithSpecificMass SyntheticScan(DiaPeptideMatch peptide, PeptideWithSetModifications sequence,
        CommonParameters commonParameters, string fullFilePath)
    {
        int scanNumber = peptide.BestPrecursorIndex + 1;
        var scan = new MsDataScan(new MzSpectrum([1.0], [1.0], false), scanNumber, 2, true, Polarity.Positive, double.NaN,
            new MzRange(0, 1), null, MZAnalyzerType.Orbitrap, double.NaN, null, null, $"scan={scanNumber}", double.NaN, null,
            null, double.NaN, null, DissociationType.HCD, null, null);
        const int charge = 2;
        return new Ms2ScanWithSpecificMass(scan, sequence.MonoisotopicMass.ToMz(charge), charge, fullFilePath, commonParameters,
            Array.Empty<IsotopicEnvelope>());
    }
}
