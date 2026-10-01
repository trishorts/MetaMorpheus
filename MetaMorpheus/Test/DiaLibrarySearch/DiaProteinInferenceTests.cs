#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer;
using EngineLayer.DiaLibrarySearch;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace Test.DiaLibrarySearch;

/// <summary>
/// Protein inference for the DIA search reuses MetaMorpheus's DDA engines unchanged (ProteinParsimonyEngine, then
/// ProteinScoringAndFdrEngine). The adapter's only job is to present each DIA peptide as the spectral match those
/// engines expect: mapped to every protein that yields it, carrying its peptide q-value and a positive score.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaProteinInferenceTests
{
    // Three proteins. LLLLLK is shared by A and B; EEEEER is unique to A, GGGGGR to B, and VVVVVK and TTTTTR to C.
    private static readonly List<Protein> Targets =
    [
        new("LLLLLKEEEEER", "A"),
        new("LLLLLKGGGGGR", "B"),
        new("VVVVVKTTTTTR", "C"),
    ];

    private static readonly CommonParameters Parameters = new(digestionParams: new DigestionParams(minPeptideLength: 5));

    private static List<Protein> TargetsAndDecoys() =>
        Targets.Concat(DecoyProteinGenerator.GenerateDecoys(Targets, DecoyType.Reverse)).ToList();

    private static Dictionary<string, List<PeptideWithSetModifications>> Index() =>
        DiaProteinInference.IndexPeptides(TargetsAndDecoys(), Parameters, [], []);

    private static DiaPeptideMatch Peptide(string sequence, double qValue, double score = 5, bool isDecoy = false, int index = 0) =>
        new(sequence, isDecoy, index, score, qValue);

    private static string DecoySequenceOf(string accession) =>
        Index().Where(e => e.Value.All(p => p.Parent.Accession == "DECOY_" + accession)).Select(e => e.Key)
            .OrderBy(s => s, StringComparer.Ordinal).First();

    [Test]
    public void APeptideMapsToEveryProteinThatYieldsIt()
    {
        var index = Index();

        Assert.That(index["LLLLLK"].Select(p => p.Parent.Accession), Is.EquivalentTo(new[] { "A", "B" }));
        Assert.That(index["EEEEER"].Select(p => p.Parent.Accession), Is.EquivalentTo(new[] { "A" }));
        Assert.That(index.Values.SelectMany(v => v).Count(p => p.Parent.IsDecoy), Is.GreaterThan(0), "decoy proteins are digested too");
    }

    /// <summary>
    /// Parsimony as in DDA: the shared peptide is explained by A, which has its own unique peptide, so B (seen only
    /// through the shared peptide) is not reported.
    /// </summary>
    [Test]
    public void ParsimonyExplainsASharedPeptideByTheProteinWithEvidenceOfItsOwn()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("LLLLLK", 0.001, index: 0),
            Peptide("EEEEER", 0.001, index: 1),
            Peptide("VVVVVK", 0.001, index: 2),
            Peptide("TTTTTR", 0.001, index: 3),
        };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        var names = result.ProteinGroups.Where(g => !g.IsDecoy).Select(g => g.ProteinGroupName).ToList();
        Assert.That(names, Is.EquivalentTo(new[] { "A", "C" }));
        var a = result.ProteinGroups.Single(g => g.ProteinGroupName == "A");
        Assert.That(a.AllPeptides.Select(p => p.FullSequence), Is.EquivalentTo(new[] { "LLLLLK", "EEEEER" }));
        Assert.That(a.UniquePeptides.Select(p => p.FullSequence), Is.EquivalentTo(new[] { "EEEEER" }));
        Assert.That(result.UnmappedPeptides, Is.EqualTo(0));
    }

    [Test]
    public void ADecoyPeptideMakesADecoyProteinGroup()
    {
        string decoy = DecoySequenceOf("C");
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("EEEEER", 0.001, score: 9, index: 0),
            Peptide(decoy, 0.002, score: 8, isDecoy: true, index: 1),
        };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        Assert.That(result.ProteinGroups.Single(g => g.IsDecoy).Proteins.Single().Accession, Is.EqualTo("DECOY_C"));
        Assert.That(result.ProteinGroups.Single(g => !g.IsDecoy).ProteinGroupName, Is.EqualTo("A"));
    }

    /// <summary>Only peptides that pass the q-value threshold count as evidence, exactly as DDA PSMs do.</summary>
    [Test]
    public void PeptidesAboveTheQValueThresholdAreNotEvidence()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("EEEEER", 0.001, index: 0),
            Peptide("VVVVVK", 0.5, index: 1),
        };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        Assert.That(result.ProteinGroups.Select(g => g.ProteinGroupName), Is.EquivalentTo(new[] { "A" }));
    }

    /// <summary>
    /// The DIA classifier's scores can be negative, and MetaMorpheus drops a spectral match scored at or below zero.
    /// A single order-preserving shift keeps every peptide, and the protein's best peptide is still the best scored.
    /// </summary>
    [Test]
    public void NegativeScoresAreShiftedNotDropped()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("EEEEER", 0.001, score: -1.5, index: 0),
            Peptide("LLLLLK", 0.001, score: -3.0, index: 1),
            Peptide("VVVVVK", 0.001, score: -7.25, index: 2),
        };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        Assert.That(result.SpectralMatches.Select(m => m.Score), Is.All.GreaterThan(0));
        Assert.That(result.SpectralMatches.OrderByDescending(m => m.Score).Select(m => m.FullSequence),
            Is.EqualTo(new[] { "EEEEER", "LLLLLK", "VVVVVK" }));
        Assert.That(result.ProteinGroups.Select(g => g.ProteinGroupName), Is.EquivalentTo(new[] { "A", "C" }));
    }

    [Test]
    public void EachSpectralMatchCarriesItsPeptideQValueAndPointsBackToItsPrecursor()
    {
        var peptides = new List<DiaPeptideMatch> { Peptide("EEEEER", 0.004, index: 41) };

        var match = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []).SpectralMatches.Single();

        Assert.That(match.FdrInfo.QValue, Is.EqualTo(0.004));
        Assert.That(match.FdrInfo.QValueNotch, Is.EqualTo(0.004));
        Assert.That(match.PeptideFdrInfo.QValue, Is.EqualTo(0.004));
        Assert.That(match.ScanNumber, Is.EqualTo(42), "one-based precursor index");
        Assert.That(match.FullFilePath, Is.EqualTo("run.raw"));
    }

    /// <summary>A peptide the database cannot yield is reported in a count, not silently lost.</summary>
    [Test]
    public void APeptideNoProteinYieldsIsCounted()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("EEEEER", 0.001, index: 0),
            Peptide("WWWWWK", 0.001, index: 1),
        };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        Assert.That(result.UnmappedPeptides, Is.EqualTo(1));
        Assert.That(result.SpectralMatches.Select(m => m.FullSequence), Is.EqualTo(new[] { "EEEEER" }));
    }

    /// <summary>
    /// A decoy peptide must land on a decoy protein, a target on a target. A sequence found only on the other side is
    /// unmapped rather than given a label it does not have.
    /// </summary>
    [Test]
    public void APeptideNeverLandsOnAProteinOfTheOtherKind()
    {
        var peptides = new List<DiaPeptideMatch> { Peptide("EEEEER", 0.001, isDecoy: true) };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        Assert.That(result.UnmappedPeptides, Is.EqualTo(1));
        Assert.That(result.ProteinGroups, Is.Empty);
    }

    // Peptide-level decoys, as a library made by reversing peptides has them: each pairs with one target sequence.
    private static readonly Dictionary<string, string> TargetOfDecoy = new()
    {
        ["KLLLLL"] = "LLLLLK",
        ["REEEEE"] = "EEEEER",
        ["KVVVVV"] = "VVVVVK",
        ["KWWWWW"] = "WWWWWK",
        ["RGGGGG"] = "GGGGGR",
    };

    private static DiaProteinInferenceResult InferWithInheritedDecoys(List<DiaPeptideMatch> peptides) =>
        DiaProteinInference.Infer(peptides, Targets, Parameters, "run.raw", [], [],
            decoy => TargetOfDecoy.GetValueOrDefault(decoy));

    /// <summary>
    /// A library whose decoys are reversed peptides has no decoy proteins to digest. Each decoy then takes its target's
    /// proteins, as decoys: the decoy of a shared peptide is shared by the same decoy proteins, so the decoy side has
    /// the target side's protein structure and protein FDR stays fair.
    /// </summary>
    [Test]
    public void ADecoyPeptideInheritsItsTargetsProteinsAsDecoys()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("KLLLLL", 0.001, isDecoy: true, index: 0),
            Peptide("REEEEE", 0.001, isDecoy: true, index: 1),
            Peptide("RGGGGG", 0.001, isDecoy: true, index: 2),
            Peptide("KVVVVV", 0.001, isDecoy: true, index: 3),
        };

        var result = InferWithInheritedDecoys(peptides);

        Assert.That(result.ProteinGroups.Select(g => g.ProteinGroupName), Is.EquivalentTo(new[] { "DECOY_A", "DECOY_B", "DECOY_C" }));
        foreach (string name in new[] { "DECOY_A", "DECOY_B" })
        {
            var group = result.ProteinGroups.Single(g => g.ProteinGroupName == name);
            Assert.That(group.AllPeptides.Select(p => p.FullSequence), Does.Contain("KLLLLL"), $"the shared decoy is in {name}");
            Assert.That(group.UniquePeptides.Select(p => p.FullSequence), Does.Not.Contain("KLLLLL"));
        }
        Assert.That(result.ProteinGroups, Is.All.Matches<ProteinGroup>(g => g.IsDecoy));
    }

    /// <summary>One decoy protein per target protein, whichever of its decoy peptides is seen; and with no other evidence, the shared decoy is explained by DECOY_A, as its target is by A.</summary>
    [Test]
    public void DecoysOfOneProteinShareOneDecoyProtein()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("KLLLLL", 0.001, isDecoy: true, index: 0),
            Peptide("REEEEE", 0.001, isDecoy: true, index: 1),
        };

        var result = InferWithInheritedDecoys(peptides);

        var a = result.ProteinGroups.Single(g => g.ProteinGroupName == "DECOY_A");
        Assert.That(a.Proteins, Has.Count.EqualTo(1));
        Assert.That(a.AllPeptides.Select(p => p.FullSequence), Is.EquivalentTo(new[] { "KLLLLL", "REEEEE" }));
    }

    [Test]
    public void ADecoyWhoseTargetNoProteinYieldsIsCounted()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("KWWWWW", 0.001, isDecoy: true, index: 0),
            Peptide("KNOPAIR", 0.001, isDecoy: true, index: 1),
            Peptide("EEEEER", 0.001, index: 2),
        };

        var result = InferWithInheritedDecoys(peptides);

        Assert.That(result.UnmappedPeptides, Is.EqualTo(2));
        Assert.That(result.ProteinGroups.Select(g => g.ProteinGroupName), Is.EquivalentTo(new[] { "A" }));
    }

    /// <summary>
    /// The DIA search has no PEP. A PEP q-value threshold would switch the protein engine to filtering on PEP, which
    /// every DIA peptide fails, and report no proteins at all. Refuse it instead.
    /// </summary>
    [Test]
    public void APepQValueThresholdIsRefused()
    {
        var pepFiltered = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5), pepQValueThreshold: 0.005);

        Assert.Throws<ArgumentException>(() =>
            DiaProteinInference.Infer([Peptide("EEEEER", 0.001)], TargetsAndDecoys(), pepFiltered, "run.raw", [], []));
    }
}
