#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Security.Cryptography;
using EngineLayer;
using Chromatography.RetentionTimeCalibration;
using EngineLayer.DiaLibrarySearch;
using MzLibUtil;
using NUnit.Framework;
using Readers.SpectralLibrary;

namespace Test.DiaLibrarySearch;

/// <summary>
/// The DIA library search's permanent contract, run on <see cref="SyntheticDiaRun"/>, where the answer is known. Every
/// later milestone deepens the engine; none may break these.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaLibrarySearchEngineTests
{
    private string _directory = "";

    /// <summary>
    /// A linear run-RT to iRT map through the true curve's endpoints (iRT -20 at 0.5 min, iRT 120 at 9.936 min). The
    /// true curve bows away from it by up to about 12 iRT, inside the engine's default iRT window.
    /// </summary>
    private static readonly IrtCalibrationModel EndpointMap = IrtCalibration.Line(
        (new RtMinutes(SyntheticDiaRun.TrueRtMinutes(-20)), new Irt(-20)),
        (new RtMinutes(SyntheticDiaRun.TrueRtMinutes(120)), new Irt(120)));

    [OneTimeSetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DiaLibrarySearchEngineTests");
        Directory.CreateDirectory(_directory);
    }

    [OneTimeTearDown]
    public void TearDown()
    {
        Directory.Delete(_directory, true);
    }

    private DiaLibrarySearchResults Search(SyntheticDiaRun run, IrtCalibrationModel map, out string libraryPath)
    {
        libraryPath = run.WriteLibrary(_directory);
        using var library = MslLibrary.Load(libraryPath);
        var engine = new DiaLibrarySearchEngine(run.Scans, library, map, new DiaLibrarySearchParameters(),
            new CommonParameters(), [], []);
        return (DiaLibrarySearchResults)engine.Run();
    }

    private DiaLibrarySearchResults Search(SyntheticDiaRun run) => Search(run, EndpointMap, out _);

    /// <summary>Three of every four targets elute. At least 90% of those are reported at 1% FDR.</summary>
    [Test]
    public void PlantedPrecursorsAreFound()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);

        var results = Search(run);

        var found = results.Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01).Select(m => m.FullSequence).ToHashSet();
        int planted = run.PlantedSequences.Count;
        Assert.That(planted, Is.GreaterThan(100), "the fixture must plant enough precursors to mean something");
        Assert.That(found.Count(run.PlantedSequences.Contains), Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * planted)));
        Assert.That(found.Count(sequence => !run.PlantedSequences.Contains(sequence)), Is.LessThanOrEqualTo(Math.Max(1, found.Count / 100)),
            "a target that never eluted is a false discovery; at 1% FDR there can be at most about 1% of them");
    }

    /// <summary>A run of noise alone has nothing to find.</summary>
    [Test]
    public void NullRunReportsNothing()
    {
        var run = SyntheticDiaRun.Build(200, _ => false);

        var results = Search(run);

        Assert.That(results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01), Is.LessThanOrEqualTo(Math.Max(1, results.Matches.Count / 100)));
    }

    /// <summary>When only the decoys elute, no target may pass: every real signal belongs to a decoy.</summary>
    [Test]
    public void ARunOfDecoysReportsNoTargets()
    {
        var run = SyntheticDiaRun.Build(200, entry => entry.IsDecoy);

        var results = Search(run);

        Assert.That(results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01), Is.Zero);
    }

    /// <summary>With no decoys, target-decoy q-values are all zero, which reads as "every target is right".</summary>
    [Test]
    public void ALibraryWithoutDecoysThrows()
    {
        var run = SyntheticDiaRun.Build(50, entry => true, withDecoys: false);

        var e = Assert.Throws<MetaMorpheusException>(() => Search(run));
        Assert.That(e!.Message, Does.Contain("decoy"));
    }

    /// <summary>
    /// A sampled search, as calibration's first pass uses, scores only precursors whose library index is a multiple of the
    /// stride. It samples targets and decoys alike, so the target-decoy competition stays fair.
    /// </summary>
    [Test]
    public void ASampledSearchScoresOnlyEveryStrideThPrecursor()
    {
        // Decoys elute too: an entry with no signal at all is never scored, and the test needs both classes
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var full = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(), new CommonParameters(), [], []).Run();
        var sampled = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(PrecursorSampleStride: 3), new CommonParameters(), [], []).Run();

        Assert.That(sampled.Matches.All(m => m.PrecursorIndex % 3 == 0));
        Assert.That(sampled.TargetCount, Is.EqualTo(full.Matches.Count(m => !m.IsDecoy && m.PrecursorIndex % 3 == 0)));
        Assert.That(sampled.DecoyCount, Is.EqualTo(full.Matches.Count(m => m.IsDecoy && m.PrecursorIndex % 3 == 0)));
        Assert.That(sampled.DecoyCount, Is.GreaterThan(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(PrecursorSampleStride: 0));
    }

    /// <summary>
    /// An offset picks which of the stride's samples is scored, so calibration can be repeated on disjoint samples: a
    /// different first-pass sample cost HF-X 20% of its precursors, and the spread over samples is what to measure.
    /// </summary>
    [Test]
    public void AnOffsetPicksADisjointSample()
    {
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var samples = Enumerable.Range(0, 3).Select(offset => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(PrecursorSampleStride: 3) { PrecursorSampleOffset = offset }, new CommonParameters(), [], []).Run()).Matches).ToList();

        for (int offset = 0; offset < 3; offset++)
            Assert.That(samples[offset].All(m => m.PrecursorIndex % 3 == offset), $"offset {offset}");
        Assert.That(samples.Sum(s => s.Count), Is.EqualTo(((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(), new CommonParameters(), [], []).Run()).Matches.Count));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(PrecursorSampleStride: 3) { PrecursorSampleOffset = 3 });
    }

    /// <summary>The search reports peptides as well as precursors, one per full sequence, with peptide-level q-values.</summary>
    [Test]
    public void PeptidesAreReportedWithTheirOwnQValues()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);

        var results = Search(run);

        Assert.That(results.Peptides.Select(p => (p.FullSequence, p.IsDecoy)).Distinct().Count(), Is.EqualTo(results.Peptides.Count));
        Assert.That(results.Peptides, Is.EquivalentTo(DiaPeptideFdr.Assign(results.Matches)));
        int found = results.Peptides.Count(p => !p.IsDecoy && p.QValue <= 0.01 && run.PlantedSequences.Contains(p.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
        Assert.That(results.ToString(), Does.Contain("Target peptides with q-value <= 0.01:"));
    }

    /// <summary>
    /// Each precursor carries a quantity: its fragments' summed peak areas. It follows abundance, so a precursor planted
    /// at twice the height quantifies at about twice the amount, and one at equal height at about the same.
    /// </summary>
    [Test]
    public void QuantityFollowsAbundance()
    {
        Func<Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry, bool> plant = entry => !entry.IsDecoy;
        var baseline = Search(SyntheticDiaRun.Build(100, plant)).Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01)
            .ToDictionary(m => m.PrecursorIndex);
        var doubled = Search(SyntheticDiaRun.Build(100, plant, abundance: e => SyntheticDiaRun.Bucket(e, 2) == 0 ? 2 : 1))
            .Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01).ToList();

        var ratios = doubled.Where(m => baseline.ContainsKey(m.PrecursorIndex))
            .Select(m => (Doubled: SyntheticDiaRun.Bucket(m.FullSequence, 2) == 0, Ratio: m.Quantity / baseline[m.PrecursorIndex].Quantity))
            .ToList();

        Assert.That(baseline.Values.Select(m => m.Quantity), Is.All.GreaterThan(0));
        Assert.That(ratios.Where(r => r.Doubled).Select(r => r.Ratio).ToList(), Has.Count.GreaterThan(20).And.All.InRange(1.8, 2.2));
        Assert.That(ratios.Where(r => !r.Doubled).Select(r => r.Ratio).ToList(), Has.Count.GreaterThan(20).And.All.InRange(0.9, 1.1));
    }

    /// <summary>
    /// Library context as features (DIA-NN feeds m/z, charge and fragment count to its classifier), so the classifier can
    /// learn when to trust a prediction. A decoy shares its target's precursor m/z and charge, so these cannot reveal the label.
    /// </summary>
    [Test]
    public void PrecursorMzChargeAndLibraryFragmentCountAreFeatures()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000,
            withCharge3Siblings: true);
        var entryOf = run.Library.ToDictionary(e => (e.FullSequence, e.ChargeState));

        var matches = Search(run).Matches;

        int mz = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "PrecursorMz");
        int charge = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "PrecursorCharge");
        int count = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "LibraryFragmentCount");
        Assert.That(new[] { mz, charge, count }, Is.All.GreaterThanOrEqualTo(0));
        Assert.That(matches.Select(m => m.Charge).Distinct().Count(), Is.EqualTo(2));
        foreach (var m in matches)
        {
            Assert.That(m.Features[mz], Is.EqualTo(m.PrecursorMz).Within(1e-3));
            Assert.That(m.Features[charge], Is.EqualTo(m.Charge));
            Assert.That(m.Features[count], Is.EqualTo(entryOf[(m.FullSequence, m.Charge)].MatchedFragmentIons.Count));
        }
    }

    /// <summary>
    /// Charge-state siblings (AlphaDIA scores elution groups): a precursor whose other charge state co-elutes at the same
    /// apex has independent evidence. The features use only the sequence and the decoy flag, never which one is a target,
    /// and decoys have sibling pairs too. A precursor with no sibling in the library gets none.
    /// </summary>
    [Test]
    public void ACoElutingChargeStateSiblingIsEvidence()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000,
            withCharge3Siblings: true);
        var withSibling = run.Library.GroupBy(e => (e.FullSequence, e.IsDecoy)).Where(g => g.Count() > 1).Select(g => g.Key).ToHashSet();

        var matches = Search(run).Matches;

        int sibling = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "SiblingCoElution");
        int delta = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "SiblingApexDeltaMinutes");
        Assert.That(new[] { sibling, delta }, Is.All.GreaterThanOrEqualTo(0));
        double Median(IEnumerable<double> values) { var s = values.Order().ToArray(); return s[s.Length / 2]; }
        var plantedPairs = matches.Where(m => run.PlantedSequences.Contains(m.FullSequence) && withSibling.Contains((m.FullSequence, m.IsDecoy))).ToList();
        var chancePairs = matches.Where(m => !run.PlantedSequences.Contains(m.FullSequence) && withSibling.Contains((m.FullSequence, m.IsDecoy))).ToList();
        var alone = matches.Where(m => !withSibling.Contains((m.FullSequence, m.IsDecoy))).ToList();
        Assert.That(plantedPairs, Has.Count.GreaterThan(30));
        Assert.That(chancePairs, Has.Count.GreaterThan(30));
        Assert.That(chancePairs.Count(m => m.IsDecoy), Is.GreaterThan(10), "decoys have siblings too");
        Assert.That(Median(plantedPairs.Select(m => m.Features[sibling])), Is.GreaterThan(0.8));
        Assert.That(Median(plantedPairs.Select(m => m.Features[delta])), Is.LessThan(0.05));
        Assert.That(Median(chancePairs.Select(m => m.Features[sibling])), Is.LessThan(0.5));
        Assert.That(alone.Select(m => m.Features[sibling]), Is.All.EqualTo(0));
        Assert.That(alone.Select(m => m.Features[delta]), Is.All.EqualTo(1.0), "no sibling counts as a full minute apart");
    }

    /// <summary>
    /// MS1 evidence that does not lean on the fragments. On the whole-proteome library, a fifth of DIA-NN's identifications we
    /// miss show a clean MS1 peak inside DIA-NN's bounds while their fragments are faint, and our only MS1 features were
    /// correlations to the fragment profile. At the apex's MS1 scan: the M0-M3 envelope against the expected isotope pattern,
    /// the M0 mass error, and how much of the window's MS1 trace maximum the apex holds. A decoy has its target's precursor,
    /// so none of these can reveal the label.
    /// </summary>
    [Test]
    public void Ms1EnvelopeMassErrorAndApexShareSeparatePlantedPrecursorsFromChance()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);

        var matches = Search(run).Matches;

        int envelope = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1EnvelopeCosine");
        int ppm = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1MassErrorPpm");
        int share = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1ApexShare");
        Assert.That(new[] { envelope, ppm, share }, Is.All.GreaterThanOrEqualTo(0));
        double Median(IEnumerable<double> values) { var s = values.Order().ToArray(); return s[s.Length / 2]; }
        var planted = matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).ToList();
        var chance = matches.Where(m => !run.PlantedSequences.Contains(m.FullSequence)).ToList();
        Assert.That(planted, Has.Count.GreaterThan(50));
        Assert.That(chance, Has.Count.GreaterThan(50));
        Assert.That(Median(planted.Select(m => m.Features[envelope])), Is.GreaterThan(0.95));
        Assert.That(Median(chance.Select(m => m.Features[envelope])), Is.LessThan(0.5));
        Assert.That(Median(planted.Select(m => m.Features[ppm])), Is.LessThan(1), "the planted M0 is at the library m/z");
        Assert.That(Median(chance.Select(m => m.Features[ppm])), Is.EqualTo(20).Within(1e-9), "no M0 found counts as the full tolerance");
        Assert.That(Median(planted.Select(m => m.Features[share])), Is.GreaterThan(0.8));
        Assert.That(matches.Select(m => m.Features[envelope]), Is.All.InRange(0.0, 1.0 + 1e-12));
        Assert.That(matches.Select(m => m.Features[share]), Is.All.InRange(0.0, 1.0));
    }

    /// <summary>
    /// Interference removal drops matches from the reported list, but the results keep them, so a miss can be told
    /// apart from a precursor never scored. Reported and removed together are every precursor a search without removal
    /// reports.
    /// </summary>
    [Test]
    public void MatchesRemovedAsInterferenceAreKeptInTheResults()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);
        string path = run.WriteLibrary(_directory);
        DiaLibrarySearchResults SearchWith(int explained)
        {
            using var library = MslLibrary.Load(path);
            return (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
                new DiaLibrarySearchParameters(InterferenceExplainedFragments: explained), new CommonParameters(), [], []).Run();
        }

        var with = SearchWith(4);
        var without = SearchWith(0);

        Assert.That(without.RemovedAsInterference, Is.Empty);
        Assert.That(with.Matches.Select(m => m.PrecursorIndex).Intersect(with.RemovedAsInterference.Select(m => m.PrecursorIndex)), Is.Empty);
        Assert.That(with.Matches.Concat(with.RemovedAsInterference).Select(m => m.PrecursorIndex),
            Is.EquivalentTo(without.Matches.Select(m => m.PrecursorIndex)));
    }

    /// <summary>
    /// The classifier's network can train on a random subsample of each fold (mzLib's maxNetworkTrainingRows): on a
    /// whole-proteome library, training on every row was two thirds of the search. A cap above the row count changes
    /// nothing; a small one reaches the rescorer.
    /// </summary>
    [Test]
    public void TheNetworkTrainingCapReachesTheClassifier()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);
        string path = run.WriteLibrary(_directory);
        List<DiaPrecursorMatch> SearchWith(int? cap)
        {
            using var library = MslLibrary.Load(path);
            return ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
                new DiaLibrarySearchParameters(MaxNetworkTrainingRows: cap), new CommonParameters(), [], []).Run()).Matches;
        }

        var all = SearchWith(null);

        Assert.That(SearchWith(10_000_000).Select(m => m.Score), Is.EqualTo(all.Select(m => m.Score)));
        Assert.That(SearchWith(50).Select(m => m.Score), Is.Not.EqualTo(all.Select(m => m.Score)));
    }

    /// <summary>
    /// Peptide length is a feature: a chance match is harder to make with more residues (more fragments, more of them at
    /// sparse high m/z), and on the whole-proteome library the decoys scoring like real identifications were shorter
    /// (median 9 residues against 11). A reversed decoy has exactly its target's length, so the feature cannot reveal the
    /// label.
    /// </summary>
    [Test]
    public void PeptideLengthIsAFeatureAndADecoyHasItsTargetsLength()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);
        var lengthOf = run.Library.ToDictionary(e => e.FullSequence, e => e.BaseSequence.Length);

        var matches = Search(run).Matches;

        int length = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "PeptideLength");
        Assert.That(length, Is.GreaterThanOrEqualTo(0));
        Assert.That(matches, Has.Count.GreaterThan(100));
        Assert.That(matches.All(m => m.Features[length] == lengthOf[m.FullSequence]));
        Assert.That(matches.Where(m => m.IsDecoy).All(m => m.Features[length] == lengthOf[new string(m.FullSequence.Reverse().ToArray())]));
    }

    /// <summary>
    /// Co-elution is measured against the smoothed best of the six most intense library fragments. A planted precursor's
    /// fragments share one elution profile, so its co-elution is near 1. A candidate seeing only noise has none.
    /// </summary>
    [Test]
    public void CoElutionMeasuresAgreementWithTheBestFragment()
    {
        // Dense noise, so candidates that were never planted still match fragments by chance and get scored
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);

        var matches = Search(run).Matches;

        int coElution = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "CoElution");
        double Median(IEnumerable<double> values) { var s = values.Order().ToArray(); return s[s.Length / 2]; }
        var planted = matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[coElution]).ToList();
        var noise = matches.Where(m => !run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[coElution]).ToList();
        Assert.That(planted, Has.Count.GreaterThan(50));
        Assert.That(noise, Has.Count.GreaterThan(50), "the fixture must score chance matches");
        Assert.That(Median(planted), Is.GreaterThan(0.9));
        Assert.That(Median(noise), Is.LessThan(0.5));
        Assert.That(matches.All(m => m.Features[coElution] is >= 0 and <= 1));
    }

    /// <summary>
    /// The search runs candidates in parallel; the result must not depend on the thread count, match for match, in
    /// order, including every feature, score and q-value.
    /// </summary>
    [Test]
    public void ParallelSearchGivesTheSameResultAsOneThread()
    {
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 600);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults Run(int threads) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(), new CommonParameters(maxThreadsToUsePerFile: threads), [], []).Run();

        var one = Run(1).Matches;
        var four = Run(4).Matches;

        Assert.That(four, Has.Count.EqualTo(one.Count));
        for (int i = 0; i < one.Count; i++)
        {
            Assert.That(four[i].PrecursorIndex, Is.EqualTo(one[i].PrecursorIndex), $"row {i}");
            Assert.That(four[i].Score, Is.EqualTo(one[i].Score), $"row {i}");
            Assert.That(four[i].QValue, Is.EqualTo(one[i].QValue), $"row {i}");
            Assert.That(four[i].Features, Is.EqualTo(one[i].Features), $"row {i}");
        }
    }

    [Test]
    public void TargetQValuesNeverDecreaseAsScoreFalls()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0);

        var targets = Search(run).Matches.Where(m => !m.IsDecoy).OrderByDescending(m => m.Score).ToList();

        Assert.That(targets, Is.Not.Empty);
        for (int i = 1; i < targets.Count; i++)
            Assert.That(targets[i].QValue, Is.GreaterThanOrEqualTo(targets[i - 1].QValue), $"at rank {i}");
        Assert.That(targets.All(m => m.QValue is >= 0 and <= 1));
    }

    [Test]
    public void TargetAndDecoyCountsAreReported()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);

        var results = Search(run);

        int targets = results.Matches.Count(m => !m.IsDecoy);
        int decoys = results.Matches.Count(m => m.IsDecoy);
        int passing = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01);
        Assert.That(results.TargetCount, Is.EqualTo(targets));
        Assert.That(results.DecoyCount, Is.EqualTo(decoys));
        Assert.That(results.ToString(), Does.Contain($"Target precursors: {targets}"));
        Assert.That(results.ToString(), Does.Contain($"Decoy precursors: {decoys}"));
        Assert.That(results.ToString(), Does.Contain($"Target precursors with q-value <= 0.01: {passing}"));
    }

    /// <summary>
    /// iRT and run minutes are different quantities. The types keep them apart, with no conversion in either
    /// direction, not even to double, so the only way across is a map. A search told the identity map (iRT = minutes)
    /// looks in the wrong place and finds almost nothing.
    /// </summary>
    [Test]
    public void IrtAndMinutesDoNotMix()
    {
        foreach (var type in new[] { typeof(Irt), typeof(RtMinutes) })
        {
            var conversions = type.GetMethods(BindingFlags.Public | BindingFlags.Static)
                .Where(method => method.Name is "op_Implicit" or "op_Explicit");
            Assert.That(conversions, Is.Empty, $"{type.Name} must not convert implicitly or explicitly");
        }

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        var identity = IrtCalibration.Line((new RtMinutes(0), new Irt(0)), (new RtMinutes(1), new Irt(1)));

        var results = Search(run, identity, out _);

        int found = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.LessThan(run.PlantedSequences.Count / 2));
    }

    /// <summary>The run-to-iRT map is applied to the run. The library is never rewritten in minutes.</summary>
    [Test]
    public void TheLibraryFileIsUnchangedBySearch()
    {
        var run = SyntheticDiaRun.Build(50, entry => !entry.IsDecoy);
        string path = run.WriteLibrary(_directory);
        byte[] before = SHA256.HashData(File.ReadAllBytes(path));

        using (var library = MslLibrary.Load(path))
            new DiaLibrarySearchEngine(run.Scans, library, EndpointMap, new DiaLibrarySearchParameters(), new CommonParameters(), [], []).Run();

        Assert.That(SHA256.HashData(File.ReadAllBytes(path)), Is.EqualTo(before));
    }

    [Test]
    public void ReportedApexIsInIrtAndInMinutes()
    {
        var run = SyntheticDiaRun.Build(100, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);

        var best = Search(run).Matches.Where(m => !m.IsDecoy && run.PlantedSequences.Contains(m.FullSequence))
            .OrderByDescending(m => m.Score).First();

        var entry = run.Library.Single(e => e.FullSequence == best.FullSequence);
        Assert.That(best.LibraryIrt.Value, Is.EqualTo(entry.RetentionTime).Within(1e-3));
        Assert.That(best.ApexRt.Value, Is.EqualTo(SyntheticDiaRun.TrueRtMinutes(entry.RetentionTime)).Within(SyntheticDiaRun.CycleMinutes));
        Assert.That(best.ApexIrt, Is.EqualTo(EndpointMap.ToIrt(best.ApexRt)));
    }
}
