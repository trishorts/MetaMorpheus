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

    /// <summary>
    /// The classifier's network can train twice, as DIA-NN's does: the second pass learns from the candidate peaks the first
    /// network picked. The setting reaches the classifier: the same precursors are scored, differently. Whether it helps is
    /// measured on real data (results/entrapment/README.md); the rescorer's own tests show the mechanism.
    /// </summary>
    [Test]
    public void TheClassifierNetworkCanTrainTwice()
    {
        // Decoys elute too, so the classifier has both classes to train on
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(int passes) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierNetworkPasses: passes), new CommonParameters(), [], []).Run();

        var once = SearchWith(1);
        var twice = SearchWith(2);

        Assert.That(twice.Matches.Select(m => m.Score), Is.Not.EqualTo(once.Matches.Select(m => m.Score)));
        Assert.That(twice.Matches.Select(m => m.PrecursorIndex), Is.EquivalentTo(once.Matches.Select(m => m.PrecursorIndex)));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(ClassifierNetworkPasses: 0));
    }

    /// <summary>
    /// The classifier's network can train on the most confident rows rather than a random sample, as DIA-NN trains after
    /// removing low-confidence identifications. The setting reaches the classifier: with a cap below the rows, the scores change.
    /// </summary>
    [Test]
    public void TheClassifierNetworkCanTrainOnItsMostConfidentRows()
    {
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(StatisticalModels.NetworkTrainingSample sample) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans,
            library, EndpointMap, new DiaLibrarySearchParameters(MaxNetworkTrainingRows: 40, ClassifierNetworkTrainingSample: sample),
            new CommonParameters(), [], []).Run();

        var random = SearchWith(StatisticalModels.NetworkTrainingSample.Random);
        var confident = SearchWith(StatisticalModels.NetworkTrainingSample.Confident);

        Assert.That(confident.Matches.Select(m => m.PrecursorIndex), Is.EquivalentTo(random.Matches.Select(m => m.PrecursorIndex)));
        Assert.That(confident.Matches.Select(m => m.Score), Is.Not.EqualTo(random.Matches.Select(m => m.Score)));
    }

    /// <summary>
    /// The confident training sample is the default: at a matched entrapment FDP of 1% it found 9.4% more precursors on
    /// PXD022589 (HF-X) and 9.9% more on PXD005573 1 h than a random sample, at no extra time.
    /// </summary>
    [Test]
    public void TheNetworkTrainsOnItsMostConfidentRowsByDefault() =>
        Assert.That(new DiaLibrarySearchParameters().ClassifierNetworkTrainingSample, Is.EqualTo(StatisticalModels.NetworkTrainingSample.Confident));

    /// <summary>
    /// A precursor's other candidate peaks are kept, with their apex and score, so that a miss where another engine chose a
    /// different peak can be told apart: was that peak among our candidates and scored lower, or never a candidate?
    /// </summary>
    [Test]
    public void APrecursorsLosingCandidatePeaksAreKept()
    {
        // Dense noise gives precursors weaker second and third candidate peaks
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 2000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(MaxApexCandidates: 3), new CommonParameters(), [], []).Run();

        Assert.That(results.LosingCandidates.Values.Sum(v => v.Count), Is.GreaterThan(0));
        foreach (var match in results.Matches)
        {
            if (!results.LosingCandidates.TryGetValue(match.PrecursorIndex, out var others))
                continue;
            Assert.That(others, Has.Count.LessThanOrEqualTo(2));
            Assert.That(others.All(o => o.Score <= match.Score && o.ApexRt.Value != match.ApexRt.Value), $"precursor {match.PrecursorIndex}");
        }
    }

    /// <summary>
    /// Sibling support can come from each other charge state's top candidate peak only. Taken from any of a sibling's
    /// candidates, a decoy gets several chances for a noise peak to sit near its apex; at equal score, decoys on HF-X had
    /// more sibling support than the targets we missed (results/entrapment/README.md, I14).
    /// </summary>
    [Test]
    public void SiblingSupportCanComeFromTheSiblingsTopCandidateOnly()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 2000,
            withCharge3Siblings: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        int sibling = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "SiblingCoElution");
        DiaLibrarySearchResults SearchWith(bool topOnly) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(SiblingTopCandidateOnly: topOnly), new CommonParameters(), [], []).Run();

        var any = SearchWith(false).Matches.ToDictionary(m => m.PrecursorIndex, m => m.Features[sibling]);
        var top = SearchWith(true).Matches.ToDictionary(m => m.PrecursorIndex, m => m.Features[sibling]);

        Assert.That(top.Keys.Intersect(any.Keys).Count(k => top[k] != any[k]), Is.GreaterThan(0));
        // The default since 2026-10-02: HF-X +3.5% at a matched entrapment FDP of 1%, PXD005573 1 h unchanged (-0.3%)
        Assert.That(new DiaLibrarySearchParameters().SiblingTopCandidateOnly, Is.True);
    }

    /// <summary>
    /// DIA-NN 1.8 scores we lacked. MinCorr is co-elution of spike-suppressed traces (each scan the minimum of itself and its
    /// neighbours). NFCorr is the unfragmented precursor m/z in MS2. ShadowCorr is the trace one isotope below each fragment:
    /// is "our" fragment another ion's M+1? All three are 0 unless asked for. A planted precursor's broad peak survives spike
    /// suppression, so its MinCorr is high.
    /// </summary>
    [Test]
    public void DiaNnScoresAreFilledOnlyWhenAskedFor()
    {
        Assert.That(DiaLibrarySearchEngine.SpikeSuppressed([0, 0, 9, 0, 0, 1, 2, 3, 2, 1]),
            Is.EqualTo(new double[] { 0, 0, 0, 0, 0, 0, 1, 2, 1, 1 }));
        int[] scores = new[] { "MinCorr", "NFCorr", "ShadowCorr" }.Select(n => Array.IndexOf(DiaPrecursorMatch.FeatureNames, n)).ToArray();
        Assert.That(scores, Has.All.GreaterThanOrEqualTo(0));

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(DiaNnScores: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.SelectMany(m => scores.Select(i => m.Features[i])), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence)).ToList();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted.Average(m => m.Features[scores[0]]), Is.GreaterThan(0.5));
        // The default since 2026-10-04: with DIA-NN peak finding at 6 candidates, +1.4% HF-X, +3.6% PXD005573 (matched paired FDP 1%)
        Assert.That(new DiaLibrarySearchParameters().DiaNnScores, Is.True);
    }

    /// <summary>
    /// DIA-NN's pSig: each of the six most intense library fragments' share of their summed signal across the peak, in library
    /// rank, so the classifier can set the observed pattern against the library's. 0 unless asked for. On, a planted
    /// precursor's shares sum to 1 and follow the library: the most intense fragment takes more than the sixth.
    /// </summary>
    [Test]
    public void FragmentSignalSharesAreFilledOnlyWhenAskedFor()
    {
        int[] shares = Enumerable.Range(1, 6).Select(k => Array.IndexOf(DiaPrecursorMatch.FeatureNames, $"SignalShare{k}")).ToArray();
        Assert.That(shares, Has.All.GreaterThanOrEqualTo(0));

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(DiaNnSignalShare: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.SelectMany(m => shares.Select(i => m.Features[i])), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence)).ToList();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted.Select(m => shares.Sum(i => m.Features[i])), Has.All.EqualTo(1).Within(1e-9));
        Assert.That(planted.Average(m => m.Features[shares[0]]), Is.GreaterThan(planted.Average(m => m.Features[shares[5]])));
        Assert.That(new DiaLibrarySearchParameters().DiaNnSignalShare, Is.False);
    }

    /// <summary>
    /// The linear model that picks each precursor's candidate peak for the network is refit a set number of times (DIA-NN about
    /// 8). 3 by default (a best-single-feature ranking, then two refits); fewer than 1 is refused. More refits still find
    /// what was planted.
    /// </summary>
    [Test]
    public void TheLinearPickCanBeRefitMoreTimes()
    {
        Assert.That(new DiaLibrarySearchParameters().ClassifierLinearIterations, Is.EqualTo(3));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(ClassifierLinearIterations: 0));

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierLinearIterations: 8), new CommonParameters(), [], []).Run();

        int found = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
    }

    /// <summary>
    /// Library fragments beyond the scored top N enter as their own features (DIA-NN's remaining-fragment correlations), so
    /// they add evidence without diluting the core scores: their co-elution with the profile, the share seen at the apex, and
    /// the co-elution weighted by library intensity. 0 unless asked for, and fewer than 0 is refused. A planted precursor's
    /// extra fragments co-elute.
    /// </summary>
    [Test]
    public void FragmentsBeyondTheTopNAreTheirOwnFeatures()
    {
        int[] extra = new[] { "ExtraCoElution", "ExtraMatchedFraction", "ExtraWeightedCoElution" }.Select(n => Array.IndexOf(DiaPrecursorMatch.FeatureNames, n)).ToArray();
        Assert.That(extra, Has.All.GreaterThanOrEqualTo(0));
        // The default since 2026-10-04: fragments 13-24. Two-seed means at a matched paired FDP of 1%, against the same build
        // without them: HF-X 53,859 -> 54,338 (+0.9%), PXD005573 1 h 39,009 -> 39,211 (+0.5%)
        Assert.That(new DiaLibrarySearchParameters().ExtraFragmentCount, Is.EqualTo(12));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(ExtraFragmentCount: -1));

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(int count) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(TopFragmentCount: 4, ExtraFragmentCount: count), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(0).Matches.SelectMany(m => extra.Select(i => m.Features[i])), Has.All.EqualTo(0));
        var planted = SearchWith(2).Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence)).ToList();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted.Average(m => m.Features[extra[0]]), Is.GreaterThan(0.5));
        Assert.That(planted.Average(m => m.Features[extra[1]]), Is.GreaterThan(0.5));
    }

    /// <summary>
    /// DIA-NN's co-elution per fragment is the best of its correlations at the full tolerance, 0.45x and 0.2x, so a noise peak
    /// at the wide tolerance cannot spoil a fragment whose real peak sits close. As a feature it is 0 unless asked for, and
    /// never below the full-tolerance co-elution of the same fragments.
    /// </summary>
    [Test]
    public void CoElutionCanTakeEachFragmentsBestTolerance()
    {
        int best = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "MaxToleranceCoElution");
        int plain = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "CoElution");
        Assert.That(best, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().MaxToleranceCoElution, Is.False);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(MaxToleranceCoElution: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[best]), Has.All.EqualTo(0));
        var on = SearchWith(true).Matches;
        Assert.That(on.Where(m => m.Features[plain] > 0).Select(m => m.Features[best] - m.Features[plain]), Has.All.GreaterThanOrEqualTo(-1e-12));
        Assert.That(on.Where(m => !m.IsDecoy && run.PlantedSequences.Contains(m.FullSequence)).Average(m => m.Features[best]), Is.GreaterThan(0.5));
    }

    /// <summary>
    /// DIA-NN's fragment rule: a fragment is scored only if it spans at least 3 residues and lies within 200-1800 m/z, so
    /// short ions shared by many peptides (y1, y2, b2) and the crowded low-m/z region do not count as evidence. On by default;
    /// a search still finds what was planted (the synthetic library's 2-residue fragments are left out).
    /// </summary>
    [Test]
    public void FragmentsCanBeLimitedAsDiaNnLimitsThem()
    {
        Omics.SpectralMatch.MslSpectralLibrary.MslFragmentIon Ion(double mz, int number) => new() { Mz = (float)mz, FragmentNumber = number, Intensity = 1 };
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(450, 3)), Is.True);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(450, 2)), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(199.9, 4)), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(1800.1, 9)), Is.False);
        // The limits can be set; DIA-NN's are the defaults
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(250, 3), minimumMz: 300, minimumResidues: 3), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(450, 3), minimumMz: 200, minimumResidues: 4), Is.False);
        Assert.That(new DiaLibrarySearchParameters().FragmentMinimumMz, Is.EqualTo(200));
        Assert.That(new DiaLibrarySearchParameters().FragmentMinimumResidues, Is.EqualTo(3));
        // A fragment charge limit can be set too (0, the default, sets none)
        var doubly = new Omics.SpectralMatch.MslSpectralLibrary.MslFragmentIon { Mz = 450, FragmentNumber = 5, Charge = 2, Intensity = 1 };
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(doubly), Is.True);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(doubly, maximumCharge: 1), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(doubly, maximumCharge: 2), Is.True);
        Assert.That(new DiaLibrarySearchParameters().FragmentMaximumCharge, Is.EqualTo(0));
        // Or relative to the precursor: a fragment cannot carry a precursor's whole charge unless its complement is neutral
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(doubly, precursorCharge: 2, chargeBelowPrecursor: true), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(doubly, precursorCharge: 3, chargeBelowPrecursor: true), Is.True);
        Assert.That(new DiaLibrarySearchParameters().FragmentChargeBelowPrecursor, Is.False);
        // Only terminal ions: FragmentNumber is an internal ion's start residue, and a diagnostic ion has no residues (MSL, 010)
        var internalIon = new Omics.SpectralMatch.MslSpectralLibrary.MslFragmentIon { Mz = 450, FragmentNumber = 4, SecondaryFragmentNumber = 9,
            ProductType = Omics.Fragmentation.ProductType.b, SecondaryProductType = Omics.Fragmentation.ProductType.y, Charge = 1, Intensity = 1 };
        var diagnostic = new Omics.SpectralMatch.MslSpectralLibrary.MslFragmentIon { Mz = 450, FragmentNumber = 5, ProductType = Omics.Fragmentation.ProductType.D, Charge = 1, Intensity = 1 };
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(internalIon), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(diagnostic), Is.False);
        // The default since 2026-10-05: two-seed means at a matched paired FDP of 1%, HF-X 54,777 -> 55,218 (+0.8%),
        // PXD005573 1 h 39,052 -> 40,558 (+3.9%)
        Assert.That(new DiaLibrarySearchParameters().DiaNnFragmentFilter, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(DiaNnFragmentFilter: true), new CommonParameters(), [], []).Run();
        int found = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
    }

    /// <summary>
    /// The MS1 tolerance (precursor traces, isotope envelope and MS1 mass error) can be set; 10 ppm by default. A tighter one
    /// keeps random MS1 peaks out of the traces. A search at 10 ppm still finds what was planted, and a precursor's MS1 mass
    /// error never exceeds the tolerance it was read with.
    /// </summary>
    [Test]
    public void TheMs1ToleranceCanBeSet()
    {
        // 5 ppm since 2026-10-06, with the calibrated MS1 offset applied: two-seed means at a matched paired FDP of 1%,
        // +1.6% HF-X, +1.3% PXD005573 against 10 ppm without it (10 since 2026-10-05, 20 before)
        Assert.That(new DiaLibrarySearchParameters().Ms1TolerancePpm, Is.EqualTo(5));
        int ms1Error = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1MassErrorPpm");

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1TolerancePpm: 10), new CommonParameters(), [], []).Run();

        Assert.That(results.Matches.Select(m => m.Features[ms1Error]), Has.All.LessThanOrEqualTo(10));
        int found = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
    }

    /// <summary>
    /// MS1 co-elution can also be read at half the MS1 tolerance, as DIA-NN scores MS1 at several tolerances: the precursor's
    /// M0 trace keeps only peaks within half the tolerance, correlated with the fragment profile. 0 unless asked for; a
    /// planted precursor's M0 sits at its library m/z, so its tight correlation stays high.
    /// </summary>
    [Test]
    public void Ms1CoElutionCanBeReadAtHalfTheTolerance()
    {
        int tight = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1TightCorrelation");
        Assert.That(tight, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().Ms1TightCorrelation, Is.False);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1TightCorrelation: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[tight]), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[tight]).Order().ToArray();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted[planted.Length / 2], Is.GreaterThan(0.5));
    }

    /// <summary>
    /// Co-elution can also be read on square-root traces, which damp the apex scans so the peak's flanks weigh more: the core
    /// fragments' mean correlation with the square-rooted profile. 0 unless asked for; a planted precursor's fragments share
    /// one elution shape, so it stays high.
    /// </summary>
    [Test]
    public void CoElutionCanBeReadOnSquareRootTraces()
    {
        int sqrt = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "SqrtCoElution");
        Assert.That(sqrt, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().SqrtCoElution, Is.False);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(SqrtCoElution: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[sqrt]), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[sqrt]).Order().ToArray();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted[planted.Length / 2], Is.GreaterThan(0.5));
    }

    /// <summary>
    /// The MS1 envelope read across the peak, not at one scan: M-1 to M3 summed over the co-elution window's MS1 scans,
    /// against the expected pattern with nothing at M-1. A clean envelope scores near 1; a peak one isotope below as large as
    /// M0 (our M0 is then probably another ion's M+1) pulls it well down.
    /// </summary>
    [Test]
    public void PeakEnvelopeCosinePenalisesAPeakOneIsotopeBelow()
    {
        double[] expected = DiaLibrarySearchEngine.ExpectedIsotopes(1500, 4);
        double[] clean = [0, .. expected.Select(e => 1000 * e)];
        double[] shadowed = [1000 * expected[0], .. expected.Select(e => 1000 * e)];

        Assert.That(DiaLibrarySearchEngine.PeakEnvelopeCosine(clean, expected), Is.EqualTo(1).Within(1e-9));
        Assert.That(DiaLibrarySearchEngine.PeakEnvelopeCosine(shadowed, expected), Is.LessThan(0.85));
        Assert.That(DiaLibrarySearchEngine.PeakEnvelopeCosine(new double[5], expected), Is.EqualTo(0));
    }

    /// <summary>
    /// The peak envelope as a feature (Ms1PeakEnvelope, on by default): 0 when switched off; on, a planted precursor's M0-M3 follow the
    /// expected pattern on every MS1 scan, so its summed envelope scores near 1.
    /// </summary>
    [Test]
    public void Ms1EnvelopeCanBeReadAcrossThePeak()
    {
        int envelope = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1PeakEnvelopeCosine");
        Assert.That(envelope, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().Ms1PeakEnvelope, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1PeakEnvelope: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[envelope]), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[envelope]).Order().ToArray();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted[planted.Length / 2], Is.GreaterThan(0.95));
    }

    /// <summary>
    /// The peak envelope again from tight peaks only (Ms1PeakEnvelopeTight): a second column reading only MS1 peaks within
    /// 0.6x the MS1 tolerance, so the network sees the envelope at two widths while lookups stay at the full one. 0 unless
    /// asked for; on, a planted precursor's isotopes sit at their m/z, so its tight envelope also scores near 1.
    /// </summary>
    [Test]
    public void Ms1EnvelopeCanAlsoBeReadFromTightPeaksOnly()
    {
        int tight = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1PeakEnvelopeTightCosine");
        Assert.That(tight, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().Ms1PeakEnvelopeTight, Is.False);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on, double ppmOffset = 0) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1PeakEnvelopeTight: on, Ms1Offset: ppmOffset == 0 ? null : Ms1OffsetModel.Constant(ppmOffset)), new CommonParameters(), [], []).Run();
        double MedianFor(DiaLibrarySearchResults results)
        {
            var values = results.Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[tight]).Order().ToArray();
            Assert.That(values, Is.Not.Empty);
            return values[values.Length / 2];
        }

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[tight]), Has.All.EqualTo(0));
        Assert.That(MedianFor(SearchWith(true)), Is.GreaterThan(0.95));
        Assert.That(MedianFor(SearchWith(true, ppmOffset: 4)), Is.LessThan(0.5), "peaks 4 ppm off are outside 3 ppm, though inside 5");
    }    /// <summary>
    /// The classifier network's hidden layers can be set (null keeps the rescorer's default, DIA-NN 2020's 25-20-15-10-5). A
    /// search with another architecture gives other scores, and a zero-unit layer is refused.
    /// </summary>
    [Test]
    public void TheClassifierNetworkLayersCanBeSet()
    {
        Assert.That(new DiaLibrarySearchParameters().ClassifierNetworkLayers, Is.Null);

        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(int[]? layers) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierNetworkLayers: layers, ClassifierNetworkMembers: 3), new CommonParameters(), [], []).Run();

        var standard = SearchWith(null);
        var wide = SearchWith([16, 8]);

        Assert.That(wide.Matches.Select(m => m.Score), Is.Not.EqualTo(standard.Matches.Select(m => m.Score)));
        Assert.Throws<ArgumentOutOfRangeException>(() => SearchWith([16, 0]), "a zero-unit layer is refused, not ignored");
    }

    /// <summary>
    /// The confident network sample can keep only targets that pass a q-value as positives (null, the default, keeps the top
    /// half-cap of targets whatever their q). With a training cap small enough for the confident sample to apply, the scores
    /// change.
    /// </summary>
    [Test]
    public void TheNetworksPositivesCanBeLimitedToPassingTargets()
    {
        Assert.That(new DiaLibrarySearchParameters().ClassifierNetworkPositiveQValue, Is.Null);

        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(double? cutoff) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierNetworkPositiveQValue: cutoff, MaxNetworkTrainingRows: 100, ClassifierNetworkMembers: 3),
            new CommonParameters(), [], []).Run();

        // Decoys are planted too, so few targets pass a strict cutoff here; a loose one bites
        Assert.That(SearchWith(0.25).Matches.Select(m => m.Score), Is.Not.EqualTo(SearchWith(null).Matches.Select(m => m.Score)));
    }

    /// <summary>
    /// MS1 evidence across the peak rather than at the apex: the M0 mass error averaged over the co-elution window, weighted by
    /// intensity, and the M+2 isotope trace's correlation with the fragment profile. Both 0 unless asked for; on, a planted
    /// precursor's M0 sits at its library m/z (small mass error) and its M+2 co-elutes.
    /// </summary>
    [Test]
    public void Ms1EvidenceCanBeReadAcrossThePeak()
    {
        int error = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1PeakMassErrorPpm");
        int m2 = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1Isotope2Correlation");
        Assert.That(new[] { error, m2 }, Has.All.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().Ms1PeakFeatures, Is.False);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1PeakFeatures: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.SelectMany(m => new[] { m.Features[error], m.Features[m2] }), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).ToList();
        double Median(IEnumerable<double> values) { var s = values.Order().ToArray(); return s[s.Length / 2]; }
        Assert.That(planted, Is.Not.Empty);
        Assert.That(Median(planted.Select(m => m.Features[error])), Is.LessThan(1));
        Assert.That(Median(planted.Select(m => m.Features[m2])), Is.GreaterThan(0.5));
    }

    /// <summary>
    /// A run's MS1 offset can be applied: every MS1 lookup moves to the library m/z shifted by the offset at that retention
    /// time, and the MS1 mass error feature is measured from the shifted m/z. Each match also keeps its raw signed apex MS1
    /// error (observed minus library m/z), from which an offset is fitted. On planted precursors, whose M0 sits at the library
    /// m/z, the raw error stays near 0 whatever the offset, and a +3 ppm offset shows up as about 3 ppm of feature error.
    /// </summary>
    [Test]
    public void AnMs1OffsetShiftsTheMs1Lookups()
    {
        Assert.That(new DiaLibrarySearchParameters().Ms1Offset, Is.Null);
        int error = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1MassErrorPpm");
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        List<DiaPrecursorMatch> PlantedWith(Ms1OffsetModel? offset) => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library,
            EndpointMap, new DiaLibrarySearchParameters(Ms1Offset: offset), new CommonParameters(), [], []).Run())
            .Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence) && double.IsFinite(m.Ms1ApexErrorPpm)).ToList();
        double Median(IEnumerable<double> values) { var s = values.Order().ToArray(); return s[s.Length / 2]; }

        var none = PlantedWith(null);
        var shifted = PlantedWith(Ms1OffsetModel.Constant(3));

        Assert.That(none, Is.Not.Empty);
        Assert.That(Median(none.Select(m => Math.Abs(m.Ms1ApexErrorPpm))), Is.LessThan(1));
        Assert.That(Median(shifted.Select(m => Math.Abs(m.Ms1ApexErrorPpm))), Is.LessThan(1), "the raw error does not depend on the offset");
        Assert.That(Median(none.Select(m => m.Features[error])), Is.LessThan(1));
        Assert.That(Median(shifted.Select(m => m.Features[error])), Is.EqualTo(3).Within(0.6));
    }

    /// <summary>
    /// The MS1 tolerance can follow the run's own mass accuracy: with Ms1ToleranceSpreadMultiple, it is that multiple of the
    /// calibrated offset's residual spread. Without the option, or without a measured spread, it stays Ms1TolerancePpm.
    /// </summary>
    [Test]
    public void TheMs1ToleranceCanFollowTheRunsMassAccuracy()
    {
        var offset = Ms1OffsetModel.Fit(Enumerable.Range(0, 400).Select(i => (i / 10.0, 2 + (i % 2 == 0 ? 1.0 : -1.0))));
        Assert.That(new DiaLibrarySearchParameters().Ms1ToleranceSpreadMultiple, Is.Null);

        Assert.That(new DiaLibrarySearchParameters(Ms1Offset: offset).EffectiveMs1TolerancePpm, Is.EqualTo(5));
        Assert.That(new DiaLibrarySearchParameters(Ms1Offset: offset, Ms1ToleranceSpreadMultiple: 3).EffectiveMs1TolerancePpm,
            Is.EqualTo(3 * 1.4826).Within(0.1));
        Assert.That(new DiaLibrarySearchParameters(Ms1ToleranceSpreadMultiple: 3).EffectiveMs1TolerancePpm, Is.EqualTo(5), "no offset fitted");
        Assert.That(new DiaLibrarySearchParameters(Ms1Offset: Ms1OffsetModel.Constant(2), Ms1ToleranceSpreadMultiple: 3).EffectiveMs1TolerancePpm,
            Is.EqualTo(5), "no spread measured");
    }

    /// <summary>
    /// A cheap gate before the full features, as DIA-NN's: a scan can be a candidate apex only if at least this many of the
    /// six most intense library fragments are seen there. Noise-only candidates are never scored, so the search scores
    /// fewer of them and still finds what was planted.
    /// </summary>
    [Test]
    public void CandidateApexesCanBeGatedOnTheirTopFragments()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 2000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(int gate) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(MinimumApexFragments: gate, DiaNnPeakFinding: false), new CommonParameters(), [], []).Run();

        var open = SearchWith(0);
        var gated = SearchWith(4);

        int Candidates(DiaLibrarySearchResults r) => r.Matches.Count + r.LosingCandidates.Values.Sum(v => v.Count);
        Assert.That(Candidates(gated), Is.LessThan(Candidates(open)));
        int found = gated.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(MinimumApexFragments: 7));
    }

    /// <summary>
    /// One more candidate by a different rule: the apex of the best single fragment's smoothed trace over the whole window,
    /// when no candidate is already there. In about 2,300 of the DIA-NN IDs we miss per file, DIA-NN's peak was never among
    /// our candidates, and more candidates by the same apex score did not help.
    /// </summary>
    [Test]
    public void TheBestFragmentsApexCanBeAnExtraCandidate()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 2000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool extra) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(FragmentApexCandidate: extra), new CommonParameters(), [], []).Run();

        var without = SearchWith(false);
        var with = SearchWith(true);

        int Candidates(DiaLibrarySearchResults r) => r.Matches.Count + r.LosingCandidates.Values.Sum(v => v.Count);
        Assert.That(Candidates(with), Is.GreaterThan(Candidates(without)));
        Assert.That(new DiaLibrarySearchParameters().FragmentApexCandidate, Is.False);
    }

    /// <summary>
    /// The classifier can normalise each fold on a sample of its training groups rather than scoring every training row.
    /// The setting reaches the classifier: the same precursors are scored, and their scores change.
    /// </summary>
    [Test]
    public void TheClassifierCanNormaliseOnASampleOfGroups()
    {
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(int? groups) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierNormalizationGroups: groups), new CommonParameters(), [], []).Run();

        var all = SearchWith(null);
        var sampled = SearchWith(20);

        Assert.That(sampled.Matches.Select(m => m.PrecursorIndex), Is.EquivalentTo(all.Matches.Select(m => m.PrecursorIndex)));
        Assert.That(sampled.Matches.Select(m => m.Score), Is.Not.EqualTo(all.Matches.Select(m => m.Score)));
        // The default since 2026-10-04: each fold's network otherwise scores every row three times over; HF-X 54,017 -> 54,126,
        // PXD005573 1 h 39,282 -> 39,193 at a matched paired FDP of 1%, classifier 3:26 -> 2:51 on HF-X
        Assert.That(new DiaLibrarySearchParameters().ClassifierNormalizationGroups, Is.EqualTo(300_000));
    }

    /// <summary>
    /// A fragment's trace point can be the most intense peak within tolerance, as DIA-NN reads it (diann.cpp 1.8, level()),
    /// rather than the nearest: a faint noise peak closer to the library m/z otherwise stands in for the real fragment.
    /// With none in tolerance there is no peak.
    /// </summary>
    [Test]
    public void AFragmentCanBeReadAsTheMostIntensePeakInTolerance()
    {
        var spectrum = new MassSpectrometry.MzSpectrum([499.990, 500.001, 500.006, 500.03], [5e5, 1e3, 2e5, 9e6], false);
        var tolerance = new MzLibUtil.PpmTolerance(20);

        Assert.That(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, 500.0, tolerance, mostIntense: false), Is.EqualTo(1));
        Assert.That(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, 500.0, tolerance, mostIntense: true), Is.EqualTo(0));
        Assert.That(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, 501.0, tolerance, mostIntense: true), Is.EqualTo(-1));
        Assert.That(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, 501.0, tolerance, mostIntense: false), Is.EqualTo(-1));
        // The default since 2026-10-04: +0.2% HF-X, +1.4% PXD005573 at a matched paired FDP of 1%
        Assert.That(new DiaLibrarySearchParameters().MostIntenseFragmentPeak, Is.True);
    }

    /// <summary>
    /// DIA-NN 1.8's candidate peaks: a scan is scored by its best fragment's summed correlation to the other top fragments.
    /// It needs at least two fragments, a score of 0.5, and the smoothed best fragment at a local maximum. Every peak within
    /// 2.0 of the best is kept. A five-fragment peak (score near 4) keeps a four-fragment one (near 3) but not a
    /// two-fragment one (near 1). A one-fragment spike is never a candidate, however intense.
    /// </summary>
    [Test]
    public void CandidatePeaksCanBeFoundAsDiaNnFindsThem()
    {
        static double[][] Traces(params (int Apex, int Fragments, double Height)[] peaks)
        {
            var traces = Enumerable.Range(0, 6).Select(_ => new double[40]).ToArray();
            foreach (var (apex, fragments, height) in peaks)
                for (int f = 0; f < fragments; f++)
                    for (int s = Math.Max(0, apex - 5); s <= Math.Min(39, apex + 5); s++)
                        traces[f][s] += height * (f + 1) * Math.Exp(-0.5 * (s - apex) * (s - apex) / 2.25);
            return traces;
        }
        var withFourFragmentPeak = Traces((10, 5, 1e4), (32, 4, 1e4));
        withFourFragmentPeak[0][25] = 1e8;
        var withTwoFragmentPeak = Traces((10, 5, 1e4), (32, 2, 1e4));
        withTwoFragmentPeak[0][25] = 1e8;

        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(withFourFragmentPeak, 3, 10), Is.EqualTo(new[] { 10, 32 }));
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(withTwoFragmentPeak, 3, 10), Is.EqualTo(new[] { 10 }));
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(withFourFragmentPeak, 3, 1), Is.EqualTo(new[] { 10 }));
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(Traces(), 3, 10), Is.Empty);

        // DIA-NN adds the MS1 trace's correlation with the best fragment to the score: two fragments that do not co-elute
        // score 0 and are no candidate, until an MS1 peak shaped like one of them lifts the score to 1
        var gaussian = Enumerable.Range(0, 40).Select(s => Math.Abs(s - 20) <= 5 ? 1e4 * Math.Exp(-0.5 * (s - 20) * (s - 20) / 2.25) : 0).ToArray();
        var inverted = gaussian.Select(v => 2e4 - v).ToArray();
        var unrelated = new[] { gaussian, inverted, new double[40], new double[40], new double[40], new double[40] };
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(unrelated, 3, 10), Is.Empty);
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(unrelated, 3, 10, ms1: gaussian), Is.EqualTo(new[] { 20 }));
        Assert.That(new DiaLibrarySearchParameters().DiaNnPeakFindingMs1, Is.False);
        // The default since 2026-10-04, with 6 candidates: neutral at 3, +1.3% HF-X and +2.8% PXD005573 at 6 (matched paired FDP 1%)
        Assert.That(new DiaLibrarySearchParameters().DiaNnPeakFinding, Is.True);
        // 10 since 2026-10-05, with the fragment rule: two-seed means +1.2% HF-X, +0.6% PXD005573 (matched paired FDP 1%)
        Assert.That(new DiaLibrarySearchParameters().MaxApexCandidates, Is.EqualTo(10));
    }

    /// <summary>
    /// A fragment's peak can be looked up from a nearby index (the fragment it shadows, one isotope up) instead of by binary
    /// search: the same peak every time, nearest or most intense, from any starting index, including targets off either end.
    /// </summary>
    [Test]
    public void APeakLookedUpFromANearbyIndexIsTheSamePeak()
    {
        var random = new Random(7);
        var tolerance = new MzLibUtil.PpmTolerance(20);
        for (int trial = 0; trial < 200; trial++)
        {
            double[] mz = Enumerable.Range(0, 1 + random.Next(300)).Select(_ => 100 + 1900 * random.NextDouble()).Distinct().OrderBy(x => x).ToArray();
            double[] intensity = mz.Select(_ => (double)random.Next(1, 50)).ToArray();
            var spectrum = new MassSpectrometry.MzSpectrum(mz, intensity, false);
            for (int k = 0; k < 20; k++)
            {
                double target = random.Next(3) == 0 ? mz[random.Next(mz.Length)] * (1 + (random.NextDouble() - 0.5) * 4e-5) : 50 + 2000 * random.NextDouble();
                int hint = random.Next(mz.Length);
                foreach (bool mostIntense in new[] { false, true })
                    Assert.That(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, target, tolerance, mostIntense, hint),
                        Is.EqualTo(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, target, tolerance, mostIntense)), $"trial {trial}, target {target}, hint {hint}");
            }
        }
    }

    /// <summary>
    /// A spectrum's m/z bin table: for each 2-Th bin, the first peak at or above the bin's lower edge. It is the start for the
    /// hinted lookup, so a fragment's peak is found without a binary search; any m/z, even off either end, gets a valid start.
    /// </summary>
    [Test]
    public void ASpectrumsBinTableStartsEachBinAtItsFirstPeak()
    {
        var random = new Random(11);
        double[] mz = Enumerable.Range(0, 500).Select(_ => 150 + 1700 * random.NextDouble()).Distinct().OrderBy(x => x).ToArray();
        var spectrum = new MassSpectrometry.MzSpectrum(mz, mz.Select(_ => 1.0).ToArray(), false);
        var bins = DiaLibrarySearchEngine.PeakBins.Of(spectrum);

        for (double edge = 150; edge < 1850; edge += DiaLibrarySearchEngine.PeakBins.Width)
        {
            int start = bins.Hint(edge + 1e-9);
            int first = Array.FindIndex(mz, x => x >= Math.Floor((edge + 1e-9) / DiaLibrarySearchEngine.PeakBins.Width) * DiaLibrarySearchEngine.PeakBins.Width);
            Assert.That(start, Is.EqualTo(first < 0 ? mz.Length - 1 : first), $"bin at {edge}");
        }
        foreach (double target in new[] { 0.0, 100.0, 5000.0 })
            Assert.That(bins.Hint(target), Is.InRange(0, mz.Length - 1));
        Assert.That(DiaLibrarySearchEngine.PeakBins.Of(new MassSpectrometry.MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false)).Hint(500), Is.EqualTo(0));
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
        Assert.That(Median(chance.Select(m => m.Features[ppm])), Is.EqualTo(new DiaLibrarySearchParameters().Ms1TolerancePpm).Within(1e-9), "no M0 found counts as the full tolerance");
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
