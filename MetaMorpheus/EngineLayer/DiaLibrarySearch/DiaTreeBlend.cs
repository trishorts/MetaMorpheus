#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.ML;
using Microsoft.ML.Data;
using Microsoft.ML.Trainers.FastTree;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// Blends the classifier's order with that of gradient-boosted trees trained on the same features. Trees alone rank worse than
/// the network, but the two make different errors, so the mean of their ranks can beat either.
/// </summary>
/// <remarks>
/// Each row is scored by trees that never saw its group: groups are split into folds, and each fold's trees learn target
/// against decoy from the other folds' top rows by score, whatever their label, as the network's confident sample does. The
/// blend only decides the order. The values stay the classifier's own, so later stages see the same score scale.
/// </remarks>
public static class DiaTreeBlend
{
    private const int Folds = 3;

    /// <summary>The fewest targets, and the fewest decoys, a fold's trees may learn from. With fewer, the scores stay as they are.</summary>
    private const int MinimumRowsPerClass = 50;

    /// <summary>One row for ML.NET, which reads its public properties.</summary>
    public sealed class Row
    {
        public bool Label { get; set; }
        public float[] Features { get; set; } = [];
    }

    /// <param name="groupKeys">Rows sharing a key, such as a peptide's charge states, stay in one fold.</param>
    /// <param name="maxTrainingRows">The most rows each fold's trees learn from, taken from the top by score.</param>
    /// <returns>The same score values as <paramref name="scores"/>, reordered, or a copy of it when there is too little to learn from.</returns>
    /// <exception cref="ArgumentException">The lists differ in length.</exception>
    public static double[] Blend(IReadOnlyList<double[]> features, IReadOnlyList<double> scores, IReadOnlyList<bool> isDecoy,
        IReadOnlyList<string> groupKeys, int maxTrainingRows, int seed, int threads = 1)
    {
        int n = scores.Count;
        if (features.Count != n || isDecoy.Count != n || groupKeys.Count != n)
            throw new ArgumentException("Features, scores, decoy flags and group keys must have one entry per row.");

        var foldOfGroup = groupKeys.Distinct().Order(StringComparer.Ordinal).Select((key, index) => (key, index))
            .ToDictionary(group => group.key, group => group.index % Folds, StringComparer.Ordinal);
        int[] fold = groupKeys.Select(key => foldOfGroup[key]).ToArray();
        int[] byScore = Enumerable.Range(0, n).OrderByDescending(i => scores[i]).ThenBy(i => i).ToArray();

        var folds = new (int[] Train, int[] Test)[Folds];
        for (int f = 0; f < Folds; f++)
        {
            int current = f;
            int[] train = byScore.Where(i => fold[i] != current).Take(maxTrainingRows).ToArray();
            int[] test = Enumerable.Range(0, n).Where(i => fold[i] == current).ToArray();
            int decoys = train.Count(i => isDecoy[i]);
            if (decoys < MinimumRowsPerClass || train.Length - decoys < MinimumRowsPerClass || !test.Any(i => isDecoy[i]))
                return scores.ToArray();
            folds[f] = (train, test);
        }

        var treeScores = new double[n];
        foreach (var (train, test) in folds)
        {
            double[] margins = TrainAndScore(features, isDecoy, train, test, seed, threads);
            double[] decoyMargins = margins.Where((_, t) => isDecoy[test[t]]).Order().ToArray();
            double median = decoyMargins[decoyMargins.Length / 2];
            double spread = decoyMargins[(int)(0.99 * (decoyMargins.Length - 1))] - median;
            // Each fold has its own trees, so put their outputs on one scale before ranking them together
            for (int t = 0; t < test.Length; t++)
                treeScores[test[t]] = (margins[t] - median) / (spread > 0 ? spread : 1);
        }

        double[] scoreRanks = Ranks(scores);
        double[] treeRanks = Ranks(treeScores);
        int[] blendedOrder = Enumerable.Range(0, n).OrderByDescending(i => scoreRanks[i] + treeRanks[i]).ThenByDescending(i => scores[i]).ThenBy(i => i).ToArray();
        var blended = new double[n];
        for (int rank = 0; rank < n; rank++)
            blended[blendedOrder[rank]] = scores[byScore[rank]];
        return blended;
    }

    private static double[] TrainAndScore(IReadOnlyList<double[]> features, IReadOnlyList<bool> isDecoy, int[] train, int[] test, int seed, int threads)
    {
        var context = new MLContext(seed);
        var schema = SchemaDefinition.Create(typeof(Row));
        schema[nameof(Row.Features)].ColumnType = new VectorDataViewType(NumberDataViewType.Single, features[0].Length);
        IEnumerable<Row> Rows(int[] indices) => indices.Select(i => new Row { Label = !isDecoy[i], Features = features[i].Select(v => (float)v).ToArray() });

        var model = context.BinaryClassification.Trainers.FastTree(new FastTreeBinaryTrainer.Options
        {
            NumberOfThreads = threads,
            NumberOfTrees = 400,
            NumberOfLeaves = 128,
            MinimumExampleCountPerLeaf = 10,
            LearningRate = 0.08,
            FeatureFraction = 0.8,
            BaggingSize = 1,
            BaggingExampleFraction = 0.8,
            LabelColumnName = nameof(Row.Label),
            FeatureColumnName = nameof(Row.Features),
            Seed = seed,
            FeatureSelectionSeed = seed,
        }).Fit(context.Data.LoadFromEnumerable(Rows(train), schema));

        return model.Transform(context.Data.LoadFromEnumerable(Rows(test), schema))
            .GetColumn<float>("Score").Select(margin => (double)margin).ToArray();
    }

    /// <summary>Each value's rank from the lowest, ties sharing their mean rank.</summary>
    private static double[] Ranks(IReadOnlyList<double> values)
    {
        int[] order = Enumerable.Range(0, values.Count).OrderBy(i => values[i]).ToArray();
        var ranks = new double[values.Count];
        for (int start = 0; start < order.Length;)
        {
            int end = start;
            while (end + 1 < order.Length && values[order[end + 1]] == values[order[start]])
                end++;
            for (int k = start; k <= end; k++)
                ranks[order[k]] = (start + end) / 2.0;
            start = end + 1;
        }
        return ranks;
    }
}
