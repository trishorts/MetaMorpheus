#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// What a run's confident targets say of the library's iRT, by sequence composition. A predicted iRT is fixed for some other column, and on both benchmark runs what is left after the run's RT->iRT
/// curve has a shape: K, R and H elute 1.5 to 3 iRT earlier than predicted. Calibration fits this on its anchors, and the search expects each precursor where the fit says, in the window it reads as
/// well as in the error it scores. Fitted on 1,000 anchors, the median error of the targets a search scores too low goes from 2.47 to 1.59 iRT on PXD022589 HF-X, and from 2.22 to 1.25 on PXD005573 1 h.
/// </summary>
/// <remarks>
/// Only targets fit it, but it goes by composition alone, so a permuted decoy moves as its target does.
/// </remarks>
public sealed class DiaIrtCompositionModel
{
    /// <summary>With fewer anchors the fit moves nothing: it has as many terms as a few dozen, and calibration usually has 2,000.</summary>
    public const int MinimumAnchors = 200;

    /// <summary>How far a search looks past the run's iRT range, since a precursor's expected iRT can lie this far from its library value.</summary>
    public const double QueryMargin = 20;

    private const string Acids = "ACDEFGHIKLMNPQRSTVWY";
    private const string Basic = "KRH";
    private const string Acidic = "DE";
    private const double Penalty = 1;

    /// <summary>Anchors further than this many standard deviations from the first fit do not take part in the second.</summary>
    private const double OutlierSds = 3;

    private static readonly int TermCount = Acids.Length + 10;

    private readonly double[]? _weights;

    private DiaIrtCompositionModel(double[]? weights) => _weights = weights;

    /// <summary>A confident target.</summary>
    /// <param name="BaseSequence">Its plain sequence.</param>
    /// <param name="FullSequence">Its sequence with modifications.</param>
    /// <param name="Charge">Its charge.</param>
    /// <param name="Residual">Its apex iRT on the run's curve, less its library iRT.</param>
    public sealed record Anchor(string BaseSequence, string FullSequence, int Charge, double Residual);

    /// <summary>How far from its library iRT this run's curve should find a precursor.</summary>
    public double Offset(string baseSequence, string fullSequence, int charge)
    {
        if (_weights is null)
            return 0;
        var terms = new double[TermCount];
        FillTerms(baseSequence, fullSequence, charge, terms);
        double offset = 0;
        for (int i = 0; i < terms.Length; i++)
            offset += _weights[i] * terms[i];
        return offset;
    }

    /// <summary>A ridge fit of the residuals on composition, length, charge and the first two residues, then again without the anchors the first fit could not explain.</summary>
    public static DiaIrtCompositionModel Fit(IReadOnlyCollection<Anchor> anchors)
    {
        ArgumentNullException.ThrowIfNull(anchors);
        if (anchors.Count < MinimumAnchors)
            return new DiaIrtCompositionModel(null);
        var rows = anchors.Select(a =>
        {
            var terms = new double[TermCount];
            FillTerms(a.BaseSequence, a.FullSequence, a.Charge, terms);
            return (Terms: terms, a.Residual);
        }).ToList();
        double[] weights = Solve(rows);
        double[] left = rows.Select(r => r.Residual - Dot(weights, r.Terms)).ToArray();
        double mean = left.Average();
        double sd = Math.Sqrt(left.Sum(v => (v - mean) * (v - mean)) / left.Length);
        var kept = rows.Where((_, i) => Math.Abs(left[i] - mean) <= OutlierSds * sd).ToList();
        return new DiaIrtCompositionModel(kept.Count >= MinimumAnchors ? Solve(kept) : weights);
    }

    private static double Dot(double[] a, double[] b)
    {
        double sum = 0;
        for (int i = 0; i < a.Length; i++)
            sum += a[i] * b[i];
        return sum;
    }

    private static double[] Solve(List<(double[] Terms, double Residual)> rows)
    {
        var xtx = new double[TermCount, TermCount];
        var xty = new double[TermCount];
        foreach (var (terms, residual) in rows)
            for (int i = 0; i < TermCount; i++)
            {
                xty[i] += terms[i] * residual;
                for (int j = 0; j < TermCount; j++)
                    xtx[i, j] += terms[i] * terms[j];
            }
        for (int i = 0; i < TermCount; i++)
            xtx[i, i] += Penalty;
        return Matrix<double>.Build.DenseOfArray(xtx).Solve(Vector<double>.Build.DenseOfArray(xty)).ToArray();
    }

    private static void FillTerms(string baseSequence, string fullSequence, int charge, double[] terms)
    {
        foreach (char residue in baseSequence)
        {
            int index = Acids.IndexOf(residue);
            if (index >= 0)
                terms[index]++;
        }
        int next = Acids.Length;
        double length = Math.Max(1, baseSequence.Length) / 10.0;
        terms[next++] = length;
        terms[next++] = length * length;
        terms[next++] = 1 / length;
        terms[next++] = charge == 2 ? 1 : 0;
        terms[next++] = charge == 3 ? 1 : 0;
        terms[next++] = baseSequence.Length > 0 && Acidic.Contains(baseSequence[0]) ? 1 : 0;
        terms[next++] = baseSequence.Length > 0 && Basic.Contains(baseSequence[0]) ? 1 : 0;
        terms[next++] = baseSequence.Length > 1 && Basic.Contains(baseSequence[1]) ? 1 : 0;
        // An oxidised methionine elutes far from a plain one, so it has a term of its own
        for (int at = fullSequence.IndexOf("Oxidation", StringComparison.Ordinal); at >= 0; at = fullSequence.IndexOf("Oxidation", at + 1, StringComparison.Ordinal))
            terms[next]++;
        next++;
        // The intercept
        terms[next] = 1;
    }
}
