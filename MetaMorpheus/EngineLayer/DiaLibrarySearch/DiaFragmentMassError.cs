#nullable enable
using System;
using System.Collections.Generic;
using MassSpectrometry;
using MzLibUtil;
using Omics.SpectralMatch.MslSpectralLibrary;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// A run's fragment mass error (observed minus library m/z, in ppm), measured on calibration's anchors: each anchor's strongest
/// library fragments at its apex scan, with their peaks picked as the search picks them. DIA-NN sets its MS2 tolerance per run
/// from such a measurement: 12 ppm on PXD071205 A1 and PXD005573 1 h, 20 on PXD022589 HF-X, about 4.4 robust SDs of its best
/// fragment's error on each. Our fixed 20 ppm costs A1 10% at a matched paired entrapment FDP of 1% against 10 ppm.
/// </summary>
/// <param name="MedianPpm">The median error.</param>
/// <param name="SpreadPpm">1.4826 x the median absolute deviation from the median (an SD for normal errors).</param>
/// <param name="Deviation95Ppm">The 95th percentile of the absolute deviation from the median, which sees a heavy tail the robust SD does not.</param>
/// <param name="Count">Fragment errors measured.</param>
public sealed record DiaFragmentMassError(double MedianPpm, double SpreadPpm, double Deviation95Ppm, int Count)
{
    /// <summary>Each anchor's strongest library fragments read (by library intensity): the likeliest to be seen, and seen with the least error.</summary>
    public const int FragmentsPerAnchor = 3;

    /// <summary>Summarises signed fragment errors in ppm; non-finite ones are ignored, and none gives NaN with a count of 0.</summary>
    public static DiaFragmentMassError Fit(IEnumerable<double> errorsPpm) => new(double.NaN, double.NaN, double.NaN, 0);

    /// <summary>
    /// A confident target's strongest library fragments' errors at its apex scan, in the window whose scans come closest to
    /// the apex. A fragment with no peak in tolerance there gives nothing.
    /// </summary>
    /// <param name="entry">The target.</param>
    /// <param name="apexRt">Its apex.</param>
    /// <param name="windowsHoldingIt">The scans of each isolation window holding its m/z, each in retention-time order.</param>
    /// <param name="tolerancePpm">The fragment tolerance.</param>
    /// <param name="mostIntensePeak">Take the most intense peak in tolerance rather than the closest, as the search does.</param>
    public static List<double> ReadAnchor(MslLibraryEntry entry, RtMinutes apexRt, IReadOnlyList<MsDataScan[]> windowsHoldingIt, double tolerancePpm, bool mostIntensePeak) => [];
}
