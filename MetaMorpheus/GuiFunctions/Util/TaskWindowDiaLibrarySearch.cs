using System.Collections.Generic;
using System.Globalization;
using TaskLayer;

namespace GuiFunctions.Util;

/// <summary>
/// What the DIA library search window shows and saves for <see cref="DiaLibrarySearchTaskParameters"/>: a text box for each
/// number, the precursor charges as a list ("2, 3"), and the decoy kind as a choice. The window loads them with
/// <see cref="Load"/>, checks what the user typed with <see cref="Error"/>, and saves it with <see cref="ForSave"/>.
/// </summary>
/// <remarks>
/// Only the form of each value is checked here. What a prediction model can take (its charges, its NCE) is the task's to
/// refuse, with a message, before it asks Koina: the window does not restate the model's limits.
/// </remarks>
public static class TaskWindowDiaLibrarySearch
{
    /// <summary>The text the window shows for each parameter, and the decoy kind it has chosen.</summary>
    public sealed record Fields(string FragmentTolerancePpm, string Ms1TolerancePpm, string ClassifierSeed,
        string LibraryPrecursorCharges, string LibraryNce, DiaLibraryDecoys LibraryDecoys);

    /// <summary>What the window shows for the loaded task's parameters.</summary>
    public static Fields Load(DiaLibrarySearchTaskParameters loaded) =>
        new(loaded.FragmentTolerancePpm.ToString(CultureInfo.InvariantCulture),
            loaded.Ms1TolerancePpm.ToString(CultureInfo.InvariantCulture),
            loaded.ClassifierSeed.ToString(CultureInfo.InvariantCulture),
            string.Join(", ", loaded.LibraryPrecursorCharges),
            loaded.LibraryNce.ToString(CultureInfo.InvariantCulture),
            loaded.LibraryDecoys);

    /// <summary>Why <paramref name="shown"/> cannot be saved, one line per field that is wrong, or null when it can.</summary>
    public static string Error(Fields shown)
    {
        var errors = new List<string>();
        if (!(Number(shown.FragmentTolerancePpm) is > 0))
            errors.Add($"The fragment tolerance must be a number of ppm above 0. You entered \"{shown.FragmentTolerancePpm}\".");
        if (!(Number(shown.Ms1TolerancePpm) is > 0))
            errors.Add($"The MS1 tolerance must be a number of ppm above 0. You entered \"{shown.Ms1TolerancePpm}\".");
        if (Integer(shown.ClassifierSeed) is null)
            errors.Add($"The classifier seed must be a whole number. You entered \"{shown.ClassifierSeed}\".");
        if (Charges(shown.LibraryPrecursorCharges) is null)
            errors.Add("The library precursor charges must be one or more whole numbers above 0, such as \"2, 3\". "
                + $"You entered \"{shown.LibraryPrecursorCharges}\".");
        if (!(Integer(shown.LibraryNce) is > 0))
            errors.Add($"The library NCE must be a whole number above 0. You entered \"{shown.LibraryNce}\".");
        return errors.Count == 0 ? null : string.Join("\n", errors);
    }

    /// <summary>The parameters the window saves: what <paramref name="shown"/> holds, which <see cref="Error"/> has passed.</summary>
    public static DiaLibrarySearchTaskParameters ForSave(Fields shown) => new()
    {
        FragmentTolerancePpm = Number(shown.FragmentTolerancePpm)!.Value,
        Ms1TolerancePpm = Number(shown.Ms1TolerancePpm)!.Value,
        ClassifierSeed = Integer(shown.ClassifierSeed)!.Value,
        LibraryPrecursorCharges = Charges(shown.LibraryPrecursorCharges)!,
        LibraryNce = Integer(shown.LibraryNce)!.Value,
        LibraryDecoys = shown.LibraryDecoys,
    };

    private static double? Number(string text) =>
        double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out double value) && double.IsFinite(value) ? value : null;

    private static int? Integer(string text) =>
        int.TryParse(text, NumberStyles.Integer, CultureInfo.InvariantCulture, out int value) ? value : null;

    /// <summary>The charges in a list separated by commas or spaces, each once, in the order given; null unless every item is a whole number above 0.</summary>
    private static List<int> Charges(string text)
    {
        var charges = new List<int>();
        foreach (string item in text.Split(new[] { ',', ' ' }, System.StringSplitOptions.RemoveEmptyEntries))
        {
            if (!(Integer(item) is int charge and > 0))
                return null;
            if (!charges.Contains(charge))
                charges.Add(charge);
        }
        return charges.Count > 0 ? charges : null;
    }
}
