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
        new("", "", "", "", "", DiaLibraryDecoys.Permuted); // red: not written yet

    /// <summary>Why <paramref name="shown"/> cannot be saved, one line per field that is wrong, or null when it can.</summary>
    public static string Error(Fields shown) => null; // red: not written yet

    /// <summary>The parameters the window saves: what <paramref name="shown"/> holds, which <see cref="Error"/> has passed.</summary>
    public static DiaLibrarySearchTaskParameters ForSave(Fields shown) => new(); // red: not written yet
}
