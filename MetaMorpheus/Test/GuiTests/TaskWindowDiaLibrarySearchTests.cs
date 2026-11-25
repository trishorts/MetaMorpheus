using GuiFunctions.Util;
using NUnit.Framework;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading;
using TaskLayer;

namespace Test.GuiTests
{
    /// <summary>
    /// What the DIA library search window shows and saves for <see cref="DiaLibrarySearchTaskParameters"/>
    /// (<see cref="TaskWindowDiaLibrarySearch"/>). The WPF window cannot be referenced from this test project, so the window
    /// calls only these, and these are what is tested.
    /// </summary>
    [TestFixture]
    public static class TaskWindowDiaLibrarySearchTests
    {
        private static DiaLibrarySearchTaskParameters Custom() => new()
        {
            FragmentTolerancePpm = 10.5,
            Ms1TolerancePpm = 7.25,
            ClassifierSeed = 42,
            LibraryPrecursorCharges = [2, 3, 4],
            LibraryNce = 30,
            LibraryDecoys = DiaLibraryDecoys.Reversed,
        };

        private static TaskWindowDiaLibrarySearch.Fields Good() => TaskWindowDiaLibrarySearch.Load(new DiaLibrarySearchTaskParameters());

        [Test]
        public static void TheWindowShowsEachParameterOfTheLoadedTask()
        {
            var shown = TaskWindowDiaLibrarySearch.Load(Custom());

            Assert.That(shown, Is.EqualTo(new TaskWindowDiaLibrarySearch.Fields("10.5", "7.25", "42", "2, 3, 4", "30", DiaLibraryDecoys.Reversed)));
        }

        [Test]
        public static void SavingWhatTheWindowShowsGivesTheLoadedParametersBack()
        {
            foreach (var loaded in new[] { new DiaLibrarySearchTaskParameters(), Custom() })
            {
                var shown = TaskWindowDiaLibrarySearch.Load(loaded);
                Assert.That(TaskWindowDiaLibrarySearch.Error(shown), Is.Null);
                var saved = TaskWindowDiaLibrarySearch.ForSave(shown);
                Assert.That(saved.FragmentTolerancePpm, Is.EqualTo(loaded.FragmentTolerancePpm));
                Assert.That(saved.Ms1TolerancePpm, Is.EqualTo(loaded.Ms1TolerancePpm));
                Assert.That(saved.ClassifierSeed, Is.EqualTo(loaded.ClassifierSeed));
                Assert.That(saved.LibraryPrecursorCharges, Is.EqualTo(loaded.LibraryPrecursorCharges));
                Assert.That(saved.LibraryNce, Is.EqualTo(loaded.LibraryNce));
                Assert.That(saved.LibraryDecoys, Is.EqualTo(loaded.LibraryDecoys));
            }
        }

        /// <summary>
        /// A parameter added to the task later must get a field in the window too, or every save from the GUI would reset it to
        /// its default. Each writable property of <see cref="Custom"/> differs from the default, and must survive the window.
        /// </summary>
        [Test]
        public static void EveryParameterTheTaskHasIsCarriedThroughTheWindow()
        {
            var loaded = Custom();
            var saved = TaskWindowDiaLibrarySearch.ForSave(TaskWindowDiaLibrarySearch.Load(loaded));
            var defaults = new DiaLibrarySearchTaskParameters();

            foreach (var property in typeof(DiaLibrarySearchTaskParameters).GetProperties().Where(p => p.CanWrite))
            {
                object custom = property.GetValue(loaded), fallback = property.GetValue(defaults);
                Assert.That(custom, Is.Not.EqualTo(fallback), $"Custom() must set {property.Name} away from its default");
                Assert.That(property.GetValue(saved), Is.EqualTo(custom), $"{property.Name} is not carried through the window");
            }
        }

        [TestCase("3,4", new[] { 3, 4 })]
        [TestCase(" 2, 3 ,4 ", new[] { 2, 3, 4 })]
        [TestCase("4 3", new[] { 4, 3 })]
        [TestCase("3, 3, 4", new[] { 3, 4 })]
        public static void ChargesAreReadFromAList(string text, int[] expected)
        {
            var shown = Good() with { LibraryPrecursorCharges = text };

            Assert.That(TaskWindowDiaLibrarySearch.Error(shown), Is.Null);
            Assert.That(TaskWindowDiaLibrarySearch.ForSave(shown).LibraryPrecursorCharges, Is.EqualTo(expected));
        }

        [TestCase("fragment", "abc", "fragment tolerance")]
        [TestCase("fragment", "0", "fragment tolerance")]
        [TestCase("ms1", "-5", "MS1 tolerance")]
        [TestCase("seed", "1.5", "seed")]
        [TestCase("charges", "", "charges")]
        [TestCase("charges", "2, x", "charges")]
        [TestCase("charges", "0", "charges")]
        [TestCase("nce", "0", "NCE")]
        [TestCase("nce", "27.5", "NCE")]
        public static void EachFieldThatIsWrongIsNamed(string field, string text, string named)
        {
            var good = Good();
            var shown = field switch
            {
                "fragment" => good with { FragmentTolerancePpm = text },
                "ms1" => good with { Ms1TolerancePpm = text },
                "seed" => good with { ClassifierSeed = text },
                "charges" => good with { LibraryPrecursorCharges = text },
                _ => good with { LibraryNce = text },
            };

            string error = TaskWindowDiaLibrarySearch.Error(shown);

            Assert.That(error, Does.Contain(named));
            Assert.That(error, Does.Contain($"\"{text}\""), "the error quotes what was typed");
        }

        [Test]
        public static void TwoWrongFieldsAreBothNamed()
        {
            var shown = Good() with { FragmentTolerancePpm = "x", LibraryNce = "y" };

            string error = TaskWindowDiaLibrarySearch.Error(shown);

            Assert.That(error, Does.Contain("fragment tolerance").And.Contain("NCE"));
        }

        [Test]
        public static void NumbersAreReadTheSameInEveryCulture()
        {
            var culture = Thread.CurrentThread.CurrentCulture;
            try
            {
                Thread.CurrentThread.CurrentCulture = new CultureInfo("de-DE"); // where "10,5" is ten and a half
                var shown = TaskWindowDiaLibrarySearch.Load(Custom());
                Assert.That(shown.FragmentTolerancePpm, Is.EqualTo("10.5"));
                Assert.That(TaskWindowDiaLibrarySearch.Error(shown), Is.Null);
                Assert.That(TaskWindowDiaLibrarySearch.ForSave(shown).FragmentTolerancePpm, Is.EqualTo(10.5));
            }
            finally
            {
                Thread.CurrentThread.CurrentCulture = culture;
            }
        }
    }
}
