using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using NUnit.Framework;
using PredictionClients.LocalModels;
using PredictionClients.MixedModels;

namespace Test
{
    /// <summary>
    /// Proves MetaMorpheus can obtain predicted fragmentation spectra "at will" from the mzLib
    /// PredictionClients package the same way it would from a Koina-supported model: pass in
    /// peptides, get back combined spectra carrying both primary (b/y, from Koina/Prosit) and
    /// internal-cleavage ions (from the local ONNX model shipped in the package).
    ///
    /// These are tagged [Explicit] because they reach the live Koina HTTP endpoint and load the
    /// local ONNX model, so they are not part of the normal CI run. Run them on demand to validate
    /// the end-to-end path inside MetaMorpheus, including that the ONNX model and onnxruntime
    /// deployed correctly via the mzLib NuGet package.
    /// </summary>
    [TestFixture]
    public static class PredictedInternalIonSpectraTests
    {
        [Test]
        public static void OnnxModel_DeployedNextToPredictionClients_Resolves()
        {
            // The internal model ships in the mzLib package as a contentFile copied to
            // <output>\LocalModels\. If this fails, the package didn't deploy the model.
            var path = InternalFragmentIntensityModel.DefaultOnnxModelPath;
            Assert.That(System.IO.File.Exists(path), Is.True,
                $"Internal-fragment ONNX model was not found at the resolved path: {path}. " +
                "Check that the mzLib package shipped it via contentFiles/copyToOutput.");
        }

        [Test, Explicit("Requires live Koina endpoint and local ONNX model")]
        public static async Task PredictPeptides_ReturnsCombinedPrimaryAndInternalSpectra()
        {
            var peptides = new List<string> { "PEPTIDEK", "ELVISLIVESK" };

            var (spectra, warning) = await CombinedLibraryModel.PredictPrimaryAndInternalSpectraAsync(peptides);

            if (warning != null)
                Assert.That(warning.Message, Does.Not.Contain("failed"),
                    $"A prediction component failed unexpectedly: {warning.Message}");

            Assert.That(spectra.Count, Is.EqualTo(2));
            foreach (var spectrum in spectra)
            {
                Assert.That(spectrum.ChargeState, Is.EqualTo(2),
                    "Default precursor charge should be 2 when none is supplied");

                bool hasPrimary = spectrum.MatchedFragmentIons.Any(f => !f.IsInternalFragment);
                bool hasInternal = spectrum.MatchedFragmentIons.Any(f => f.IsInternalFragment);

                Assert.That(hasPrimary, Is.True,
                    $"{spectrum.Name}: expected primary (b/y) ions from the Koina model");
                Assert.That(hasInternal, Is.True,
                    $"{spectrum.Name}: expected internal-cleavage ions from the local ONNX model");
            }
        }
    }
}
