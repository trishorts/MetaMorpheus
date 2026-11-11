using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;

namespace EngineLayer
{
    /// <summary>
    /// What MetaMorpheus records beside an experimental design it wrote from an SDRF (QuantProject QP-S41, ruled
    /// 2026-11-07): the design file and the SDRF, each by its sha256, and how the condition was built, including any
    /// grouping the user stated (ruled 2026-11-11). It is how a later run tells a design MetaMorpheus wrote from one
    /// the user wrote or edited: while the design's sha256 is the one recorded here, the design is MetaMorpheus's, and
    /// a run writes it again for its own spectra files; once edited, it is the user's and is never touched.
    ///
    /// A plain text file of tab-separated keys and values beside the design, so a person reading the folder sees why
    /// the design is there.
    /// </summary>
    public sealed class SdrfDesignMarker
    {
        /// <summary>The marker's file name, beside the design.</summary>
        public const string FileName = "DesignFromSdrf.txt";

        private const string Heading = "# MetaMorpheus wrote the experimental design beside this file from an SDRF. While the design is as written, " +
                                       "MetaMorpheus writes it again for each run; once edited, the design is yours and is never touched.";

        public SdrfDesignMarker(string designFileName, string designSha256, string sdrfPath, string sdrfSha256, string condition)
        {
            DesignFileName = designFileName;
            DesignSha256 = designSha256;
            SdrfPath = sdrfPath;
            SdrfSha256 = sdrfSha256;
            Condition = condition;
        }

        /// <summary>The design's file name: ExperimentalDesign.tsv or TmtDesign.txt.</summary>
        public string DesignFileName { get; }

        /// <summary>The design's sha256 as MetaMorpheus wrote it, lowercase hex.</summary>
        public string DesignSha256 { get; }

        /// <summary>The SDRF the design was made from, as a full path.</summary>
        public string SdrfPath { get; }

        /// <summary>The SDRF's sha256 when the design was made, lowercase hex.</summary>
        public string SdrfSha256 { get; }

        /// <summary>How the condition was built, in words, with any grouping the user stated.</summary>
        public string Condition { get; }

        /// <summary>The marker beside a design in <paramref name="designDirectory"/>.</summary>
        public static string PathIn(string designDirectory) => Path.Combine(designDirectory, FileName);

        /// <summary>The line results.txt carries for a design written from an SDRF.</summary>
        public string Describe() =>
            $"Experimental design ({DesignFileName}) written by MetaMorpheus from the SDRF {SdrfPath} (sha256 {SdrfSha256}); condition: {Condition}.";

        /// <summary>
        /// True when the design this marker names is in <paramref name="designDirectory"/> and is byte for byte what
        /// MetaMorpheus wrote. False when it is missing, unreadable or edited.
        /// </summary>
        public bool DescribesUneditedDesignIn(string designDirectory)
        {
            string design = Path.Combine(designDirectory, DesignFileName);
            try
            {
                return File.Exists(design) && string.Equals(Sha256Of(design), DesignSha256, StringComparison.OrdinalIgnoreCase);
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                return false;
            }
        }

        /// <summary>Writes the marker into <paramref name="designDirectory"/>, replacing any there.</summary>
        public void WriteTo(string designDirectory)
        {
            var lines = new List<string>
            {
                Heading,
                "design\t" + DesignFileName,
                "design sha256\t" + DesignSha256,
                "sdrf\t" + SdrfPath,
                "sdrf sha256\t" + SdrfSha256,
                "condition\t" + Condition,
            };
            File.WriteAllText(PathIn(designDirectory), string.Join("\n", lines) + "\n", new UTF8Encoding(false));
        }

        /// <summary>
        /// The marker in <paramref name="designDirectory"/>, or null when there is none or it cannot be read as one: a
        /// marker that is missing a value is no marker, so the design beside it is treated as the user's.
        /// </summary>
        public static SdrfDesignMarker ReadFrom(string designDirectory)
        {
            string path = PathIn(designDirectory);
            Dictionary<string, string> values;
            try
            {
                if (!File.Exists(path))
                    return null;
                values = File.ReadAllLines(path)
                    .Where(l => !l.StartsWith("#") && l.Contains('\t'))
                    .Select(l => l.Split('\t', 2))
                    .GroupBy(kv => kv[0], StringComparer.Ordinal)
                    .ToDictionary(g => g.Key, g => g.First()[1], StringComparer.Ordinal);
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                return null;
            }

            string Value(string key) => values.TryGetValue(key, out var v) && !string.IsNullOrWhiteSpace(v) ? v : null;
            string design = Value("design");
            string designSha = Value("design sha256");
            string sdrf = Value("sdrf");
            string sdrfSha = Value("sdrf sha256");
            string condition = Value("condition");
            if (design == null || designSha == null || sdrf == null || sdrfSha == null || condition == null)
                return null;
            return new SdrfDesignMarker(design, designSha, sdrf, sdrfSha, condition);
        }

        /// <summary>A file's sha256, lowercase hex, read from the bytes on disk.</summary>
        public static string Sha256Of(string path)
        {
            using var stream = File.OpenRead(path);
            return Convert.ToHexString(SHA256.HashData(stream)).ToLowerInvariant();
        }
    }
}
