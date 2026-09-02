using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    // SCRATCH -- diagnostic only, delete after use. Re-measures acquisition structure for PXD053713.
    // M9 established HCD-parent/electron-child pairing for PXD017646; these are a different lab and
    // instrument method, so it must be measured again rather than inherited. Getting it wrong is
    // silent: the c/z ions localization needs are never matched and it reads as hard data.
    [TestFixture]
    public class ScratchScanStructure
    {
        private const string OutPath =
            "C:/Users/trish/AppData/Local/Temp/claude/E--CodeReview-phred/71179d18-0c21-49c0-918d-710f0a897173/scratchpad/pxd053713_scans.txt";
        private static void Log(string s) => File.AppendAllText(OutPath, s + Environment.NewLine);

        [Test]
        public static void DumpPxd053713()
        {
            var dir = "F:/phred-data/PXD053713/raw";
            foreach (var path in Directory.GetFiles(dir, "*.raw").OrderBy(p => p).Take(3))
            {
                var scans = MsDataFileReader.GetDataFile(path).LoadAllStaticData().GetAllScansList();
                var byNum = scans.ToDictionary(s => s.OneBasedScanNumber);
                Log("=== " + Path.GetFileName(path) + "   total " + scans.Count);

                var buckets = new Dictionary<string, int>();
                foreach (var s in scans.Where(s => s.MsnOrder == 2))
                {
                    string d = s.DissociationType?.ToString() ?? "null";
                    string parent = "missing";
                    if (s.OneBasedPrecursorScanNumber.HasValue &&
                        byNum.TryGetValue(s.OneBasedPrecursorScanNumber.Value, out var p))
                    {
                        parent = "MS" + p.MsnOrder +
                                 (p.MsnOrder == 2 ? "/" + (p.DissociationType?.ToString() ?? "null") : "");
                    }
                    string k = d + "  ->  precursor is " + parent;
                    buckets[k] = buckets.TryGetValue(k, out var c) ? c + 1 : 1;
                }
                Log("   MS1 " + scans.Count(s => s.MsnOrder == 1));
                foreach (var kv in buckets.OrderByDescending(k => k.Value))
                    Log(string.Format("   {0,-44} n={1}", kv.Key, kv.Value));
            }
        }
    }
}
