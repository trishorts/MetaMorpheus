using System;
using System.Linq;
using ILGPU;
using ILGPU.Runtime;
using ILGPU.Runtime.Cuda;
using MzLibUtil;

namespace EngineLayer.ParallelSearch.Scoring
{
    /// <summary>
    /// Hands out an <see cref="ISpectralScorer"/> to each search partition (thread) and owns any
    /// shared, expensive backend state.
    ///
    /// Threading model: the GPU accelerator and the resident experimental spectra are created ONCE
    /// (uploaded to VRAM once) and shared across all partitions — the GPU scorer is a single
    /// thread-safe instance returned to every thread. The CPU provider instead hands each thread its
    /// own lightweight scorer (per-thread match buffer, no shared state).
    /// </summary>
    public interface ISpectralScorerProvider : IDisposable
    {
        /// <summary>Scorer for the calling partition. Do NOT dispose the returned scorer; the provider owns it.</summary>
        ISpectralScorer GetScorer();
        string BackendDescription { get; }
    }

    /// <summary>CPU provider: a fresh (cheap) CpuSpectralScorer per partition.</summary>
    public sealed class CpuScorerProvider : ISpectralScorerProvider
    {
        private readonly SpectralScoringData _data;
        private readonly Tolerance _tolerance;
        public CpuScorerProvider(SpectralScoringData data, Tolerance tolerance) { _data = data; _tolerance = tolerance; }
        public ISpectralScorer GetScorer() => new CpuSpectralScorer(_data, _tolerance);
        public string BackendDescription => "CPU (binary search per fragment)";
        public void Dispose() { }
    }

    /// <summary>
    /// Detects CUDA availability at runtime (direct ILGPU reference on this branch). Cached for the
    /// process lifetime; never throws.
    /// </summary>
    public static class GpuDeviceDetector
    {
        private static readonly object _lock = new();
        private static bool _probed;
        private static bool _available;
        private static string _description = "not probed";

        public static bool IsCudaAvailable { get { Probe(); return _available; } }
        public static string Description { get { Probe(); return _description; } }

        private static void Probe()
        {
            lock (_lock)
            {
                if (_probed) return;
                _probed = true;
                try
                {
                    using var ctx = Context.Create(b => b.Cuda());
                    var dev = ctx.Devices.FirstOrDefault(d => d.AcceleratorType == AcceleratorType.Cuda);
                    _available = dev != null;
                    _description = _available
                        ? $"CUDA device: {dev.Name} ({dev.MemorySize / 1024 / 1024} MB)"
                        : "no CUDA device found";
                }
                catch (Exception ex)
                {
                    _available = false;
                    _description = $"CUDA unavailable: {ex.GetType().Name}: {ex.Message}";
                }
            }
        }
    }

    /// <summary>
    /// Selects the best scorer provider: GPU (CUDA via ILGPU) when available, else CPU. Falls back to
    /// CPU if GPU initialization throws.
    /// </summary>
    public static class SpectralScorerFactory
    {
        public static ISpectralScorerProvider Create(SpectralScoringData data, Tolerance productTolerance, bool preferCpu = false)
        {
            if (!preferCpu && GpuDeviceDetector.IsCudaAvailable)
            {
                try { return new GpuScorerProvider(data, productTolerance); }
                catch (Exception ex)
                {
                    Console.Error.WriteLine($"[ParallelSearch][GPU] init failed, falling back to CPU: {ex.Message}");
                }
            }
            return new CpuScorerProvider(data, productTolerance);
        }

        public static string DescribeBackend(bool preferCpu = false)
        {
            if (preferCpu) return "CPU (forced)";
            return GpuDeviceDetector.IsCudaAvailable
                ? $"GPU ({GpuDeviceDetector.Description})"
                : $"CPU ({GpuDeviceDetector.Description})";
        }
    }
}
