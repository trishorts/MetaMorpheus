using System;
using Chemistry;
using ILGPU;
using ILGPU.Runtime;
using ILGPU.Runtime.Cuda;
using MzLibUtil;

namespace EngineLayer.ParallelSearch.Scoring
{
    /// <summary>
    /// GPU (CUDA via ILGPU, FP64) implementation of <see cref="ISpectralScorer"/>.
    ///
    /// The experimental spectra (<see cref="SpectralScoringData"/>) are uploaded to VRAM ONCE in the
    /// constructor and stay resident for the scorer's lifetime — reused across every transient
    /// database. Each <see cref="ScoreBatch"/> uploads only the per-batch theoretical fragments and
    /// work items, runs a block-per-work-item kernel that caches the candidate scan's experimental
    /// fragment slice in shared memory (so the peptide's fragment binary-searches hit shared memory
    /// instead of scattering across global memory), and returns the matched experimental index per
    /// (work item, fragment). The host then rebuilds the exact same MatchedFragmentIon inputs as the
    /// CPU path from the same arrays, so results are BYTE-IDENTICAL.
    ///
    /// Only PPM product tolerance is supported (the kernel replicates mzLib PpmTolerance.Within
    /// exactly); other tolerance types fall back to CPU via the provider.
    ///
    /// A single instance is shared across all search threads; <see cref="ScoreBatch"/> is serialized
    /// with a lock (one GPU device).
    /// </summary>
    public sealed class GpuSpectralScorer : ISpectralScorer
    {
        private const int GroupSize = 64;
        private const int MaxSlice = 512;   // experimental fragments per scan cached in shared memory

        private readonly object _gpuLock = new();
        private readonly SpectralScoringData _data;
        private readonly double _ppm;

        private readonly Context _context;
        private readonly Accelerator _accelerator;

        // resident experimental spectra (uploaded once)
        private readonly MemoryBuffer1D<double, Stride1D.Dense> _dExpMass;
        private readonly MemoryBuffer1D<int, Stride1D.Dense> _dExpCharge;
        private readonly MemoryBuffer1D<int, Stride1D.Dense> _dScanFragOff;
        private readonly MemoryBuffer1D<int, Stride1D.Dense> _dScanPrecCharge;

        // per-batch buffers (grown as needed)
        private MemoryBuffer1D<int, Stride1D.Dense> _dWorkScan;
        private MemoryBuffer1D<int, Stride1D.Dense> _dWorkSlot;
        private MemoryBuffer1D<int, Stride1D.Dense> _dPepFragOff;
        private MemoryBuffer1D<double, Stride1D.Dense> _dFragMass;
        private MemoryBuffer1D<int, Stride1D.Dense> _dOrder;
        private MemoryBuffer1D<int, Stride1D.Dense> _dOut;

        private int[] _order = new int[1024];
        private FragmentMatch[] _matchBuf = new FragmentMatch[128];

        private readonly Action<KernelConfig,
            ArrayView<int>, ArrayView<int>, ArrayView<int>, ArrayView<int>, ArrayView<double>,
            ArrayView<int>, ArrayView<double>, ArrayView<int>, ArrayView<int>,
            double, int, ArrayView<int>> _kernel;

        public string BackendDescription { get; }

        public GpuSpectralScorer(SpectralScoringData data, double ppmTolerance)
        {
            _data = data ?? throw new ArgumentNullException(nameof(data));
            _ppm = ppmTolerance;

            _context = Context.Create(b => b.Cuda());
            var device = System.Linq.Enumerable.First(_context.Devices); // CUDA (detector already confirmed)
            _accelerator = device.CreateAccelerator(_context);
            BackendDescription = $"GPU CUDA (ILGPU FP64): {_accelerator.Name}";

            _dExpMass = _accelerator.Allocate1D(data.FragmentMonoMasses);
            _dExpCharge = _accelerator.Allocate1D(data.FragmentCharges);
            _dScanFragOff = _accelerator.Allocate1D(data.ScanFragmentOffsets);
            _dScanPrecCharge = _accelerator.Allocate1D(data.ScanPrecursorCharges);

            _kernel = _accelerator.LoadStreamKernel<
                ArrayView<int>, ArrayView<int>, ArrayView<int>, ArrayView<int>, ArrayView<double>,
                ArrayView<int>, ArrayView<double>, ArrayView<int>, ArrayView<int>,
                double, int, ArrayView<int>>(MatchKernel);
        }

        public void ScoreBatch(ScoringBatch batch, IScoringSink sink)
        {
            if (batch == null) throw new ArgumentNullException(nameof(batch));
            if (sink == null) throw new ArgumentNullException(nameof(sink));
            int count = batch.WorkItemCount;
            if (count == 0) return;

            lock (_gpuLock)
            {
                // scan-sorted processing order (locality: neighbor blocks -> neighbor scans)
                if (_order.Length < count) _order = new int[Math.Max(count, _order.Length * 2)];
                for (int i = 0; i < count; i++) _order[i] = i;
                Array.Sort(_order, 0, count, new ScanComparer(batch.WorkScanIndex));

                // output stride = max theoretical fragments per peptide in this batch
                int stride = 1;
                for (int s = 0; s < batch.PeptideCount; s++)
                {
                    int len = batch.PeptideFragmentOffsets[s + 1] - batch.PeptideFragmentOffsets[s];
                    if (len > stride) stride = len;
                }

                // upload (resident spectra already on device); whole host arrays — kernel reads only
                // the valid ranges. Device buffers are >= host array length.
                _dWorkSlot = Ensure(_dWorkSlot, batch.WorkPeptideSlot.Length); _dWorkSlot.CopyFromCPU(batch.WorkPeptideSlot);
                _dWorkScan = Ensure(_dWorkScan, batch.WorkScanIndex.Length); _dWorkScan.CopyFromCPU(batch.WorkScanIndex);
                _dPepFragOff = Ensure(_dPepFragOff, batch.PeptideFragmentOffsets.Length); _dPepFragOff.CopyFromCPU(batch.PeptideFragmentOffsets);
                _dFragMass = EnsureD(_dFragMass, batch.FragmentNeutralMasses.Length); _dFragMass.CopyFromCPU(batch.FragmentNeutralMasses);
                _dOrder = Ensure(_dOrder, _order.Length); _dOrder.CopyFromCPU(_order);

                long outLen = (long)count * stride;
                _dOut = Ensure(_dOut, outLen);

                _kernel(new KernelConfig(count, GroupSize),
                    _dOrder.View, _dWorkSlot.View, _dWorkScan.View, _dPepFragOff.View, _dFragMass.View,
                    _dScanFragOff.View, _dExpMass.View, _dExpCharge.View, _dScanPrecCharge.View,
                    _ppm, stride, _dOut.View);
                _accelerator.Synchronize();

                int[] outArr = _dOut.GetAsArray1D(); // length = buffer length (>= outLen)

                // drain in ORIGINAL work-item order (matches CPU scorer -> PSM update order -> byte-identical)
                for (int w = 0; w < count; w++)
                {
                    int slot = batch.WorkPeptideSlot[w];
                    int fragCount = batch.PeptideFragmentOffsets[slot + 1] - batch.PeptideFragmentOffsets[slot];
                    long baseOut = (long)w * stride;

                    int matchCount = 0;
                    for (int local = 0; local < fragCount; local++)
                    {
                        int gidx = outArr[baseOut + local];
                        if (gidx < 0) continue;
                        if (_matchBuf.Length <= matchCount) Array.Resize(ref _matchBuf, _matchBuf.Length * 2);
                        int ch = _data.FragmentCharges[gidx];
                        _matchBuf[matchCount++] = new FragmentMatch(
                            local, _data.FragmentMonoMasses[gidx].ToMz(ch), _data.FragmentIntensities[gidx], ch);
                    }
                    if (matchCount > 0)
                        sink.AcceptWorkItem(w, _matchBuf, matchCount);
                }
            }
        }

        // Block per work item: cache the candidate scan's experimental fragment slice in shared
        // memory, then each of the peptide's fragments binary-searches the cached slice. Output is the
        // matched global experimental index per (work item, local fragment), or -1.
        private static void MatchKernel(
            ArrayView<int> order,
            ArrayView<int> workSlot, ArrayView<int> workScan,
            ArrayView<int> pepFragOff, ArrayView<double> fragMass,
            ArrayView<int> scanFragOff, ArrayView<double> expMass, ArrayView<int> expCharge,
            ArrayView<int> scanPrecCharge,
            double ppm, int stride,
            ArrayView<int> outPerLocalIdx)
        {
            int w = order[Grid.IdxX];
            int scan = workScan[w];
            int sLo = scanFragOff[scan];
            int slen = scanFragOff[scan + 1] - sLo;
            int precCharge = scanPrecCharge[scan];
            int precAbs = precCharge < 0 ? -precCharge : precCharge;

            var slice = SharedMemory.Allocate<double>(MaxSlice);
            bool useShared = slen <= MaxSlice;
            if (useShared)
                for (int i = Group.IdxX; i < slen; i += Group.DimX)
                    slice[i] = expMass[sLo + i];
            Group.Barrier();

            int slot = workSlot[w];
            int fStart = pepFragOff[slot];
            int fcount = pepFragOff[slot + 1] - fStart;

            for (int local = Group.IdxX; local < fcount; local += Group.DimX)
            {
                long outPos = (long)w * stride + local;
                outPerLocalIdx[outPos] = -1;

                double m = fragMass[fStart + local];
                if (m != m || slen == 0) continue;

                int L = 0, R = slen;
                while (L < R)
                {
                    int mid = L + ((R - L) >> 1);
                    double v = useShared ? slice[mid] : expMass[sLo + mid];
                    if (v < m) L = mid + 1; else R = mid;
                }
                int ip = L, j;
                if (ip == 0) j = 0;
                else if (ip == slen) j = slen - 1;
                else
                {
                    double vL = useShared ? slice[ip - 1] : expMass[sLo + ip - 1];
                    double vR = useShared ? slice[ip] : expMass[sLo + ip];
                    j = ((vR - m) < (m - vL)) ? ip : ip - 1;
                }

                double em = useShared ? slice[j] : expMass[sLo + j];
                int gidx = sLo + j;
                double diff = (em - m) / m * 1e6; if (diff < 0) diff = -diff;
                int ec = expCharge[gidx]; int ecAbs = ec < 0 ? -ec : ec;
                if (diff <= ppm && ecAbs <= precAbs)
                    outPerLocalIdx[outPos] = gidx;
            }
        }

        private MemoryBuffer1D<int, Stride1D.Dense> Ensure(MemoryBuffer1D<int, Stride1D.Dense> buf, long len)
        {
            if (buf != null && buf.Length >= len) return buf;
            buf?.Dispose();
            return _accelerator.Allocate1D<int>(Math.Max(len, 1024));
        }

        private MemoryBuffer1D<double, Stride1D.Dense> EnsureD(MemoryBuffer1D<double, Stride1D.Dense> buf, long len)
        {
            if (buf != null && buf.Length >= len) return buf;
            buf?.Dispose();
            return _accelerator.Allocate1D<double>(Math.Max(len, 1024));
        }

        private sealed class ScanComparer : System.Collections.Generic.IComparer<int>
        {
            private readonly int[] _scan;
            public ScanComparer(int[] scan) { _scan = scan; }
            public int Compare(int a, int b) => _scan[a].CompareTo(_scan[b]);
        }

        public void Dispose()
        {
            _dExpMass?.Dispose(); _dExpCharge?.Dispose(); _dScanFragOff?.Dispose(); _dScanPrecCharge?.Dispose();
            _dWorkScan?.Dispose(); _dWorkSlot?.Dispose(); _dPepFragOff?.Dispose(); _dFragMass?.Dispose();
            _dOrder?.Dispose(); _dOut?.Dispose();
            _accelerator?.Dispose();
            _context?.Dispose();
        }
    }

    /// <summary>Owns the shared GPU scorer (one accelerator, resident spectra). PPM tolerance only.</summary>
    public sealed class GpuScorerProvider : ISpectralScorerProvider
    {
        private readonly GpuSpectralScorer _scorer;
        public GpuScorerProvider(SpectralScoringData data, Tolerance productTolerance)
        {
            if (productTolerance is not PpmTolerance ppm)
                throw new NotSupportedException("GPU scorer supports PPM product tolerance only.");
            _scorer = new GpuSpectralScorer(data, ppm.Value);
        }
        public ISpectralScorer GetScorer() => _scorer;
        public string BackendDescription => _scorer.BackendDescription;
        public void Dispose() => _scorer.Dispose();
    }
}
