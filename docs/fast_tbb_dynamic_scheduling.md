# Optional TBB dynamic scheduling for Fast

`FastSolverConfig.parallel_chunksize` opts one solve into Numba parallel-loop
chunking. Its default is `None`, so existing callers keep their current Numba
setting. A positive chunk size selects dynamic scheduling only when Numba's
TBB threading layer is active; `workqueue` remains the portable default and
does not provide dynamic scheduling.

Numba requires the TBB package to use this backend. Install it with `pip
install tbb` or `conda install tbb`, then select `NUMBA_THREADING_LAYER=tbb`
before importing Numba or `stiffgwpy_fast`. The setting is per calling thread
and Fast restores the previous value when a solve returns or raises. The
accepted range is `[0, 2**32 - 1]`, matching Numba 0.67's unsigned C setter.
See [Numba scheduling](https://numba.readthedocs.io/en/stable/user/parallel.html)
and [threading-layer setup](https://numba.readthedocs.io/en/stable/user/threading-layer.html).

## Evidence

On the measured Windows 10 / Python 3.11.9 / Numba 0.67.0 host, the exact
candidate patch was applied to a clean detached worktree. A 100-pair run across
five trusted profiles and 21 points in the existing `cr=1` region gave
bitwise-identical result, root, frequency-support, and spectrum hashes. With
two Numba threads, TBB dynamic chunk size 1 had P50 ratios of `0.785–0.897`
(median `0.849`) and P95 ratios of `0.784–0.880` (median `0.838`) versus
workqueue static scheduling. No case exceeded the 5% regression limit.

This certifies only the measured two-thread matrix. Keep the default disabled
and benchmark the exact workload and thread count before opting in. This
scheduling option changes runtime only; it does not expand the independently
validated accuracy domain.

## Reproduction

Set the backend before starting Python, then run the benchmark driver:

```powershell
$env:NUMBA_THREADING_LAYER = 'tbb'
$env:FAST_TBB_BENCH_THREADS = '2'
$env:FAST_TBB_BENCH_PAIRS = '100'
$env:FAST_TBB_BENCH_OUTPUT = 'docs/fast_tbb_dynamic_candidate_clean_commit_100pairs.json'
python scripts/benchmark_fast_tbb_dynamic_chunks_region.py
```

The exact-candidate result is
[`fast_tbb_dynamic_candidate_clean_commit_100pairs.json`](fast_tbb_dynamic_candidate_clean_commit_100pairs.json).
