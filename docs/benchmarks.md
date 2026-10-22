# Benchmarks

Status: current summary of the audited scope; results remain commit- and protocol-bound

Date: 2026-10-22

Code version: the corrected formal fast profile below used
`6707572e73d576a0ccf294172579cedc26e5e394`; the earlier stage profile used
`cdbdb387479e11cd87a38d76a0e58348ca13e70f`. The historical Round 68 run used
`c7d4766f15c1cda4693d306597723de50abd4554`. The PTA cache comparison used
source commit `c185df58938a05c208410a5cd47faafbc5d7fa2c`. Each run artifact
records its own commit.

## 2026-10-22 corrected formal fast profile

Review of the profile driver found that it called the low-level solver without
passing the selected preset's outer tolerance. The low-level API defaulted to
`1e-4`, while the formal `fast` preset specifies `1e-6`. The driver now uses
the selected preset tolerance unless `--tol` is explicitly supplied, and each
record includes the tolerance, frequency grid, quadrature, Numba thread layer,
and resource telemetry.

This rerun uses the formal `fast` tolerance `1e-6`, PCHIP, the `goal` grid,
exact kink split, `h=0.005`, `col_step=8`, `z_tail=5`, and `phase_max=0.25`.
Outer full-solve reuse is enabled and the parallel chunksize is `0` (static
scheduling). Each case ran in a fresh process for 26 calls, with the first
process/JIT call reported separately and the remaining 25 used for warm
statistics. The raw JSON stores explicit warm-only median and p95 fields as
well as the cold first-call duration. Resources were fixed to 2 Numba
`workqueue` threads, CPU affinity `[0, 1]`, one BLAS/OMP/NumExpr thread, and
one process. The JSON records the complete solver configuration and these
environment values.

| Case | First process call | Warm median | Warm p95 | Tensor median | Frequencies | Propagation calls |
|---|---:|---:|---:|---:|---:|---:|
| default | 292.5 ms | 5.987 ms | 6.882 ms | 2.905 ms | 76 | 1 |
| high-T | 301.2 ms | 10.049 ms | 10.848 ms | 5.397 ms | 76 | 2 |
| stiff | 306.2 ms | 9.179 ms | 9.819 ms | 5.653 ms | 76 | 2 |
| high-kappa | 295.8 ms | 9.458 ms | 10.112 ms | 5.172 ms | 77 | 2 |

All 26 records per case have identical hashes for frequencies, spectrum,
`DN_gw`, `g2`, and `w2`. This checks repeatability for each measured case; it
does not add Oracle accuracy coverage or establish a cross-machine speed ratio.
Raw records: [`default`](profile_fast_breakdown_6707572_20261022_default.json),
[`high-T`](profile_fast_breakdown_6707572_20261022_highT.json),
[`stiff`](profile_fast_breakdown_6707572_20261022_stiff.json), and
[`high-kappa`](profile_fast_breakdown_6707572_20261022_high_kappa.json).

Reproduce a case with the fixed-resource wrapper; repeat for each case name:

```powershell
python scripts/run_fixed_profile.py --threads 2 --affinity-count 2 --reps 26 --case default --kink-split --json docs/profile_fast_breakdown_6707572_20261022_default.json
```

## Earlier 2026-10-21 stage profile (historical low-level tolerance)

The 2026-10-21 run used PCHIP, the goal frequency grid, exact kink split, and
`phase_max=0.25`, with 2 Numba `workqueue` threads and CPU affinity `[0, 1]`.
However, the driver called the low-level solver without an explicit `tol`, so
it used `1e-4` rather than the formal preset's `1e-6`. Keep these values as a
dated lower-tolerance stage profile; do not read them as formal `fast` preset
runtime or compare them directly with the corrected profile above. Each case
ran in a fresh process for 26 calls; the first process call is reported
separately and the warm statistics use the remaining 25.

| Case | First process call | Warm median | Warm p95 | Tensor median | Background median | Frequencies | Propagation calls |
|---|---:|---:|---:|---:|---:|---:|---:|
| default | 6,370.5 ms | 6.790 ms | 8.002 ms | 3.456 ms | 0.605 ms | 76 | 1 |
| high-T | 486.2 ms | 10.284 ms | 12.654 ms | 5.976 ms | 0.725 ms | 76 | 2 |
| stiff | 316.1 ms | 10.095 ms | 11.734 ms | 6.434 ms | 0.655 ms | 76 | 2 |
| high-kappa | 338.9 ms | 10.084 ms | 11.221 ms | 5.861 ms | 0.729 ms | 77 | 2 |

The first process call includes process startup and any JIT compilation needed
for that run; it is not included in the warm statistics. All 26 records in
each case have identical digests for frequencies, spectrum, `DN_gw`, `g2`, and
`w2`. These are local single-machine measurements, with no matched LSODA
comparison, and do not establish a cross-machine speed ratio.

Raw per-case records: [`default`](profile_fast_breakdown_cdbdb38_20261021_default.json),
[`high-T`](profile_fast_breakdown_cdbdb38_20261021_highT.json),
[`stiff`](profile_fast_breakdown_cdbdb38_20261021_stiff.json), and
[`high-kappa`](profile_fast_breakdown_cdbdb38_20261021_high_kappa.json).

## Tensor sample-index optimization A/B (2026-10-21)

The kernel now tracks the next coarse output column and advances that index
only when the integration reaches it. This removes a modulo operation from
each inner integration step; transfer arithmetic and output-slot selection are
unchanged. The baseline kernel was loaded from commit
`68b17acf8e049bf17ed100422da0f2f3553f15e9`; the optimized implementation is
commit `ae460e9a1aec46a902033ad287d5f424413c7972`.

The full-path A/B alternated baseline and candidate 80 times in each of the
six fast benchmark regimes, using two Numba `workqueue` threads, static
chunks, and CPU affinity `[0, 1]`. All five output-array digests matched for
every run, and convergence/failure states matched. Paired percentage changes
are candidate relative to baseline; confidence intervals use 20,000 bootstrap
resamples of paired round indices.

| Case | Total paired median change | Bootstrap 95% interval |
|---|---:|---:|
| default | −2.52% | −4.79% to +1.73% |
| high-T | −5.48% | −7.98% to −3.66% |
| stiff | −6.97% | −8.65% to −4.78% |
| high-kappa | −5.14% | −7.61% to −3.23% |
| low-T | −7.89% | −9.10% to −5.04% |
| low-r | −6.63% | −11.07% to −1.75% |
| Round-wise geometric mean | −5.44% | −7.18% to −3.55% |

The tensor stage's round-wise geometric-mean paired median improved by 9.29%
(95% interval 7.51% to 10.30% faster). Additional five-repeat full-path
comparisons at `col_step=1` and `col_step=4` also matched all output digests;
those short compatibility checks are not speed claims. The six-regime A/B
adds no oracle-accuracy points or broad parameter-space certification. Raw
samples, hashes, and the benchmark script are in
[`benchmark_kernel_assembly_stride_20261021.json`](benchmark_kernel_assembly_stride_20261021.json)
and [`benchmark_kernel_assembly_stride.py`](../scripts/benchmark_kernel_assembly_stride.py).

## Thread-safe backend concurrent-call A/B (2026-10-21)

The wrapper serializes solver calls when Numba selects `workqueue`, which is
not thread-safe for concurrent parallel calls. With an initialized TBB/OpenMP
layer, explicit immutable solver configs can use independent thread masks and
run concurrently. The benchmark compared two full-path calls (`default` and
`highT`) in a two-thread Python pool, each solver using two Numba threads and
four-CPU affinity `[0, 1, 2, 3]`. It alternated the previous serialized-lock
policy against the backend-aware policy for 80 paired rounds on TBB.

| Batch policy | Median wall time for both calls |
|---|---:|
| Serialized | 13.178 ms |
| Concurrent TBB | 11.618 ms |

The paired median wall-time change was −14.93% (bootstrap 95% interval
−18.17% to −9.79%). All five output arrays matched bit-for-bit for both cases
in every round, and both calls converged without failure. This is a local
two-call throughput comparison, not a general MCMC speed ratio or a single-call
latency result. It adds no oracle-accuracy coverage. Raw samples and source
provenance are in
[`benchmark_backend_concurrency_20261021.json`](benchmark_backend_concurrency_20261021.json);
the reproduction script is
[`benchmark_backend_concurrency.py`](../scripts/benchmark_backend_concurrency.py).

## Reference-spectrum quadrature breakpoints (2026-10-22)

The interpolation-error estimate integrates a full-grid PCHIP and a
half-density PCHIP. Previously, both adaptive quadratures split at every
full-grid knot. The coarse spline only has every-other knot, so the intervening
fine-grid nodes created unnecessary intervals. For even input lengths,
`x[::2]` also omitted the final endpoint and extrapolated the coarse spline over
the last interval. The candidate passes each interpolant's own knots to
`quad` and includes the final input point when needed.

The benchmark used 40 alternating pairs per case: three deterministic
synthetic spectra and a 76-frequency spectrum generated by the formal-fast
default case. It times `integrate_spectrum` alone; the fast-spectrum case is
input data, not an end-to-end solver measurement.

| Input | QUADPACK function calls, previous → candidate | Candidate/previous median | Paired bootstrap 95% interval |
|---|---:|---:|---:|
| Power law, 401 points | 16,800 → 12,600 | 0.7555 | 0.7507–0.7578 |
| Peaked spectrum, 247 points | 10,332 → 7,749 | 0.7513 | 0.7476–0.7549 |
| Irregular grid, 240 points | 10,038 → 7,539 | 0.7618 | 0.7492–0.7825 |
| Formal-fast default spectrum, 76 points | 3,150 → 2,373 | 0.7657 | 0.7504–0.7708 |

The function-call count fell about 25% in every case. The integral and main
quadrature-error outputs matched exactly. On odd-length inputs, the
interpolation-error estimate changed by at most `2.62e-13` absolute. On the
240-point even-length case, including the omitted final endpoint changed that
estimate by `6.83e-7` absolute; this corrects the previous extrapolated coarse
integral. These measurements apply to this helper on the recorded Windows/SciPy
environment and do not establish a full reference or fast solver speedup. Raw
samples, input digests, and environment details are in
[`benchmark_reference_quad_breakpoints_20261022.json`](benchmark_reference_quad_breakpoints_20261022.json);
the reproduction script is
[`benchmark_reference_quad_breakpoints.py`](../scripts/benchmark_reference_quad_breakpoints.py).

## Historical Round 68 full-path profile

The only user-facing fast profile uses `frequency_quadrature=pchip` with
`h=0.005`, `col_step=8`, `z_tail=5`, `phase_max=0.25`, exact kink split, and
the goal frequency grid. The 2026-09-23 Round 68 run at solver commit
`c7d4766f15c1cda4693d306597723de50abd4554` remains a matched historical
baseline.

Protocol: 25 repeats, Numba=2/workqueue, `FAST_THREADS=2`, one BLAS thread,
one worker/process, Windows CPU affinity `[0, 1]`, Python 3.11.9, NumPy 2.4.4,
SciPy 1.17.1, Numba 0.67.0. Warm runtime excludes the cold/JIT call.

| Case | First sample | Warm median | Warm p95 | Frequencies | Status |
|---|---:|---:|---:|---:|---|
| default | 545.2 ms | 7.189 ms | 7.916 ms | 76 | converged |
| low-T | 8.3 ms | 7.116 ms | 7.782 ms | 76 | converged |
| high-T | 10.5 ms | 10.475 ms | 10.890 ms | 76 | converged |
| stiff | 11.4 ms | 10.260 ms | 10.936 ms | 76 | converged |
| low-r | 6.8 ms | 5.923 ms | 6.313 ms | 76 | converged |
| high-kappa | 9.7 ms | 9.929 ms | 11.062 ms | 77 | converged |

Only the default row is the process-cold first call and includes JIT startup;
the other first-sample values were measured after that process initialized
the kernels. All six rows report no failure reason. The raw run is
[`benchmark_head_matrix_round68_20260923.json`](benchmark_head_matrix_round68_20260923.json);
stage profiles are `profile_fast_breakdown_round68_20260923_*.json`. This
profile has no matched-resource LSODA measurement, so this page does not claim
a current LSODA speedup ratio. The stable `<4 ms/point` target remains unmet.

## Later PTA likelihood measurement

Commit `0a03660702493caaf8f0bd4034187324e3617ce0` moved fixed IPTA KDE and
NANOGrav spline preparation into each likelihood instance's `initialize()`
hook. The benchmark script loads the pre-cache classes from commit `1f648a5`,
initializes both implementations, then alternates 25 direct `log_likelihood`
calls per implementation after one warm-up call. It records these medians:

| Adapter | Before median | Cached median | Cached initialization | Serialized cache per instance |
|---|---:|---:|---:|---:|
| IPTA | 34.033 ms | 27.568 ms | 66.933 ms | 4,321,682 bytes |
| NANOGrav | 4.574 ms | 0.203 ms | 7.198 ms | 5,601,185 bytes |

All 25 timed outputs matched exactly between implementations. A separate
10-case nuisance/theory sweep per adapter also matched exactly, with maximum
absolute difference zero. The run records its environment and results in
[`benchmark_pta_likelihood_cache_20261021.json`](benchmark_pta_likelihood_cache_20261021.json);
the JSON retains all timed samples and outputs so the medians and exactness
checks can be recomputed. To reproduce it from the repository root, run
`python scripts/benchmark_pta_likelihood_cache.py`. By default, the script
writes a new timestamped JSON file; `--output` selects a specific path.
These local adapter timings exclude initialization from each evaluation and
do not measure end-to-end Cobaya MCMC runtime or change the scope of the Round
68 solver profile.

## Scoped accuracy audit (dated 2026-09-12)

This four-regime PCHIP/oracle baseline is summarized in
[`fast_v02_audit_report.md`](fast_v02_audit_report.md) and uses its recorded
2026-09-12 evidence rather than a rerun at the Round 68 or later commits. The
newer two-edge `DN_re` results are summarized in
[`accuracy.md`](accuracy.md) and the
[`Oracle C assessment`](oracle_c_wkb_assessment.md):

| Measure | Observed | Scope / decision |
|---|---:|---|
| Oracle C WKB `DN_gw` relative error, four named points | median `1.1e-5`, max `1.66e-4` | Passes the `<2e-4` gate on these four points |
| Independent reference, same grid, `z_tail=8`, six named points | median `2.93e-4`, max `2.96e-4` | Above the gate; tail-convention limited in this comparison |
| Default point vs independent reference at `z_tail=10` | `4.96e-6` | Passes the `<2e-4` gate for this point |
| Nested native-frequency refinement, six named points | relative changes `7.4e-12`–`2.8e-9` | Native-grid refinement converged in this test |
| Stability screen | 0 numerical failures; 3 explicit physical guards / 24 points | No silent fallback in the audited sample |
| Formal-fast full parameter-space accuracy | — | `NOT VERIFIED` |

The scoped accuracy comparison supports PCHIP for the tested cases. The
independent reference at `z_tail=8` has a frozen-handoff tail systematic of the
same order as its difference from fast; the deeper `z_tail=10` default-point
comparison and Oracle C anchor agree closely. See the audit report for the
full attribution and remaining uncertainty.

## PCHIP and Simpson comparison

The paired 50-repeat comparison in
[`fast_quadrature_default_switch_assessment.md`](fast_quadrature_default_switch_assessment.md)
reports PCHIP / Simpson warm-median ratios of `1.019–1.035` at 2 threads and
`0.987–1.069` at 16 threads for default, low-T, high-T, and stiff cases. PCHIP
is the current default; Simpson remains explicitly selectable. The Oracle C
relative errors for PCHIP are `1.07e-5` (default), `1.66e-4` (low-T),
`9.88e-6` (high-T), and `1.15e-5` (stiff).

## Dated measurements and comparability

The 2026-09-10 fixed-resource benchmark recorded `4.26 ms/point` warm median
at 20 threads. It is a dated snapshot and is not directly comparable with the
two-thread Round 68 or 2026-10-21 profiles. The 2026-09-03 before/after report
is also a historical comparison; its protocol and limitations are in
[`performance_comparison_20260903.md`](performance_comparison_20260903.md).

Later `profile_fast_breakdown_*` artifacts include diagnostics through
2026-10-21. Earlier runs use different thread/resource settings and profiling
protocols, so they do not replace a like-for-like benchmark or establish a
cross-version speed claim. The separate 2026-09-24 MCMC report measures
sampling-step time and effective samples per second, not isolated solver
runtime; its per-scenario speed ratios are about 2.9–4.9x and are reported in
[`mcmc_sagenet_compare/report.md`](mcmc_sagenet_compare/report.md). Browse the
solver artifacts by family in
[`experiment_catalog.md`](experiment_catalog.md).

## `DN_re` edge accuracy and latency (2026-10-22)

The full-grid Oracle C follow-up found that the standard fast profile
(`h=0.005`) differs from the WKB-corrected reference by `2.74e-3` and
`2.49e-3` at `edge_dnre_lo` and `edge_dnre_hi`. A matched Oracle rerun at
`h=0.0003125` brings the differences to `1.65e-4` and `1.63e-4`, below the
scoped `2e-4` gate. The more accurate candidate is much slower:

| Point | Standard (`h=0.005`) | Fine uniform (`h=0.0003125`) | Transition-refine (`h=0.0025`) | Fine uniform vs WKB | Transition-refine vs WKB |
|---|---:|---:|---:|---:|---:|
| `edge_dnre_lo` | `5.69 ms` | `62.7 ms` | `46.5 ms` | `1.65e-4` | `1.08e-6` |
| `edge_dnre_hi` | `4.53 ms` | `71.7 ms` | `73.7 ms` | `1.63e-4` | `1.59e-4` |

Each configuration uses two warmup calls and ten timed calls, with a fresh
model per solve, two Numba `workqueue` threads, and one BLAS thread. The
comparison tool verifies exact agreement in `DN_eff` and the native frequency
grid before calculating error; configurations with mismatched anchors fail
closed. These results cover two points only, do not establish a general speed
ratio, and have not changed the formal fast preset. They were measured in
separate invocations with the same declared thread budget, without CPU affinity
pinned; use the raw samples when assessing host timing noise. Raw anchors, comparisons,
and reproduction commands are in
[`oracle_c_wkb_assessment.md`](oracle_c_wkb_assessment.md) and
[`oracle_c_wkb_corner_20261022/`](oracle_c_wkb_corner_20261022/).

The transition-refined candidate uses the vectorized H2 implementation and
has errors of `3.97e-6` / `1.64e-4` versus deep `z=10`; it is about `8.2x` /
`16.3x` slower than standard fast. It beats the uniform fine-step runtime on
`edge_dnre_lo` only. Its matched Oracle records are in
[`oracle_c_wkb_corner_20261022/transition_h00025/`](oracle_c_wkb_corner_20261022/transition_h00025/).

Exposing the refinement radius and factor in `FastSolverConfig` made a smaller
zone reproducible. At `h=0.0025`, radius `0.25`, and factor `8`, the matched
edge errors are `5.17e-6`/`1.62e-4` versus deep `z=10`, with warm medians of
`35.6`/`51.7 ms`. This is `1.29x`/`1.42x` faster than radius `2.0`, but remains
about `6.3x`/`11.4x` slower than standard fast. The preset is unchanged. Raw
Oracle and repeated fast timing records are in
[`oracle_c_wkb_corner_20261022/transition_h00025_r025/`](oracle_c_wkb_corner_20261022/transition_h00025_r025/).

| Base `h` | Refine radius/factor | Local step | Low/high error vs deep | Low/high median | Gate |
|---:|---:|---:|---:|---:|---|
| `0.0025` | `0.125 / 8` | `0.0003125` | `5.61e-6 / 1.62e-4` | `33.9 / 49.8 ms` | pass |
| `0.005` | `0.125 / 16` | `0.0003125` | `1.85e-5 / 3.84e-6` | `19.0 / 26.7 ms` | pass |
| `0.02` | `0.125 / 64` | `0.0003125` | `1.16e-4 / 8.25e-5` | `7.63 / 9.49 ms` | pass |
| `0.04` | `0.125 / 64` | `0.000625` | `4.49e-4 / 3.17e-4` | `4.99 / 6.05 ms` | fail |

The `h=0.02` point is the fastest passing candidate measured so far. It remains
about `1.34x/2.09x` slower than standard fast, so no production promotion is
made. The `h=0.04` result shows that matching only the window radius is not
enough: doubling its local step loses the accuracy gate. Per-run configurations,
reference residuals, exact anchor checks, and samples are linked from
[`oracle_c_wkb_assessment.md`](oracle_c_wkb_assessment.md).

The `h=0.04`, radius-`0.25`, factor-`128` comparison keeps the local step at
`0.0003125` but has strongly different edge costs: repeated warm medians are
`47.15 ms` and `7.75 ms` (two warmups and ten measured calls). Its low-edge
error remains above the gate, while the high edge passes. Treat this as a
rejected two-point candidate, not a general speed result; the exact samples
and environment are in the linked raw JSON.

At `h=0.02/factor=64`, reducing the radius from `0.125` to `0.0625` lowers the
two-edge medians from `7.63/9.49 ms` to `6.49/8.74 ms` and passes both accuracy
gates with slightly smaller errors. This is a matched two-point result for the
validation path; the formal fast preset is still faster and unchanged.
Reducing the radius again to `0.03125` changes the medians only slightly, to
`6.41/8.09 ms`, in a separate unpinned run. Treat this as a candidate result,
not a stable speedup over radius `0.0625`.

A same-process comparison pinned to CPUs `[0, 1]` alternated radii `0.015625`
and `0.0078125` for 30 measured pairs per edge. Candidate/base median ratios
were `0.989`/`0.979`, while bootstrap 95% intervals were
`0.977–1.013`/`0.936–1.014`. Both intervals include `1`, so a stable speed gain
is not established. Raw paired samples are in
[`benchmark_transition_radius_paired_20261022.json`](benchmark_transition_radius_paired_20261022.json).

## Transition-refined `H2_vec` batching (2026-10-22)

`H2_vec` used to call the SciPy density spline once per background node while
building a transition-refined grid. Batching the same spline queries preserves
the old scalar results exactly and reduces full-solver latency on that path:

| Point | Scalar-reference median | Batched median | Speedup | Output arrays |
|---|---:|---:|---:|---|
| `edge_dnre_lo` | `93.4 ms` | `22.8 ms` | `4.09x` | exact match |
| `edge_dnre_hi` | `144.1 ms` | `37.3 ms` | `3.86x` | exact match |

The A/B used two warmup calls and ten alternating measured calls per method,
fresh models, two Numba `workqueue` threads, one BLAS thread, and unpinned CPU
affinity. It covers the validation `transition-refine` path (`h=0.005`, goal grid, PCHIP); it does not
change the formal exact-kink fast path or add new Oracle accuracy coverage.
The scalar implementation is retained in the benchmark script as the
pre-vectorization reference. Raw samples and the reproduction driver are
[`benchmark_h2_vectorization_20261022.json`](benchmark_h2_vectorization_20261022.json)
and [`benchmark_h2_vectorization.py`](../scripts/benchmark_h2_vectorization.py).
