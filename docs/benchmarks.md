# Benchmarks

Status: current summary of the audited scope; results remain commit- and protocol-bound

Date: 2026-09-24

Code version: see `docs/validation/validation_manifest.json` → `commit` and each run artifact's `commit`

## Latest fixed-resource full-path profile

The only user-facing fast profile uses `frequency_quadrature=pchip` with
`h=0.005`, `col_step=8`, `z_tail=5`, `phase_max=0.25`, exact kink split, and
the goal frequency grid. The newest full-path profile is the 2026-09-23 Round
68 run at solver commit `c7d4766f15c1cda4693d306597723de50abd4554`; the
current repository HEAD adds MCMC/reporting work but does not change
`stiffgwpy_fast/` solver source after that commit.

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

## Scoped accuracy audit (dated 2026-09-12)

The latest detailed PCHIP/oracle audit is summarized in
[`fast_v02_audit_report.md`](fast_v02_audit_report.md), generated from its
recorded 2026-09-12 evidence rather than re-run at the Round 68 commit:

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
Round 68 profile above, which used two threads. The 2026-09-03 before/after
report is also a historical comparison; its protocol and limitations are in
[`performance_comparison_20260903.md`](performance_comparison_20260903.md).

Later `profile_fast_breakdown_*` artifacts include diagnostics through
2026-09-23. They use different thread/resource settings and profiling
protocols, so they do not replace a like-for-like benchmark or establish a
cross-version speed claim. The separate 2026-09-24 MCMC report measures
sampling-step time and effective samples per second, not isolated solver
runtime; its per-scenario speed ratios are about 2.9–4.9x and are reported in
[`mcmc_sagenet_compare/report.md`](mcmc_sagenet_compare/report.md). Browse the
solver artifacts by family in
[`experiment_catalog.md`](experiment_catalog.md).
