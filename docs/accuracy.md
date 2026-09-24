# Accuracy and evidence boundaries

Status: current summary of scoped evidence; no full parameter-space release certification

Date: 2026-09-24

The formal fast profile is the single `goal-kink-hybrid` mode. The independent
continuous-sigma DOP853 pipeline is the precision reference; LSODA is a
regression/runtime path, not a truth oracle. Every numerical result below is
bounded by the artifact's commit, reference convention, parameter points, and
frequency grid.

## Latest scoped oracle evidence

The latest detailed accuracy audit is the dated
[`fast_v02_audit_report.md`](fast_v02_audit_report.md), generated from its
recorded 2026-09-12 evidence. It is not a full parameter-space certification or
a fresh precision rerun at the current repository HEAD.

| Comparison | Recorded result | Scope and interpretation |
|---|---:|---|
| PCHIP fast vs Oracle C WKB `DN_gw` | `1.07e-5`–`1.66e-4` | Four named points; below the `2e-4` gate on those points |
| Fast vs reference, same native grid, reference `z_tail=8` | median `2.93e-4`, max `2.96e-4` | Six named points; includes the reference frozen-tail convention |
| Fast vs reference, default point, reference `z_tail=10` | `4.96e-6` | One-point deep-tail attribution, not a parameter-space result |
| Nested native-frequency refinement | relative change `7.4e-12`–`2.8e-9` | Six named points; supports grid convergence in that test |

The same audit reports Simpson errors above the PCHIP results and records
PCHIP as the formal default. Its `z_tail=8` reference difference is limited by
the reference's own frozen-handoff tail sensitivity; the one-point `z_tail=10`
comparison is an attribution check, not a substitute for broader validation.

## Latest runtime evidence is separate

The 2026-09-23 Round 68 full-path profile used 25 repeats, two Numba
`workqueue` threads, one BLAS thread, and CPU affinity `[0, 1]`. All six
profiled cases converged; the default warm median/p95 was `7.189/7.916 ms`.
This profile did not run a matched-resource LSODA comparison and did not add a
new oracle accuracy sweep. It therefore establishes neither a speedup factor
nor new accuracy coverage. Full per-case values and provenance are in
[`benchmarks.md`](benchmarks.md).

## Likelihood and posterior evidence

The 2026-09-03 importance-sampling validation compares the fast spectrum with
the continuous-sigma reference on 240 posterior-bulk points and 11 likelihood
bins. It reports max per-bin dex error `3.10e-4`, max absolute `Delta logL`
`7.30e-3`, ESS `4167`, and a `log10(r)` posterior shift of `-0.0011 sigma`.
These are results of the recorded mock-data/importance-reweighting procedure,
not an independent reference-engine MCMC chain; see
[`mcmc_posterior/posterior_validation.md`](mcmc_posterior/posterior_validation.md).

The separate 2026-09-24 LVK report runs four chains for each method and
scenario. It contains convergence diagnostics and posterior-spectrum samples;
SageNet+ does not meet the report's diagnostic references in two scenarios.
The low-reheating scenario has 0% LVK-band coverage for all methods and cannot
support a data-fit conclusion. See
[`mcmc_sagenet_compare/report.md`](mcmc_sagenet_compare/report.md).

## Local error budget

`estimate_local_error` distinguishes telemetry measured during the solve,
terms calibrated from a fiducial point, and uncertified defaults when no solve
telemetry is attached. A `certification_status` of
`certified-fiducial-calibrated` is limited to that calibration model; it does
not imply a uniform error guarantee over all parameter values.

## Accepted execution optimization evidence

The dated H2 endpoint-cache A/B recorded identical spectrum and background
output digests over its 13 named/Sobol/edge points and 50 repeats. The measured
high-kappa warm median improved by 6.54%; the low-T case regressed within the
reported noise range. This supports that specific cache change and protocol,
not a universal runtime improvement. Consult the linked artifact before
reusing the result.
