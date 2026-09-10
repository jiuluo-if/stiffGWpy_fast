# Benchmarks

Status: current
Date: 2026-09-10
Code version: see manifest `commit`

Current formal fast benchmark (Windows, fixed CPU 0-19, Numba
`workqueue`, 20 threads, 25 repeats) reports default `5.36 ms/point` warm
median (`5.73 ms` p95; cold JIT `0.254 s`). The six-point medians are
`4.89/7.28/7.54/4.59/7.28 ms` for low-T/high-T/stiff/low-r/high-kappa.
These are separate cold/warm measurements and do not yet satisfy the `<4 ms`
acceptance target.

The current default same-grid independent reference comparison reports
Simpson/PCHIP `DN_gw` relative errors `7.14e-4/2.94e-4`; PCHIP is therefore
kept opt-in. The 76/80/90/110 sweep is non-monotonic, so no grid promotion is
claimed. See the HEAD-stamped JSON artifacts and `findings.md` for details.

The detailed before/after comparison, breakdown, thread scaling and numerical AB
are in `docs/performance_comparison_20260903.md`. The plain-grid oracle envelope
is unchanged: signal relative median `1.867e-2`, max `7.019e-2`.

Older runtime-vs-physical-error Pareto artifacts were removed from the active
tree; their old `~1000x` headline is not a current claim and remains traceable
through Git history.
