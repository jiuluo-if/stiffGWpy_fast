# Benchmarks

Status: current
Date: 2026-09-03
Code version: see manifest `commit`

Current benchmark (Windows, `FAST_THREADS=4`): plain-grid default A is
`4.442 ms/point` warm median (`5.105 ms` p95; cold JIT `0.325 s`), and
production is `21.772 ms/point` warm median (`22.149 ms` p95; cold JIT
`0.226 s`). A recent LSODA A-point run is `22.137 s`; the independent reference
remains a historical `~360–383 s/point` anchor. These are separate cold/warm
measurements, not a claim that the physical accuracy tiers changed.

The detailed before/after comparison, breakdown, thread scaling and numerical AB
are in `docs/performance_comparison_20260903.md`. The plain-grid oracle envelope
is unchanged: signal relative median `1.867e-2`, max `7.019e-2`.

Older runtime-vs-physical-error Pareto artifacts were removed from the active
tree; their old `~1000x` headline is not a current claim and remains traceable
through Git history.
