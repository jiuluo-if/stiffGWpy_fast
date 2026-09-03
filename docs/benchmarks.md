# Benchmarks

Status: current
Date: 2026-09-03
Code version: see manifest `commit`

Runtime (this host, warm, 4 threads): plain-grid ≈0.37 s/point, production
(matched z8) ≈3.7–4.1 s/point, reference (oracle) ≈360 s/point.  The LSODA
anchor is ≈18.6 s/point (z5), so production is ≈4.5x faster at matched
accuracy; plain-grid is ≈50x faster but is **not** accuracy-certified
(signal rel median 1.9e-2).

Runtime-vs-physical-error Pareto: `docs/archive/reference/pareto_default.json`,
`deep_oracle_default.json` (superseded / historical).  The pre-fix warm
`~1000x` plain-grid headline is archived and is not a current claim — it holds
only for the coarsest pre-fix grid at p95 ~7e-2 spectrum error.
