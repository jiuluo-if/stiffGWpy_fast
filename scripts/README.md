# Repository scripts

[中文说明](README_zh.md) | English

This directory contains validation, benchmark, profiling, packaging, and
historical diagnostic scripts. Run commands from the repository root so that
package imports and artifact paths resolve as documented.

## Script groups

| Prefix or file | Purpose |
|---|---|
| `validate_*.py` | Replay fast/reference comparisons, parameter edges, sweeps, and mode checks |
| `build_*.py` | Assemble validation matrices and manifests from existing artifacts |
| `bench_*.py`, `benchmark_*.py` | Measure runtime, accuracy, or a named diagnostic hypothesis |
| `profile_*.py` | Attribute end-to-end runtime to solver stages |
| `smoke_*.py`, `verify_distribution.py` | Check installed-wheel resources and package boundaries |
| `plot_*.py` | Render plots from saved validation or benchmark data |
| `_resource_budget.py`, `_smoke_probe.py` | Shared support for thread limits and isolated smoke checks |

Many `benchmark_*_spike.py` files are historical, standalone candidate tests;
they do not select or modify a production solver path. Read a script's CLI
arguments (or run `--help` when available) and check its output path first.
Preserve dated artifacts and
use a new output name for a new experiment unless replacement is intentional.
The SageNet comparison driver also requires a separate SageNet checkout and
runtime; SageNet is not installed by this package's extras.

## Reproducible runs

Use the documented environment, thread budget, seeds, warmup policy, and repeat
count for the question being tested. Separate cold/JIT startup from warm timing.
Record the source commit and keep numerical accuracy comparisons separate from
speed comparisons. See [`docs/reproducibility.md`](../docs/reproducibility.md),
[`docs/benchmarks.md`](../docs/benchmarks.md), and
[`docs/experiment_catalog.md`](../docs/experiment_catalog.md).
