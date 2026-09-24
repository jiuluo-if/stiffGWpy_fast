# stiffgwpy_fast — documentation index

[中文说明](README_zh.md) | English

Status: current documentation map
Date: 2026-09-24
Code version: see `docs/validation/validation_manifest.json` → `commit`

## Research and experiment navigation

Start with [`experiment_catalog.md`](experiment_catalog.md) for the grouped
research threads, their current evidence status, and links to the supporting
reports and raw artifacts. It separates accepted findings from rejected
prototypes, exploratory runs, and still-unverified release claims.

Raw JSON/JSONL files are run records. Their filename families and directories
are grouped in the catalogue so that related rounds can be compared without
mistaking a smoke run or a candidate profile for a production result. Keep the
recorded run paths intact: experiment output files and validation manifests
refer to those paths.

The README is the top-level user guide.  These documents give the substance
behind the two user-facing fast profiles.  Every accuracy number is read back
from `docs/validation/validation_manifest.json`, which is itself generated from
the committed validation artifacts (`docs/paramsweep_*`, `docs/mcmc_posterior/`)
by `scripts/build_two_mode_manifest.py` (a read-only replay — no physics re-run).

The independent reference truth used by the posterior mock is kept at
`mcmc_posterior/oracle_truth.json`; it is required test data, not a temporary
benchmark output.

| Document | Content | Status |
|---|---|---|
| `physics.md` | background + tensor-mode + `Delta N_eff` closure | current |
| `numerical_method.md` | fast two-profile scheme, transition refine, tail, adaptive grid | current |
| `accuracy.md` | layered accuracy limits and dated validation evidence | see manifest and experiment catalogue |
| `parameter_validation.md` | parameter schema and date-bound sweep records | see manifest and experiment catalogue |
| `cobaya.md` | adapter options, mode mapping, `eval_freqs` | current |
| `benchmarks.md` | current scoped audit plus date-bound benchmark records | see report scope |
| `performance_comparison_20260903.md` | optimization before/after, breakdown, AB and gates | current |
| `reproducibility.md` | drivers, gates, environment metadata | current |
| `experiment_catalog.md` | research themes, evidence status, raw experiment families | current |
| `archive/` | superseded snapshots and local-only diagnostic logs | historical |

The 2026-09-11 `baseline_54d65e3.md` snapshot is retained at
[`archive/baselines/baseline_54d65e3.md`](archive/baselines/baseline_54d65e3.md).
It records an earlier two-thread baseline and is not a current release claim.
