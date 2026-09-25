# stiffgwpy_fast — documentation index

[中文说明](README_zh.md) | English

Status: current documentation map
Date: 2026-09-25
Code version: each evidence artifact records its own commit; the manifest is dated `2026-09-17` and is not current-HEAD certification.

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

The root README is the top-level user guide. `fast` is the only user-facing
profile; historical `plain-grid` and `production` results apply only to their
recorded artifacts. Read accuracy values with each run's commit and scope: the
manifest is a read-only summary of earlier validations, not a current-HEAD
certification.

The posterior mock's independent reference truth is stored in
`mcmc_posterior/oracle_truth.json`. It is required test data, not a temporary
benchmark output.

| Document | Content | Status |
|---|---|---|
| `physics.md` | background + tensor-mode + `Delta N_eff` closure | current |
| `numerical_method.md` | single fast profile, exact kink split, tail, goal grid | current |
| `accuracy.md` | layered accuracy limits and dated validation evidence | see manifest and experiment catalogue |
| `parameter_validation.md` | parameter schema and date-bound sweep records | see manifest and experiment catalogue |
| `cobaya.md` | adapter options, mode mapping, `eval_freqs` | current |
| `benchmarks.md` | current scoped audit plus date-bound benchmark records | see report scope |
| `performance_comparison_20260903.md` | optimization before/after, breakdown, AB and gates | dated historical snapshot |
| `reproducibility.md` | drivers, gates, environment metadata | current |
| `experiment_catalog.md` | research themes, evidence status, raw experiment families | current |
| `archive/` | superseded snapshots and local-only diagnostic logs | historical |
| `mcmc_sagenet_compare/` | 2026-09-24 LVK MCMC speed, diagnostics, and posterior-spectrum comparison | dated experiment |
| `mcmc_posterior/` | importance-sampling posterior validation against the continuous-sigma reference | dated validation |

The 2026-09-11 `baseline_54d65e3.md` snapshot is retained at
[`archive/baselines/baseline_54d65e3.md`](archive/baselines/baseline_54d65e3.md).
It records an earlier two-thread baseline and is not a current release claim.
