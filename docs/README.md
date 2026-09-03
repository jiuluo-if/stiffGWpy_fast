# stiffgwpy_fast — documentation index

Status: current for this revision
Date: 2026-09-03
Code version: see `docs/validation/validation_manifest.json` → `commit`

The README is the top-level user guide.  These documents give the substance
behind the two user-facing fast profiles.  Every accuracy number is read back
from `docs/validation/validation_manifest.json`, which is itself generated from
the committed validation artifacts (`docs/paramsweep_*`, `docs/mcmc_posterior/`)
by `scripts/build_two_mode_manifest.py` (a read-only replay — no physics re-run).

| Document | Content | Status |
|---|---|---|
| `physics.md` | background + tensor-mode + `Delta N_eff` closure | current |
| `numerical_method.md` | fast two-profile scheme, transition refine, tail, adaptive grid | current |
| `accuracy.md` | layered accuracy, honestly reported limits | current |
| `parameter_validation.md` | parameter schema, axis + space sweep, classification | current |
| `cobaya.md` | adapter options, mode mapping, `eval_freqs` | current |
| `benchmarks.md` | current cold/warm runtime and method comparison pointer | current |
| `performance_comparison_20260903.md` | optimization before/after, breakdown, AB and gates | current |
| `reproducibility.md` | drivers, gates, environment metadata | current |
| Git history | superseded historical audits / benchmarks (LSODA-era, pre-fix) | archived |
