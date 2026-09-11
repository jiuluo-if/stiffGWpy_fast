# stiffGWpy_fast

[中文说明](README_zh.md) | English

The Chinese overview is maintained in [README_zh.md](README_zh.md); keep the
two language entry pages synchronized when updating features, parameters, or
validation conclusions.

The PyPI distribution name is `stiffgwpy_fast`; the Python import name is
`stiffgwpy_fast`. They are intentionally distinct from the older `stiffgwpy`
project.

**LCDM + stiff matter + primordial stochastic gravitational-wave background (SGWB)**
cosmology code, with one **user-facing fast profile**, an independent
continuous-sigma high-accuracy reference pipeline
(the precision oracle), and the original LSODA path kept only for regression and
runtime benchmarking.

> **Documentation rule.**  Accuracy claims are read back from
> `docs/validation/validation_manifest.json`; current performance measurements
> are maintained in `docs/benchmarks.md` and
> `docs/performance_comparison_20260903.md`.  `docs/` holds the canonical
> per-topic documents; superseded historical audits remain available in Git
> history rather than in the active `docs/` tree.

---

## What it computes

For a given cosmology (`LCDM_SG`), `stiffGWpy`:

1. builds the background (radiation + neutrinos + stiff matter + Lambda),
2. solves the tensor-mode perturbation equation per frequency channel, and
3. iterates `Delta N_eff` (the SGWB contribution to extra radiation) to
   self-consistency.

The outputs are the today spectrum `Omega_GW(f)`, the integrated
`Delta N_eff` (`DN_gw`), the number of radiation species `kappa_r`, the
expansion `g2`/`w2` curves, and derived quantities used by the Cobaya theory.

## Physics model

Parameters (all enter the background / tensor-source physics):

| Parameter | Meaning |
|---|---|
| `Omega_bh2`, `Omega_ch2`, `H0` | baryon/CDM densities and Hubble constant |
| `DN_eff` | constant extra radiation (Delta N_eff) |
| `A_s`, `r`, `n_t` | scalar amplitude, tensor-to-scalar ratio, tensor tilt |
| `cr` | >0: enforce the single-field consistency relation; <=0: use `n_t`/`DN_re` |
| `T_re` | reheating temperature [GeV] |
| `DN_re` | e-folds of matter-like reheating |
| `kappa10` | `rho_stiff / rho_photon` at 10 MeV |

The instantaneous-reheating transition is a kink in `sigma(N)`.  The formal
`fast` profile now combines an exact `N_re` transfer split, phase-capped
sub-stepping, and a goal-oriented sparse frequency grid.  Historical
`production`/`transition_refine` names remain only as validation compatibility
entries; they are not separate user-facing production tiers.  The precision
anchor remains the independent `stiffgwpy_fast.reference` pipeline.

## Fast profile

The current branch's single `fast` preset is:

| Config key | `fast` |
|---|---|
| `h` | 0.005 |
| `col_step` | 8 |
| `z_tail` | 5.0 |
| `phase_max` | 0.25 |
| `freq_grid` | goal (typically 70–120 points) |
| `kink_split` | on |
| `frequency_quadrature` | pchip (default; `simpson` stays explicitly selectable) |
| outer tol | 1e-6 |

The goal grid reserves nodes around the reheating feature and preserves
`eval_freqs` as native solve nodes.  The deep-subhorizon part is handed to the
analytic WKB envelope at `z_tail`; the local error budget is exposed through
`stiffgwpy_fast.fast_sgwb.estimate_local_error`.

Fresh HEAD evidence uses the formal `fast` path (`h=.005`, `col_step=8`,
`z_tail=5`, `phase_max=.25`, exact kink split, goal grid). On the six-point
A,B,B,A paired matrix at fixed 20-thread `workqueue`, the default warm
median/p95 is `4.77/5.58 ms/point`; this does not yet meet the `<4 ms` target. On the
same native 76-node grid, the independent reference comparison gives a
`DN_gw` relative error of `2.94e-4`, which is a combined tail/transfer and
frequency-quadrature residual; against the independent Oracle C WKB anchor the
PCHIP quadrature residual alone is `1.07e-5`–`1.66e-4` at the four named
points. PCHIP is therefore the single formal `fast` default, and `simpson`
stays explicitly selectable for audits. `eval_freqs` is separated from
integration support nodes, and its six-point DN invariant measured zero change.
See `findings.md` and `progress.md` for exact artifacts and rejected
candidates.

## Reference / oracle

`stiffgwpy_fast/reference.py` is an independent, higher-order implementation of the
same physics: continuous `sigma(N)` (kink as an exact breakpoint), adaptive
`DOP853` per frequency mode, shape-preserving PCHIP + adaptive Gauss-Kronrod
quadrature with error estimates.  It is **the** accuracy anchor.  LSODA is used
only for regression checks and runtime benchmarking, never as a precision
truth value.

**Oracle independence audit (honest).**  The reference also uses a frozen
analytic tail at `z_tail`.  Measuring its own sensitivity at the default point:

| oracle choice | `DN_gw` relative change |
|---|---|
| `z_tail` 7 → 8 | 4.2e-4 |
| `z_tail` 8 → 10 | 3.0e-4 |
| `z_tail` 14 (deep/no-tail) | *infeasible* — the ODE becomes deep-subhorizon stiff |

So the reference itself carries a ~3e-4 `z_tail`-frozen-tail sensitivity, and the
production engine's ~4e-4 residual is at the same level.  This is reported as an
honest limit, not hidden behind a gate.  `reference.oracle_variants()` runs
oracle A/B/C and reports `CONSISTENT` / `ORACLE-SENSITIVE`.

## Installation

```bash
pip install stiffgwpy-fast          # 从 PyPI 安装最新版本
pip install .                 # 安装运行依赖：numpy scipy astropy pyyaml numba
pip install .[cobaya]         # 安装串行 Cobaya MCMC 接口
pip install .[cobaya,mpi]     # 安装 Cobaya 和可选 mpi4py
pip install .[dev]            # 安装 pytest、ruff、build、matplotlib
```

Packaging and PyPI publishing instructions are in
[PACKAGING.md](PACKAGING.md). The release archive intentionally excludes
`docs/`, tests, validation scripts, CI files, and research-only configuration.

## Local privacy and repository hygiene

The repository ignores local credentials and machine-generated files, including
`.pypirc`, `.env*`, key/certificate files, caches, coverage reports, build
directories, MCMC chains, checkpoints, and temporary parameter-sweep outputs.
Keep the local PyPI token outside version control and never paste it into an
issue, README, command log, or commit. Before committing, check
`git status --short` and use `git check-ignore -v <path>` for any local file
that may contain credentials or private experiment output.

The committed `docs/` tree contains only selected validation documentation and
small reproducibility summaries. Reproducible but large or unfinished local
research outputs are intentionally excluded from Git and from PyPI archives.

## Python usage

```python
from stiffgwpy_fast import LCDM_SG

m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
m.SGWB_iter()  # 默认使用唯一 fast goal-kink-hybrid preset
print(m.DN_gw[-1])

# 显式选择唯一正式 fast 档位
m2 = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
m2.SGWB_iter(engine='fast', accuracy_mode='fast')

# 显式选择原始 LSODA 回归路径
m3 = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
m3.SGWB_iter(engine='lsoda')
```

`accuracy_mode='fast'` is the only formal user-facing fast profile.
`production`/`transition_refine` are deprecated compatibility aliases and are
mapped to `fast` by the high-level API; direct validation code may retain its
internal presets. Explicit `h`/`col_step`/`z_tail`/`freq_res`/`tol` override
the selected preset. For the high-level API, `SGWB_iter()` defaults to the
`fast` engine and its goal-kink-hybrid preset. Pass `accuracy_mode=None` only when you
intentionally want a snapshot of the legacy manual module settings. Lower-level
calls can pass an immutable `fast_sgwb.FastSolverConfig` per invocation.

```python
from stiffgwpy_fast import fast_sgwb
b = fast_sgwb.estimate_local_error(m)   # point-local 11-category budget
print(b['DN_gw_error'], b['certification_status'])
```

## Cobaya usage

```yaml
theory:
  stiffgwpy_fast.cobaya.stiffGW.stiffGW:
    engine: fast
    fallback: True
    accuracy_mode: fast            # 唯一正式 fast 档位
    fast_threads: 8
```

**Mode mapping is unambiguous.**  In `stiffGW.yaml` the fast knobs default to
`0` (a sentinel meaning "use the selected accuracy_mode").  The resolution order
is strictly:

```
accuracy_mode  ->  preset defaults  ->  explicit user overrides only
```

so a default `engine: fast` YAML runs the combined `fast` settings.  The YAML
default values cannot silently mask the selected preset. `accuracy_mode: fast`
maps to the goal-kink-hybrid path; `accuracy_mode: production` is retained only
for validation compatibility.

**`eval_freqs` (likelihood bins as native nodes).**  Set
`eval_freqs: [list|path-to-file]` in the theory YAML to force-add
`log10(f/Hz)` native solve nodes.  This is the path by which likelihood
frequency bins reach the fast solver as native nodes
(`SGWB_iter_fast(..., eval_freqs=...)`), removing the interpolation error at
steep spectral features (Layer C measured the per-bin dex interpolation error at
≤3.1e-4 across 11 PTA-like bins).  By default (`eval_freqs: null`) the solver
uses its own grid and the likelihood interpolates over the returned spectrum;
use native nodes when the bin spacing approaches the spectral features. Native
evaluation nodes are excluded from the bolometric integration support grid, so
changing `eval_freqs` does not change self-consistent `DN_gw`.

## Accuracy

The current single fast path is still under active branch-level optimization.
On the exact native 76-node default grid, the independent continuous-sigma
oracle gives a `DN_gw` relative error of `2.94e-4` for the default PCHIP
quadrature (the legacy Simpson option is `7.14e-4`), while the independent
Oracle C WKB anchor bounds the quadrature-only residual at `1.07e-5`–`1.66e-4`.
The 76/80/90/110-node sweep is non-monotonic, so these numbers are not a
universal parameter-space certification. `reference` remains the precision
oracle.

## Full parameter validation

Parameter schema (11 physical params, ranges in `scripts/validate_two_modes.py`):

* **Single-parameter axis edges:** `docs/paramsweep_z8b/` (16 points on
  r/n_t/cr/T_re/DN_re/kappa10 axis edges + transition interiors).
* **Param space (LHS, plain-grid screen):** 400 points,
  **254 success / 146 shared-`Delta_Neff` guard / 0 numerical failure**
  (`docs/validation/param_sweep_plain.json`).  The 36% guard fraction is a
  physical rejection (total `N_eff > 5`), reported explicitly, never hidden.
* **Param space (Sobol, production):** 240 points, **212 ok / 28 guard**
  (`docs/paramsweep_ref/fast_sweep.jsonl`).

Rejections are categorised `PHYSICAL_INVALID` / `PHYSICAL_GUARD` /
`NUMERICAL_FAILURE` / `FAST_ERROR` / `ORACLE_ERROR` (see the manifest and
`scripts/validate_two_modes.py`).  A guard rejection is a solver-honest physical
rejection, not a numerical failure.

## Benchmark

Current measurements are summarized below; cold JIT and warm runtime are kept
separate.  See `docs/performance_comparison_20260903.md` for the full method
comparison, stage breakdown, AB evidence, and thread scaling.

| | runtime/point | vs LSODA |
|---|---|---|
| fast (goal-kink-hybrid) | 4.77 ms warm median; 5.58 ms p95; 0.222 s cold | interim result; `<4 ms` not yet met |
| reference (oracle) | ≈360–383 s/point historical | anchor only |

The speedup entries use the recent A-point LSODA measurement (`22.137 s`) and
candidate warm median; older `0.37 s` / `3.7–4.1 s` / `~1000x` figures are
pre-JIT historical measurements and are not current claims.

## MCMC validation

Importance-sampling posterior validation (Layer C) from 9000 fast-production
draws with a fixed seed: ESS **4167** (gate 2000 PASS), `log10 r` posterior shift
**-0.0011 sigma** (gate <0.1 sigma PASS), per-bin dex max **3.1e-4**,
`|Delta logL|` max **7.3e-3** (gate 0.1 PASS).  See
`docs/mcmc_posterior/posterior_validation.md`.

Honest limit: a full two-chain reference-engine MCMC was **not** run because a
reference solve is ≈360 s/point on this host.  The former bounded real-Cobaya
scaffold chains were unfinished adapter-plumbing outputs and have been removed
from the working tree.  Posterior-shift conclusions therefore rest on
importance reweighting, not on an independent reference chain.

## Limitations

* The current same-grid default `DN_gw` error is `2.94e-4` with PCHIP (the
  legacy Simpson option is `7.14e-4`); that residual is a combined tail and
  frequency-quadrature term, so the branch has not yet met the final `<2e-4`
  DN gate or the first `<=4 ms/point` speed target.
* Fast execution is over 100x faster than the recent LSODA A-point runtime,
  but this is an interim optimization result, not an accuracy certification.
* MCMC validation rests on importance reweighting, not an independent
  reference chain (reference is ~360 s/point).
* The oracle uses a frozen tail; a deep/no-tail certification (z_tail ≥ 14) is
  computationally infeasible because the mode equation becomes deep-subhorizon
  stiff.

## Reproducibility

The contributor/release gate is summarized in [CONTRIBUTING.md](CONTRIBUTING.md)
and the staged findings are tracked in [docs/engineering_audit.md](docs/engineering_audit.md).

```bash
python scripts/validate_two_modes.py --phase convergence          # fast-vs-fast convergence
python scripts/validate_two_modes.py --phase param_sweep --n 400  # LHS screen
python scripts/build_two_mode_manifest.py                          # manifest (read-only replay)
python -m pytest                                                   # tests (slow deselected)
python -m pytest -m cobaya                                         # Cobaya adapter gate
python scripts/validate_manifest.py                                # committed artifact schema
python -m build --wheel                                            # wheel build gate
python scripts/smoke_installed_wheel.py dist/stiffgwpy_fast-*.whl        # installed-resource smoke
```

Every driver records git-commit + environment metadata.  The regression suite is
103 passed (6 slow LSODA gates deselected by default; counts may change as tests evolve).

## Directory structure

```
stiffgwpy_fast/            pip package
  stiff_SGWB.py       LCDM_SG + engine dispatch (fast|lsoda|reference)
  fast_sgwb.py        fast solver, ACCURACY_MODES, FAST_PROFILES, estimate_local_error
  reference.py        independent continuous-sigma oracle (engine='reference')
  freq_adaptive.py    curvature-adaptive frequency grid
  exact_background.py continuous-sigma expansion integrals / kink-refined grid
  cobaya/             Cobaya theory adapter + likelihoods
tests/                pytest suite (99 unit + 6 slow gates)
scripts/              validation drivers (validate_two_modes, build_two_mode_manifest, ...)
docs/                 selected docs + compact validation artifacts; local chains/outputs are ignored
Git history           superseded historical audits / benchmarks
```

## License

GPL-3.0 (see `LICENSE.md`).
