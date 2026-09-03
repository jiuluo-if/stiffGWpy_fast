# stiffGWpy_fast

**LCDM + stiff matter + primordial stochastic gravitational-wave background (SGWB)**
cosmology code, with a fast solver exposed as exactly **two user-facing fast
profiles**, an independent continuous-sigma high-accuracy reference pipeline
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

The instantaneous-reheating transition is a kink in `sigma(N)`.  The fast
production profile treats that kink as an ODE integration breakpoint so it is
never smeared by a grid spline; the plain-grid profile does not, which is the
main source of its bias (see below).

## Two fast modes

There are exactly **two user-facing fast profiles**.  The extra names
(`debug`, `deep`, `reference`) are validation/benchmark variants of the same
internal solver — they are not advertised as additional production tiers and are
not meant for the MCMC thermal path.  The true precision anchor is the
independent `stiffgwpy.reference` pipeline.

| Config key | `fast` (plain-grid) | `production` (transition-refine) |
|---|---|---|
| `h` | 0.02 | 0.01 |
| `col_step` | 8 | 4 |
| `z_tail` | 5.0 | 8.0 |
| `freq_res` | 1.0 | 1.0 |
| `transition_refine` | off | on |
| `phase_max` | 0.0 | 0.5 |
| `freq_grid` | construct | adaptive |
| outer tol | 1e-6 | 1e-7 |

### fast plain-grid

> *Maximum practical speed under a documented accuracy envelope.*

Uses a plain fixed/construct frequency grid and the fixed-step Magnus scheme
with no transition refinement, no phase-aware horizon-crossing sub-stepping and
no adaptive frequency grid.  It is a **fast approximate solver with a clear
scientific error budget — not a "gross error" mode** — but that budget is
`NOT VERIFIED` against the continuous-sigma oracle.

- **Where it is acceptable:** exploratory coverage of the signal shape, runtime
  screening, degenerate-probe scans, and the cheap outer tier of a
  parameter-space screen where the final accuracy gate is provided by
  `production` + the oracle.
- **Where it is NOT acceptable:** scientific results or MCMC.  Measured at 9
  matched z8 points against the continuous-sigma reference, the plain-grid
  spectrum relative error median is **1.9e-2** (p95 6.9e-2) and the integrated
  `DN_gw` relative error median is **9.1e-3** (p95 2.7e-2).  The dominant term
  is the fixed-`sigma`-grid bias across the reheating kink, **not** a tuning
  artifact.
- **Runtime (this host, warm, 4 threads):** ≈4.442 ms/point (default A median,
  15 repeats), p95 ≈5.105 ms.  Cold JIT `0.325 s` is reported separately.

### fast transition-refine / production

> *Default scientific-production solver for Cobaya / MCMC.*

The production profile is transition-aware: the reheating kink is an exact
integration breakpoint, horizon crossing uses `phase_max`-capped phase-aware
sub-stepping, the frequency grid is curvature-adaptive, the deep-subhorizon tail
is handed off to an analytic WKB/adiabatic solution at `z_tail`, and every solve
carries a point-local a-posteriori error estimate
(`stiffgwpy.fast_sgwb.estimate_local_error`).

- **Matched z8 accuracy vs the oracle (9 points):** signal/transition spectrum
  relative error max **7.1e-4** (dex max 3.1e-4), gate <1e-3 **PASS**;
  integrated `DN_gw` relative error median **4.3e-4** (p95 1.2e-3) — the <1e-4
  gate is **NOT met**.
- **Axis-edge suite (16 z8 points):** 14 solvable (13 PASS, 1 outlier
  `edge_r_hi` at 1.6e-3 signal-rel) and 2 explicit shared-`Delta_Neff` guard
  rejections (physical, never silent).
- **Parameter space (240 Sobol, production):** 212 ok / 28 explicit guard
  rejections; the artifact runtime median ≈5.3 s/point is a pre-JIT historical
  measurement; current point benchmarks are in `docs/benchmarks.md`.
- **Runtime (this host, warm, 4 threads):** ≈21.772 ms/point (default A median;
  p95 ≈22.149 ms).  Cold JIT `0.226 s` is reported separately.

## Reference / oracle

`stiffgwpy/reference.py` is an independent, higher-order implementation of the
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
pip install .                 # runtime deps: numpy scipy astropy pyyaml numba
pip install .[cobaya]         # + Cobaya MCMC interface (cobaya, mpi4py)
pip install .[dev]            # + pytest, ruff, build, matplotlib
```

## Python usage

```python
from stiffgwpy import LCDM_SG

m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
m.SGWB_iter(engine='fast', accuracy_mode='production', fallback=True)
print(m.DN_gw[-1])

# speed-first pass (plain-grid)
m2 = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
m2.SGWB_iter(engine='fast', accuracy_mode='fast')
```

`accuracy_mode` accepts the two user-facing names (`fast`, `production`) and the
aliases `plain_grid`/`plain-grid`/`transition_refine`/`transition-refine`.
Explicit `h`/`col_step`/`z_tail`/`freq_res`/`tol` override the preset.

```python
from stiffgwpy import fast_sgwb
b = fast_sgwb.estimate_local_error(m)   # point-local 11-category budget
print(b['DN_gw_error'], b['certification_status'])
```

## Cobaya usage

```yaml
theory:
  stiffgwpy.cobaya.stiffGW.stiffGW:
    engine: fast
    fallback: True
    accuracy_mode: production      # fast (plain-grid) | production (transition-refine)
    fast_threads: 8
```

**Mode mapping is unambiguous.**  In `stiffGW.yaml` the fast knobs default to
`0` (a sentinel meaning "use the selected accuracy_mode").  The resolution order
is strictly:

```
accuracy_mode  ->  preset defaults  ->  explicit user overrides only
```

so a default `engine: fast` YAML runs the `production` science settings and you
can never accidentally run `production` as a coarser scan, or have the YAML
default values mask the selected preset.  `accuracy_mode: fast` truly maps to
plain-grid; `accuracy_mode: production` truly maps to transition-refine.

**`eval_freqs` (likelihood bins as native nodes).**  Set
`eval_freqs: [list|path-to-file]` in the theory YAML to force-add
`log10(f/Hz)` native solve nodes.  This is the path by which likelihood
frequency bins reach the fast solver as native nodes
(`SGWB_iter_fast(..., eval_freqs=...)`), removing the interpolation error at
steep spectral features (Layer C measured the per-bin dex interpolation error at
≤3.1e-4 across 11 PTA-like bins).  By default (`eval_freqs: null`) the solver
uses its own grid and the likelihood interpolates over the returned spectrum;
use native nodes when the bin spacing approaches the spectral features.

## Accuracy

Two-layer accuracy guarantees, all vs the continuous-sigma oracle:

| | plain-grid | transition-refine |
|---|---|---|
| spectrum rel (signal, matched z8) | median 1.9e-2 / max 7.0e-2 | max 7.1e-4 |
| spectrum dex (signal, matched z8) | max 8.2e-3 | max 3.1e-4 |
| integrated `DN_gw` rel (median) | 9.1e-3 | 4.3e-4 |
| gate `<1e-3` spectrum | **NOT met** | **PASS** |
| gate `<1e-4` integrated `DN_gw` | **NOT met** | **NOT met** |

Integrated `DN_gw <1e-4` is an honest limit: the residual is at the level of the
oracle's own frozen-tail `z_tail` sensitivity (~3e-4), not a tuning artifact.

## Full parameter validation

Parameter schema (11 physical params, ranges in `scripts/validate_two_modes.py`):

* **Single-parameter axis edges:** `docs/paramsweep_z8b/` (16 points on
  r/n_t/cr/T_re/DN_re/kappa10 axis edges + transition interiors).
* **Param space (LHS, plain-grid screen):** 400 points,
  **255 success / 145 shared-`Delta_Neff` guard / 0 numerical failure**
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
| plain-grid | 4.442 ms warm median; 0.325 s cold | ≈4754x vs current LSODA A run; NOT accuracy-certified |
| transition-refine (production) | 21.772 ms warm median; 0.226 s cold | ≈1017x vs recent LSODA A run; accuracy limits unchanged |
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
reference solve is ≈360 s/point on this host.  The bounded real-Cobaya chains
under `docs/mcmc_posterior/chains/` are ~30-row scaffold runs (adapter plumbing
only, documented as not converged).  Posterior-shift conclusions therefore rest
on importance reweighting, not on an independent reference chain.

## Limitations

* Integrated `DN_gw` relative <1e-4 is NOT met (median 4.3e-4); the residual is
  the reference's own frozen-tail `z_tail` sensitivity, not a knob that can be
  tuned away.
* Fast execution is now over 100x faster than the recent LSODA A-point runtime,
  but this is an execution optimization, not an accuracy certification; the
  plain-grid oracle envelope remains signal median 1.867e-2 / max 7.019e-2.
* The axis-edge suite has a single 1.6e-3 spectrum outlier (`edge_r_hi`,
  r=7.9e-2) — production is not uniformly <1e-3 everywhere in the box; use the
  local error budget / escalation there.
* MCMC validation rests on importance reweighting, not an independent
  reference chain (reference is ~360 s/point).
* The oracle uses a frozen tail; a deep/no-tail certification (z_tail ≥ 14) is
  computationally infeasible because the mode equation becomes deep-subhorizon
  stiff.

## Reproducibility

```bash
python scripts/validate_two_modes.py --phase convergence          # fast-vs-fast convergence
python scripts/validate_two_modes.py --phase param_sweep --n 400  # LHS screen
python scripts/build_two_mode_manifest.py                          # manifest (read-only replay)
python -m pytest                                                   # tests (slow deselected)
python -m pytest -m cobaya                                         # Cobaya adapter gate
python -m build --wheel                                            # wheel build gate
```

Every driver records git-commit + environment metadata.  The regression suite is
99 passed (6 slow LSODA gates deselected by default).

## Directory structure

```
stiffgwpy/            pip package
  stiff_SGWB.py       LCDM_SG + engine dispatch (fast|lsoda|reference)
  fast_sgwb.py        fast solver, ACCURACY_MODES, FAST_PROFILES, estimate_local_error
  reference.py        independent continuous-sigma oracle (engine='reference')
  freq_adaptive.py    curvature-adaptive frequency grid
  exact_background.py continuous-sigma expansion integrals / kink-refined grid
  cobaya/             Cobaya theory adapter + likelihoods
tests/                pytest suite (99 unit + 6 slow gates)
scripts/              validation drivers (validate_two_modes, build_two_mode_manifest, ...)
docs/                 canonical docs + validation artifacts (validation/manifest, paramsweep_*, mcmc_posterior)
Git history           superseded historical audits / benchmarks
```

## License

GPL-3.0 (see `LICENSE.md`).
