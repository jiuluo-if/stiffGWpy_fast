# stiffgwpy

**LCDM + stiff matter + primordial stochastic gravitational-wave background (SGWB)** cosmology code,
**LCDM + stiff matter + primordial stochastic gravitational-wave background (SGWB)** cosmology code,
with an accelerated fast solver (`stiffgwpy.fast_sgwb`), the original LSODA path, and an independent
continuous-`sigma(N)` high-accuracy reference pipeline (`stiffgwpy/reference.py`). The reference
pipeline is the accuracy oracle for every validation claim below; **LSODA is never the truth anchor**.

> **Status (2026-09-03).** Three-layer fast-vs-reference certification complete (host-measured;
> per-gate table: `docs/audit_acceptance.md`).
> - Single-point physics (9 matched z8 cases): signal/transition-region `Omega_GW` relative error max
>   **7.09e-4** (gate <1e-3 PASS); integrated `Delta_Neff` residual median **4.3e-4** — the <1e-4
>   gate is honestly NOT met at the current Magnus/grid architecture level.
> - Parameter space (240 Sobol, production z8): 212 ok / 28 explicit shared-background guard
>   rejections; posterior-bulk likelihood-bin dex max **3.10e-4**; |Delta logL| max **7.30e-3**.
> - Posterior (importance sampling, fixed seed): ESS **4167** (gate >=2000 PASS); `log10 r` posterior
>   shift **-0.0011 sigma** (gate <0.1 sigma PASS).
> - Cobaya adapter (production-grade): engines `fast` / `lsoda` / `reference`, `fallback`,
>   likelihood-aware `auto_escalate` (upgrade to the continuous-sigma reference), per-evaluation
>   status `FAST` / `FAST_ESCALATED` / `REFERENCE` / `LSODA_FALLBACK`; no silent fallback.
> - Runtime (honest): production z8 ≈ **4.1 s/point vs the 18.56 s LSODA anchor ≈ 4.5x**, with an
>   integrated-`Delta_Neff` residual ≈6x closer to the reference at the default point. The ~1000x
>   warm numbers below hold only for the coarser plain-grid mode, which does not meet the 1e-3 gate.

## Install

```bash
# from source
pip install .

# from GitHub
pip install git+https://github.com/jiuluo-if/stiffGWpy_fast.git

# with Cobaya MCMC interface
pip install .[cobaya]

# development extras (pytest, ruff, ...)
pip install .[dev]
```

Runtime dependencies: `numpy`, `scipy`, `astropy`, `pyyaml`, `numba`.

## Quick start

### Fast solver

```python
from stiffgwpy import LCDM_SG

m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
m.SGWB_iter(engine='fast', fallback=True)   # fast path; falls back to LSODA on failure
```

`engine='fast'` is the unified entry point added in this revision (default remains `engine='lsoda'`).
The low-level `fast_sgwb.SGWB_iter_fast(m)` also still works. On success the fast path fills the core
outputs (`log10OmegaGW`, `DN_gw`, `kappa_r`, `hubble`, `g2`, `w2`); it does **not** produce the auxiliary
arrays `N_hc`, `Th`, `Oj`, `Ogw`, `Opgw` of the original solver, so it is **not** a drop-in replacement
for code that reads those attributes.

The unified entry point accepts the audit-tuned fast-path knobs directly:

```python
m.SGWB_iter(engine='fast', accuracy_mode='production')   # named preset, see below
m.SGWB_iter(engine='fast', h=0.005, col_step=4, z_tail=7.0,
            freq_res=1.0, threads=8, tol=1e-7)           # explicit overrides
```

Explicit `h`/`col_step`/`threads` and non-default `z_tail`/`freq_res`/`tol` override the preset.
These settings are process-global module state (same semantics as the setters below).

### Original LSODA path (unchanged)

```python
m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
m.SGWB_iter()            # engine='lsoda', original slow path
```

Import everything through the package (`stiffgwpy.*`); the legacy top-level modules
`stiff_SGWB` / `functions` / `global_param` / `LCDM_stiff_Neff` were removed to avoid a
duplicate import channel (old scripts should switch to `from stiffgwpy import LCDM_SG`).

## Performance and runtime (measured 2026-09-03)

Matched-accuracy production numbers — the honest comparison (default point unless noted):

- **fast production z8 ≈ 4.1 s per self-consistent point**; a 240-point Sobol sweep at the same
  settings gives median 5.34 s, p95 8.82 s, max 14.2 s (this host, 8 threads).
- The documented LSODA anchor is **18.56 s/point** (z5 settings) → **≈4.5x** faster, while the
  integrated-`Delta_Neff` error vs the continuous-sigma reference drops from LSODA's **-1.68e-3** to
  fast production's -7.6e-4..+3.0e-4 at the default point (≈6x closer to the reference).
- A 100x+ speedup vs LSODA holds only for the plain-grid coarser mode (0.012-0.24 s/pt,
  `Delta_Neff` ≈ -1.6e-3 relative vs the reference), which does not pass the 1e-3 physics gate.

Warm-cache coarsest-mode measurements (independent re-measurement, 2026-08-29; repeated
in-process calls — the regime behind the historical `10^3x` headline):

- Median-case speedup **906-3099x** across the 12-case grid (cross-case median **2162x**;
  best-value speedup 1465-3483x).
- First call in a fresh process: module import 1.3-1.6 s + first JIT solve 5.4-5.7 s, i.e.
  end-to-end 6.8-7.2 s, only ~3.2-3.4x faster than the original on the same case. Use persistent
  processes when the fast solver matters.
- Original LSODA: ≈7-23 s per self-consistent run depending on parameters.

Thread scaling is modest (8 threads ≈ 2–3x over 1 thread on an 8-core quota); most of the gain comes
from the new numerical scheme, not from raw thread count.

## Accuracy and validation status (2026-09-03)

Physics-first error budget: `docs/audit_error_budget.md`. Independent high-accuracy reference
solver: `stiffgwpy/reference.py` + `docs/audit_reference.md`. Requirement-by-requirement acceptance
audit with the current per-gate status table: `docs/audit_acceptance.md`. Cross-cutting
speed-vs-accuracy Pareto: `docs/audit_speed_accuracy.md`.

The reference pipeline is a *different*, higher-order implementation of the same physics:
continuous `sigma(N)` (the instantaneous-reheating kink is an exact breakpoint, never smeared by a
grid), adaptive high-order `DOP853` per frequency mode, shape-preserving PCHIP + adaptive
Gauss-Kronrod quadrature with error estimates. It is the accuracy anchor for every conclusion in
this README. It is intentionally slow — ≈360 s/point at the matched z8 242-mode setting (rtol=1e-9)
and ≈935 s for the 484-mode deep oracle (rtol=1e-10) on this host — and is for certification and
benchmark points, not for MCMC.

Validation is organised in three layers; fast and reference always solve the SAME grid and the SAME
`z_tail` so the residual isolates engine error.

**Layer A — single-point physics (9 matched z8 cases).** `docs/paramsweep_z8/validation_summary.md`
(+ `validation_summary.json`): default / stiff / low-T / high-T / radiation-dominated / tiny-r /
transition-sensitive / cr0-blue / extreme points.
- signal band `log10(f/Hz) in [-6, +1]`: relative error max **7.09e-4**, dex max 3.08e-4
  (gate <1e-3: PASS); transition band: relative max 7.09e-4 (PASS).
- integrated `Delta_Neff` rel abs: median 4.34e-4, p95 1.18e-3, max 1.46e-3 (low-T point; a
  DN-of-DN artifact: DN ~ 5.2e-8, both engines agree to 5.5e-10 absolute) — gate <1e-4: **NOT met**.

**Layer B — parameter-space validation (240 Sobol, production z8).** `docs/paramsweep_ref/fast_sweep.jsonl`:
212/240 `ok`; the 28 rejections are the explicit shared-background `Delta_Neff` guard in extreme
(r, DN_re, kappa10) corners (documented, never silent). Runtime median 5.34 s/point; adaptive
frequency grid median 236 nodes (p95 266, max 286); per-mode WKB-handoff defect median 6.84e-4;
point-local 11-category error estimate (`fast_sgwb.estimate_local_error`) median combined relative
`Delta_Neff` error 4.0e-4 (saturating at ~1.0 only where DN -> 0; absolute error <= 1e-5,
physically unobservable).

**Layer C — likelihood-aware posterior validation.** `docs/mcmc_posterior/posterior_validation.md`:
240 posterior-bulk points solved by BOTH engines at the same 11 likelihood bins, with the bins as
native solve nodes (`SGWB_iter_fast(..., eval_freqs=...)`): per-bin dex max **3.10e-4**,
|Delta logL| max **7.30e-3** (gate 0.1: PASS). IS posterior from 9000 fast-production draws:
ESS **4167** (gate 2000: PASS). e^{Delta logL}-reweighting of the same 240 points gives the
reference-consistent posterior: `log10 r` shift **-0.0011 sigma**, n_t +0.0002 sigma (gate <0.1
sigma: PASS; n_t is prior-dominated under cr=1 — documented, not certified).
Physics-first fixes behind these numbers (why each is more correct in its regime):

- The dominant historical error — the shared fixed-step `sigma`-grid bias through the reheating
  kink (fast ≈ -1.3% vs the reference at the default point) — is fixed by treating the transition
  as an ODE integration breakpoint: production `transition_refine` uses a kink-aware grid with the
  kink inside a refined sub-step (never crossed by a spline/grid), plus `phase_max=0.5` caps the
  per-substep phase advance `e^z dh` around horizon crossing (z = ln(k/aH) event positions enter
  the step control analytically). Richardson check: |dDeltaN(phase_max 0.5 -> 0.125)| = 7.5e-6
  relative. Matched-z8 residual: median 4.3e-4 (was ~1%).
- Deep-subhorizon modes hand off to an analytic frozen-tail (WKB/adiabatic) solution at `z_tail`;
  the per-mode adiabaticity defect `handoff_eps` is returned per solve as a verifiable local error
  estimate (median 6.8e-4 over the sweep).
- Curvature-adaptive frequency sampling (`freq_adaptive.py`) is the production grid: refine
  `log10 f` where `Omega_GW` is locally curved (spectral knee, stiff enhancement, high-f cutoff),
  stay sparse in smooth regions (median 236 nodes/point, Layer B). The sub-horizon-today
  (very low f) region is physically ill-defined for a static `Omega_GW` (sign-changing `Ogw - Oj`)
  and carries no `Delta_Neff` weight, so refinement targets the re-entered spectral region.
- Every solve carries an 11-category a-posteriori local error budget (background/model,
  sigma-transition, ODE integration, horizon-crossing, WKB handoff, interpolation, frequency grid,
  quadrature, tail approximation, floating-point/cancellation, self-consistency):
  `estimate_local_error` combines systematic (max of model/transition) + RSS of the rest — not a
  fixed benchmark constant.
- Analytic limits, energy/scaling consistency, float robustness and WKB-frequency checks are green
  (`tests/test_physics_limits.py`). Full suite: **105 passed** (including 6 slow-marked gates).

Honest limits (reported because they are real, not threshold-tuned):

- Integrated-`Delta_Neff` relative < 1e-4 is NOT met: median 4.3e-4 at matched z8 (freq_res=1) and
  -2.94e-4 at the deep-oracle default point (freq_res=2, rtol=1e-10). Per-mode residuals
  (2.6-3.1e-4 dex) are internally converged vs h / z_tail / freq_res / phase_max (<1e-4 dex); the
  remaining term is the frozen-z Magnus + grid architecture envelope, not a tuning artifact.
- Production 100x vs LSODA is NOT met at the matched-accuracy setting (≈4.5x; see Performance).
- Layer C rests on the IS posterior (ESS 4167) + exact e^{Delta logL} reweighting; full two-chain
  KS/Wasserstein/KL/covariance comparisons were not run because a reference-engine MCMC chain costs
  ≈350-935 s/point on this host. The bounded real-Cobaya chains under `docs/mcmc_posterior/chains/`
  are ~30-row scaffold runs (adapter plumbing only; documented as not converged).
- Historical fixed-grid-era claims in earlier audit documents (~1% model bias, same-grid engine dex
  p95 7e-2, posterior comparisons still open) describe the pre-fix architecture and are superseded
  by this section and by `docs/audit_acceptance.md` §0.

## Model and parameters

`LCDM_SG` solves the background (radiation + neutrinos + stiff matter + Lambda), the tensor-mode
perturbation equation per frequency channel, and iterates `Delta N_eff` (the SGWB contribution to
extra radiation) with a bisection loop until the outer `1e-4` relative convergence criterion.

Main base parameters (can be given as kwargs, dict or YAML file):

| Parameter | Meaning |
|---|---|
| `Omega_bh2`, `Omega_ch2`, `H0` | baryon/CDM densities and Hubble constant |
| `DN_eff` | constant extra radiation (Delta N_eff) |
| `A_s`, `r`, `n_t` | primordial scalar amplitude, tensor-to-scalar ratio, tensor tilt |
| `cr` | >0: enforce the single-field consistency relation; <=0: use given `n_t`/`DN_re` |
| `T_re` | reheating temperature [GeV] |
| `DN_re` | e-folds of matter-like reheating |
| `kappa10` | `rho_stiff / rho_photon` at 10 MeV |

## Tuning the fast solver and accuracy modes

Set these environment variables **before importing** `stiffgwpy`:

```bash
export FAST_THREADS=8     # OpenMP threads; default = numba's own default (>=1, <= detected cores)
export FAST_COL_STEP=4    # coarse-column stride, speed/accuracy trade-off (1-8, default 4)
export FAST_H=0.01        # fixed step / expansion-grid spacing (1e-4 .. 0.1)
export FAST_Z_TAIL=5.0    # analytic deep-subhorizon tail threshold (2.0 .. 15.0)
export SGWB_POOL_SIZE=4   # LSODA reference-path frequency-parallel workers (>=1);
                          # under MPI (world size > 1) the default is 1 unless set
```

or at runtime:

```python
from stiffgwpy import fast_sgwb
fast_sgwb.set_threads(8)     # validated against numba's detected core count
fast_sgwb.set_col_step(4)    # validated to 1..8
fast_sgwb.set_h(0.01)        # validated to 1e-4..0.1
fast_sgwb.set_z_tail(5.0)    # validated to 2.0..15.0
```

All setters validate their input; invalid values raise `ValueError` instead of failing later.

Four named accuracy modes plus a backward-compatible `ultra-fast` alias and a `deep` validation
variant are selected through `SGWB_iter(accuracy_mode=...)` or the Cobaya theory yaml. The
authoritative table lives in `fast_sgwb.ACCURACY_MODES`:

| Mode | h | col_step | z_tail | freq_res | freq grid | phase_max | outer tol | use |
|---|---|---|---|---|---|---|---|---|
| `debug` | 0.005 | 1 | 10 | 2.0 | uniform | 0.25 | 1e-8 | highest grid accuracy + diagnostics |
| `reference` | 0.00125 | 1 | 10 | 2.0 | uniform | 0.1 | 1e-8 | tightest fixed-grid mode (slowest) |
| `production` | 0.01 | 4 | 8 | 1.0 | adaptive | 0.5 | 1e-7 | default science / certified fast config |
| `deep` | 0.01 | 4 | 10 | 1.0 | adaptive | 0.25 | 1e-7 | deep-tail validation variant |
| `fast` | 0.02 | 8 | 5 | 1.0 | uniform | 0.0 | 1e-6 | fast exploratory scans |
| `ultra-fast` | 0.02 | 8 | 5 | 1.0 | uniform | 0.0 | 1e-6 | alias of `fast` |

`production` (and `reference`/`debug`) enable `transition_refine` (kink-aware grid: the reheating
transition is an integration breakpoint, σ(N) never crosses a discontinuity on a spline/grid) plus
`phase_max`-capped Magnus sub-stepping (per-substep phase `e^z dh` bounded). At matched z8 this
preset delivers the Layer A/B/C numbers above at ≈4-5 s/point; `fast`/`ultra-fast` keep the plain
grid for quick exploration (0.012-0.24 s/pt, `Delta_Neff` ≈ -1.6e-3 relative vs the reference —
below the 1e-3 physics gate).

`fast_sgwb.estimate_local_error(m)` returns a point-local a-posteriori error budget for the last
fast solve: 11 physics categories (background_model / sigma_transition / ode_integration /
horizon_crossing / wkb_handoff / interpolation / frequency_grid / quadrature /
tail_approximation / floating_point / self_consistency), each tagged `local` (computed from that
solve's telemetry: `handoff_eps`, `freq_grid_error`, `phase_max_used`, `z_tail_used`, quadrature
Richardson, cancellation ratio, bisection bracket) or `calibrated` (measured default-point anchors
scaled to the solve's settings), plus combined `DN_gw_error` and `Delta_Neff_abs_error`. The model
exposes the same after a fast solve (`m.error_estimates`, `m.DN_gw_error`,
`m.local_error_budget`).

`SGWB_iter(engine='fast', sigma_exact=True)` recomputes the expansion integrals F/Phi/S2 from the
continuous piecewise-exact `sigma` (kink as an exact breakpoint) instead of the fixed-grid spline.
The certified production path instead uses kink-aware `transition_refine` stepping (above); both
remove the grid-smeared-kink bias that used to dominate (~1%).

`SGWB_iter(engine='fast', accuracy_mode=..., auto_escalate=True, error_tol=...)` re-solves with
the tightened fast `reference` grid when the calibrated local error exceeds `error_tol` (default
5e-3); with `escalate_to_reference=True` it runs the independent continuous-sigma reference
pipeline instead. The engine-level gate is likelihood-aware when `likelihood_sigma` is given:
|Delta logL| = 0.5*(Delta_Neff_abs_error/likelihood_sigma)^2 is compared against `dlogl_tol`
(default 1e-3). Honest caveat: a fast-vs-fast escalation cannot see continuous-sigma model bias;
only the independent `reference.py` pipeline removes it (~360-935 s/point), so it is not for the
MCMC thermal path.

## Cobaya

`pip install .[cobaya]` provides the `stiffgwpy.cobaya.stiffGW` theory. The production adapter is
configurable through the theory yaml (all knobs validated at call time):

```yaml
stiffGW:
  engine: fast              # lsoda (default, conservative) | fast | reference (continuous-sigma)
  fallback: True            # engine=fast: rerun with LSODA on numerical failure (tagged LSODA_FALLBACK)
  fast_threads: 8           # OpenMP threads (clamped to the machine budget at call time)
  h: 0.01                   # step size (0 = module/preset default)
  col_step: 4               # column stride (0 = module/preset default)
  z_tail: 7.0               # analytic-tail threshold (0 = preset default; Layer A-C certification used 8.0)
  freq_res: 1.0             # frequency-grid density (1.0 = default)
  accuracy_mode: production # reference | production | deep | ultra-fast | '' (none)
  auto_escalate: False      # escalate when the estimated error exceeds the budget
  error_tol: 0.005          # relative Delta_Neff error gate (engine-level, no likelihood sigma)
  likelihood_sigma: null    # effective sigma of the likelihood on Delta_Neff (likelihood-aware gate)
  dlogl_tol: 0.001          # |Delta logL| budget in natural log-units
  escalate_to_reference: False  # False: retighten fast grid; True: run continuous-sigma reference
  reference_rtol: 1e-11     # engine=reference: continuous-sigma ODE tolerance
  reference_z_tail: 5.0     # engine=reference: analytic-tail threshold
```

MPI note: the fast engine spawns no subprocesses (Numba OpenMP threads only), so it is MPI-safe.
The LSODA reference path spawns `SGWB_POOL_SIZE` worker processes per solve; under MPI (world
size > 1) it automatically falls back to 1 worker per rank unless `SGWB_POOL_SIZE` is set.

`engine: reference` also works directly through the native API: `m.SGWB_iter(engine='reference')`
runs the independent continuous-sigma high-accuracy pipeline (`reference.py`) and exposes the same
derived outputs; it is slow (~360-935 s/point at z8 on this host) and is for certification and
benchmark points, not MCMC.

For observability, `theory.engine_stats` exposes cumulative `fast_evals`, `fast_failures`,
`fast_guard_rejections`, `fast_physical_rejections`, `lsoda_evals`, `lsoda_fallbacks`,
`reference_evals`, `escalations`, `fallback_fraction`, `escalation_fraction`, plus per-evaluation
status (`last_eval_status`, `eval_status_counts` over FAST / FAST_ESCALATED / REFERENCE / LSODA /
LSODA_FALLBACK), `dlogl_estimated`, `DN_gw_error` and the last fast failure reason. The adapter
logs the same summary from `close()` when Cobaya finishes a run and warns above a 5% fallback or
escalation fraction. A shared `Delta_Neff > 5` physical rejection is explicit
(`fast_failure_reason='shared_Neff_guard'`) and is not retried through LSODA; numerical failures
(exceptions, non-finite results, iteration failures) use `fallback=True` and are tagged
`LSODA_FALLBACK` — there is no silent fallback.

## Reproduce the benchmark / validation

```bash
python scripts/bench_fast.py                      # warm/cold wall-clock (12 cases, env+commit metadata)
python scripts/validate_fast.py                   # 12-case accuracy gates; exits non-zero on violation
python scripts/validate_random.py                 # stratified-random parameter-space gates
python scripts/check_random_freq.py               # random 10-frequency spot check + plot
python scripts/benchmark_reference.py --point default --freq-full  # physics-first benchmark vs reference
python scripts/validate_fast_vs_reference.py      # Layer A: 9 matched single points, fast vs continuous-sigma reference (z8)
python scripts/importance_posterior.py --help     # Layer C: IS posterior validation phases (draw/posterior/pointwise/report)
python scripts/cobaya_posterior_fast_vs_reference.py  # bounded real-Cobaya scaffold (adapter plumbing)
python -m pytest                                  # regression tests (slow gates deselected by default)
python -m pytest -m slow                          # opt in to the 6 long slow gates
```

Artifacts: Layer A/B summaries in `docs/paramsweep_z8/`; the 240-point production sweep in
`docs/paramsweep_ref/fast_sweep.jsonl`; Layer C posterior artifacts in `docs/mcmc_posterior/`
(`is_report.json`, `is_pointwise.json`, `is_draws.npz`, `posterior_validation.md`); deep-oracle
default anchor in `docs/reference/deep_oracle_default.json`.

The certified Layer C driver (`scripts/importance_posterior.py`) uses 9000 Gaussian-proposal
importance draws with a fixed seed (20260903), enforces ESS >= 2000 on the fast posterior, and
converts it to the reference-consistent posterior by exact e^{Delta logL} reweighting of 240
posterior-bulk points solved by both engines (gates: |Delta logL| < 0.1, posterior shift
< 0.1 sigma). `scripts/cobaya_posterior_fast_vs_reference.py` exercises the real Cobaya adapter
with identical priors/initial points/seeds but is a short scaffold run (documented as not
converged; validates adapter plumbing only).

All scripts record environment and git-commit metadata; `--json` emits machine-readable records.

## Repository layout

```
stiffgwpy/            pip package (import stiffgwpy)
  stiff_SGWB.py       main model class LCDM_SG (engine='lsoda'|'fast'|'reference')
  fast_sgwb.py        accelerated solver + ACCURACY_MODES + estimate_local_error
  exact_background.py continuous-sigma expansion integrals (sigma_exact path)
  reference.py        independent high-accuracy reference solver (physics-first anchor)
  freq_adaptive.py    curvature-adaptive frequency sampling (production grid)
  _metrics.py         physically meaningful comparison metrics
  functions.py        FD integrals, LSODA solver
  global_param.py     constants + thermal-history splines (th.dat)
  LCDM_stiff_Neff.py  base cosmology class
  fd_table.npz        precomputed Fermi-Dirac lookup table
  cobaya/             Cobaya theory adapter (engine fast|lsoda|reference + telemetry)
tests/                pytest regression tests (unit + 6 slow gates)
scripts/              benchmark / validation drivers (incl. Layer A-C certification)
docs/                 audit reports + validation artifacts (paramsweep_z8/, mcmc_posterior/)
base_param.yml        example parameter file
pyproject.toml        PEP 621 build config
```

## License

GPL-3.0 (see `LICENSE.md`).
