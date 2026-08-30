# stiffgwpy

**LCDM + stiff matter + primordial stochastic gravitational-wave background (SGWB)** cosmology code,
with an **experimental approximate fast solver** (`stiffgwpy.fast_sgwb`) alongside the original LSODA implementation.

> **Status (2026-08-29).** The accelerated solver is a new approximate numerical implementation — Numba JIT,
> fixed-step Magnus-type stepping, coarse time-column integration with PCHIP refinement, frequency parallelism
> and precomputed lookup tables. Independent component-level tests confirm it is numerically solid, and repeated
> warm in-process calls are typically **~10^3x faster** than the original LSODA path. However, it has **not yet**
> completed full-parameter-space accuracy certification, full API compatibility or Cobaya integration. Treat it as
> an experimental fast solver, and keep the original LSODA path for cross-validation and as a fallback.

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

### Fast solver (experimental)

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

## Performance (independent re-measurement, 2026-08-29)

- **Warm cache, repeated in-process calls** (the MCMC regime): median-case speedup
  **906–3099x** across the 12-case grid (cross-case median **2162x**; best-value speedup 1465–3483x).
- **First call in a fresh process**: module import 1.3–1.6 s + first JIT solve 5.4–5.7 s, i.e.
  end-to-end **6.8–7.2 s**, only **~3.2–3.4x** faster than the original on the same case.
  Use persistent processes when the fast solver matters.
- Original LSODA: ≈7–23 s per self-consistent run depending on parameters.

Thread scaling is modest (8 threads ≈ 2–3x over 1 thread on an 8-core quota); most of the gain comes
from the new numerical scheme, not from raw thread count.

## Accuracy and validation status

Quick summary of the audit (MCMC speedup, max errors, max ΔlogL, failure
rate, recommended modes): `docs/audit_summary.md`. Independent review of the
external audit text vs. the current repo: `docs/audit_text_review.md`.
Physics-first error budget and an independent high-accuracy reference pipeline:
`docs/audit_error_budget.md`, `docs/audit_reference.md`,
`stiffgwpy/reference.py`.
Requirement-by-requirement acceptance audit (VERIFIED / PARTIALLY VERIFIED /
NOT VERIFIED per the 14 audit sections): `docs/audit_acceptance.md`.

The reference pipeline (`stiffgwpy/reference.py`) is a *different*, higher-order
implementation of the same physics (continuous `sigma(N)` evaluator, so the
instantaneous-reheating kink is not smeared by the fixed-step grid; high-order
adaptive `DOP853` ODE; shape-preserving PCHIP + adaptive Gauss-Kronrod
quadrature with error estimates). It is the new accuracy anchor, replacing
LSODA as the "truth" target. On the default point it shows that fast *and*
LSODA both underestimate the integrated `Delta N_eff` by ~1% because both share
the fixed-step `sigma` grid through the reheating kink (fast −1.32%, LSODA
−0.95% vs the continuous-`sigma` reference). It is intentionally slow
(~30 s / point, 246 frequencies, no parallelism) and is meant for benchmark
points, pathological points and convergence certification, not for MCMC.

Curvature-adaptive frequency sampling (`stiffgwpy/freq_adaptive.py`) refines
the `log10 f` grid where `Omega_GW` is locally curved (spectral knee, stiff
peak, high-f cutoff) and stays sparse in smooth regions. On the default point
it reproduces the fine-grid bolometric integral to ~0.13% (203 points) versus a
~2% error for a 60-point coarse grid, at a similar point count to a 220-point
uniform grid. The sub-horizon-today (very low f) region is physically
ill-defined for a static `Omega_GW` (sign-changing `Ogw - Oj`) and carries no
`Delta N_eff` weight, so the refinement targets the re-entered spectral region.
`grid_independent_freqs` builds the frequency set from continuous background
quantities (not the grid `f_hor`), so the sampling is invariant to the
`sigma`-grid resolution (verified identical for the fixed-grid and variable-grid
models), which is a prerequisite for transition-refined `sigma` without
polluting `Delta N_eff`.

Verified positive results (covered by regression tests in `tests/`):

- `gen_fast` expansion vs original: max abs diff ≈ 1.7e-11 (`sigma`), 1.35e-13 (`f_hor`).
- Custom Simpson weights vs SciPy: worst rel diff ≈ 1e-12; PCHIP refinement vs SciPy ≈ 1.8e-15.
- Fast output is bitwise identical across 1/2/4/9 threads.
- Final `Delta N_eff` on representative cases: rel diff vs original ≈ 5e-5–8e-5.

Known limits (please do not overclaim):

- The `1e-4` tolerance is the **outer self-consistency stopping criterion** of the `Delta N_eff`
  bisection loop, not an ODE error bound nor a uniform accuracy guarantee for every output.
- Spectrum agreement: ≈ 4e-4 dex (absolute); the physical linear-`Omega_GW` relative difference is
  ≈ 8e-4–1e-3 on the 12-case grid and up to ≈ 1.3e-3 on random points.
- Full `DN_gw` evolution curves differ by up to **1%–37%** at the largest relative difference,
  concentrated in the early near-zero region.
- Convergence studies (`h`, `COL_STEP`, `z_tail`, `freq_res`) are complete (see
  `docs/audit_phase2.md`). The deterministic 1030-point parameter-space sweep is
  complete (Sobol 400 + LHS 350 + edge 200 + extreme 80): after de-duplicating
  append-only retries, 782 points are physical (`ok`) and 248 are shared
  fast/LSODA guard rejections. Three transient LSODA `MemoryError` records from
  an overloaded initial run were successfully reproduced as `ok` in the serial
  retry; they remain documented as an infrastructure risk, not solver accuracy.
  On
  `ok` points `DN_gw_last_rel` has max 2.0e-4 and spectrum dex p95 7.05e-2 (max
  0.237); these are same-grid engine diagnostics, not continuum truth claims.
- Cobaya **posterior** comparisons (LSODA chain vs fast chain, `Delta logL`, posterior shift) are
  still **open**; `scripts/mcmc_compare.py` now reports ESS, covariance, KS/Wasserstein/KL and
  refuses to label posterior equivalence certified until both chains exceed the configured
  `--min-effective-samples` threshold (default 2000; output schema `mcmc_compare_v2`).
  Requirement-by-requirement evidence is tracked in `docs/audit_completion_matrix.md`.

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

Four named accuracy modes plus a backward-compatible `ultra-fast` alias (see
`docs/audit_modes.md` for the validation evidence) can be selected through
`SGWB_iter(accuracy_mode=...)` or the Cobaya theory yaml:

| Mode | h | col_step | z_tail | freq_res | outer tol | target use |
|---|---|---|---|---|---|---|
| `debug` | 0.005 | 1 | 10 | 2.0 | 1e-8 | highest grid accuracy + diagnostics |
| `reference` | 0.00125 | 1 | 10 | 2.0 | 1e-8 | tightest fixed-grid mode (slowest) |
| `production` | 0.01 | 4 | 7 | 1.0 | 1e-7 | default for science runs |
| `fast` | 0.02 | 8 | 5 | 1.0 | 1e-6 | fast exploratory scans |
| `ultra-fast` | 0.02 | 8 | 5 | 1.0 | 1e-6 | alias of `fast` |

`fast_sgwb.ACCURACY_MODES` holds the tables; `fast_sgwb.apply_accuracy_mode(name)` applies one and
returns it. `fast_sgwb.get_settings()` snapshots the current module settings.
`fast_sgwb.estimate_error(mode)` returns a calibrated error budget per mode
(`DN_gw_error`, `spectrum_error`, `quadrature_error`, `integration_error`,
`ODE_error`, `tail_error`, `model_bias_error`); the model exposes it as
`m.error_estimates` after a fast solve with an `accuracy_mode`.

`SGWB_iter(engine='fast', sigma_exact=True)` re-computes the expansion integrals
(`F`/`Phi`/`S2`) from the continuous piecewise-exact `sigma` (reheating kink as
an exact breakpoint) instead of the fixed-grid cubic spline. On the default
point this halves the continuous-sigma-vs-grid `model_bias` (0.94% -> 0.42% at
`z_tail=5`). The residual ~0.2-0.4% is the transition-region step-phase
oscillation; the reference pipeline (`stiffgwpy/reference.py`) removes it with
an adaptive high-order ODE, while the fixed-step fast kernel still needs
transition-aware stepping (see `docs/audit_reference.md` §7).

`SGWB_iter(engine='fast', accuracy_mode=..., auto_escalate=True,
error_tol=...)` escalates to `reference` when the calibrated integration error
exceeds `error_tol` (default 5e-3). Note the honest caveat: every fast-grid mode,
*including* `reference`, carries a `model_bias` (the continuous-sigma vs
fixed-grid-sigma bias) that a fast-vs-fast convergence check cannot see. Only
the independent `stiffgwpy/reference.py` pipeline (continuous `sigma(N)`)
removes it, at ~30 s/point, so it is not for MCMC.

## Cobaya

`pip install .[cobaya]` provides the `stiffgwpy.cobaya.stiffGW` theory. The production adapter is
configurable through the theory yaml (all knobs validated at call time):

```yaml
stiffGW:
  engine: fast          # lsoda (default) | fast | reference (slow: continuous-sigma)
  fallback: True        # engine=fast: rerun with LSODA on failure
  fast_threads: 8       # OpenMP threads (0 = module default)
  h: 0.01               # step size (0 = module default)
  col_step: 4           # column stride (0 = module default)
  z_tail: 7.0           # analytic-tail threshold (0 = module default)
  freq_res: 1.0         # frequency-grid density
accuracy_mode: production   # reference | production | ultra-fast | '' (none)
  auto_escalate: False       # escalate to 'reference' when integration error > error_tol
  error_tol: 0.005           # reference-escalation tolerance (relative Delta N_eff)
  reference_rtol: 1e-11      # engine='reference': continuous-sigma ODE tolerance
  reference_z_tail: 5.0      # engine='reference': analytic-tail threshold
```

MPI note: the fast engine spawns no subprocesses (Numba OpenMP threads only), so it is MPI-safe.
The LSODA reference path spawns `SGWB_POOL_SIZE` worker processes per solve; under MPI (world
size > 1) it automatically falls back to 1 worker per rank unless `SGWB_POOL_SIZE` is set.

`engine: reference` also works directly through the native API (not only Cobaya):
`m.SGWB_iter(engine='reference')` runs the independent continuous-sigma
high-accuracy pipeline (`reference.py`) and exposes the same derived outputs; it is
slow (~30-140 s/point) and is for certification/benchmark points, not MCMC.

For observability, `theory.engine_stats` exposes cumulative `fast_evals`,
`fast_failures`, deterministic `fast_guard_rejections`, `lsoda_fallbacks`, and
`fallback_fraction`, plus the last failure reason. The adapter logs the same
summary from `close()` when Cobaya finishes a run and warns above a 5% fallback
fraction. A shared `DN_eff>5` physical rejection is not retried through LSODA;
exceptions, non-finite results, and iteration failures still use `fallback=True`
as the numerical safety mechanism.

## Reproduce the benchmark / validation

```bash
python scripts/bench_fast.py           # warm/cold wall-clock + speedup (12 cases, env+commit metadata)
python scripts/validate_fast.py        # 12-case accuracy gates; exits non-zero on violation
python scripts/validate_random.py      # stratified-random 11-parameter-space gates (P1 starting point)
python scripts/check_random_freq.py    # random 10-frequency spot check + plot
python scripts/benchmark_reference.py --point default --freq-full --with-lsoda  # physics-first benchmark vs the high-accuracy reference
python -m pytest                       # regression tests (slow LSODA gates are deselected by default)
python -m pytest -m slow                # opt in to the two long LSODA reference gates
```

For chain pilots, `scripts/mcmc_compare.py --skip-pointwise` omits the expensive
same-point LSODA diagnostic and is explicitly non-certifying; production
validation must leave this flag off and require the ESS gate.
Certification also requires every finite per-parameter posterior shift to be
within `--max-posterior-shift` (default `0.1`).

For reproducible engine twins, pass `--initial-point <json>`; the file is
applied as deterministic Cobaya `ref` values to both chains while leaving the
prior unchanged. The comparison record stores the point and a SHA-256 hash of
the shared YAML.

All scripts record environment and git-commit metadata; `--json` emits machine-readable records.

## Repository layout

```
stiffgwpy/            pip package (import stiffgwpy)
  stiff_SGWB.py       main model class LCDM_SG (engine='lsoda'|'fast')
  fast_sgwb.py        experimental accelerated solver
  exact_background.py continuous-sigma expansion integrals (sigma_exact path)
  reference.py        independent high-accuracy reference solver (physics-first anchor)
  freq_adaptive.py    curvature-adaptive frequency sampling
  _metrics.py         physically meaningful comparison metrics
  functions.py        FD integrals, LSODA solver
  global_param.py     constants + thermal-history splines (th.dat)
  LCDM_stiff_Neff.py  base cosmology class
  fd_table.npz        precomputed Fermi-Dirac lookup table
  cobaya/             optional Cobaya theory/likelihood interfaces
tests/                pytest regression tests (unit + slow gates)
scripts/              benchmark / validation scripts
docs/                 report + raw validation data
base_param.yml        example parameter file
pyproject.toml        PEP 621 build config
```

## License

GPL-3.0 (see `LICENSE.md`).
