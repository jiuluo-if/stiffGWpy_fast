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
rate, recommended modes): `docs/audit_summary.md`.

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
  `docs/audit_phase2.md`); the 1000-point parameter-space sweep is **91% complete**
  (918/1030 points, stopped on request before the extreme-corner set): `DN_gw_last_rel`
  max 2.0e-4 (inside the `validate_fast` 1e-3 gate), dex tail p95 7.7e-2. All 221
  non-ok points are safety-guard aborts shared by **both** engines (production
  `engine='lsoda'` also returns None there; the sweep reference loop has no guard
  and returns divergent non-physical values `DN_eff > 5`, up to 4.2e13); within
  the physical region (`DN_eff <= 5`, n=697) fast failure rate is 0.00%
  and within the MCMC-relevant region (`DN_eff <= 2`, n=624) spectrum dex error
  is p95 0.019 / max 0.128 with `DN_gw_last_rel` max 2.0e-4 — the large dex tails
  (0.1-0.24) sit only next to the guard boundary (see `docs/audit_phase3.md`).
- Cobaya **posterior** comparisons (LSODA chain vs fast chain, `Delta logL`, posterior shift) are
  still **open**; the adapter (engine/fallback/threads/h/col_step/z_tail/freq_res/accuracy mode)
  is implemented and unit-tested.

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

Three named accuracy modes (audit phase 8 recommendations; see `docs/audit_modes.md` for the
validation evidence) can be selected through `SGWB_iter(accuracy_mode=...)` or the Cobaya theory yaml:

| Mode | h | col_step | z_tail | freq_res | outer tol | target use |
|---|---|---|---|---|---|---|
| `reference` | 0.00125 | 1 | 10 | 2.0 | 1e-8 | closest to the LSODA reference (slowest) |
| `production` | 0.01 | 4 | 7 | 1.0 | 1e-7 | default for science runs |
| `ultra-fast` | 0.01 | 8 | 5 | 1.0 | 1e-6 | fast exploratory scans |

`fast_sgwb.ACCURACY_MODES` holds the tables; `fast_sgwb.apply_accuracy_mode(name)` applies one and
returns it. `fast_sgwb.get_settings()` snapshots the current module settings.

## Cobaya

`pip install .[cobaya]` provides the `stiffgwpy.cobaya.stiffGW` theory. The production adapter is
configurable through the theory yaml (all knobs validated at call time):

```yaml
stiffGW:
  engine: fast          # lsoda (default) | fast
  fallback: True        # engine=fast: rerun with LSODA on failure
  fast_threads: 8       # OpenMP threads (0 = module default)
  h: 0.01               # step size (0 = module default)
  col_step: 4           # column stride (0 = module default)
  z_tail: 7.0           # analytic-tail threshold (0 = module default)
  freq_res: 1.0         # frequency-grid density
  accuracy_mode: production   # reference | production | ultra-fast | '' (none)
```

MPI note: the fast engine spawns no subprocesses (Numba OpenMP threads only), so it is MPI-safe.
The LSODA reference path spawns `SGWB_POOL_SIZE` worker processes per solve; under MPI (world
size > 1) it automatically falls back to 1 worker per rank unless `SGWB_POOL_SIZE` is set.

## Reproduce the benchmark / validation

```bash
python scripts/bench_fast.py           # warm/cold wall-clock + speedup (12 cases, env+commit metadata)
python scripts/validate_fast.py        # 12-case accuracy gates; exits non-zero on violation
python scripts/validate_random.py      # stratified-random 11-parameter-space gates (P1 starting point)
python scripts/check_random_freq.py    # random 10-frequency spot check + plot
python -m pytest                       # regression tests (fast unit tests + 2 slow LSODA gates)
```

All scripts record environment and git-commit metadata; `--json` emits machine-readable records.

## Repository layout

```
stiffgwpy/            pip package (import stiffgwpy)
  stiff_SGWB.py       main model class LCDM_SG (engine='lsoda'|'fast')
  fast_sgwb.py        experimental accelerated solver
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
