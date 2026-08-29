# stiffgwpy

**LCDM + stiff matter + primordial stochastic gravitational-wave background (SGWB)** cosmology code,
with a drop-in accelerated solver that is **~1900-3900x faster** than the original LSODA implementation
while keeping the same numerical scheme and the same `1e-4` convergence tolerance.

- Original implementation: per-frequency `scipy.integrate.solve_ivp(method='LSODA')` + Simpson integration
  (≈7-22 s per self-consistent run).
- Accelerated implementation (`stiffgwpy.fast_sgwb`): numba JIT kernels, fixed-step analytic-rotation
  (Magnus-type) stepping with an analytic deep-subhorizon tail, precomputed Simpson weight matrix,
  PCHIP grid refinement and OpenMP parallelism (≈2-6 ms warm).

## Install

```bash
# from source
pip install .

# from GitHub
pip install git+https://github.com/jiuluo-if/stiffGWpy_fast.git

# with Cobaya MCMC interface
pip install .[cobaya]

# development extras
pip install .[dev]
```

Runtime dependencies: `numpy`, `scipy`, `astropy`, `pyyaml`, `numba`.

## Quick start

### New channel (recommended)

```python
from stiffgwpy import LCDM_SG
from stiffgwpy import fast_sgwb

m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
fast_sgwb.SGWB_iter_fast(m)      # fills the same attributes as m.SGWB_iter()

print(m.SGWB_converge)
print(m.log10OmegaGW)            # spectrum
print(m.DN_gw)                   # Delta N_eff evolution
```

### Old channel (unchanged, still works)

```python
from stiff_SGWB import LCDM_SG    # legacy top-level shim
import global_param               # legacy shim, th.dat resolved from package data

m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
m.SGWB_iter()                     # original slow LSODA path (kept)
```

The legacy top-level modules `stiff_SGWB`, `functions`, `global_param`,
`LCDM_stiff_Neff` are thin shims that re-export the packaged modules, so old
scripts keep working unchanged.

## Model and parameters

`LCDM_SG` solves the background (radiation + neutrinos + stiff matter + Lambda),
the tensor-mode perturbation equation per frequency channel, and iterates
`Delta N_eff` (the SGWB contribution to extra radiation) with a bisection loop
until a `1e-4` relative convergence criterion.

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

See `docs/benchmark_report.md` for the full performance/accuracy analysis
(12-combination parameter grid, all ≥1000x, all outputs within the `1e-4`
algorithm tolerance).

## Tuning the fast solver

Set these environment variables before importing `stiffgwpy`:

```bash
export FAST_THREADS=32     # OpenMP threads (default 32)
export FAST_COL_STEP=4     # coarse-column stride, speed/accuracy trade-off (1-8)
```

or at runtime:

```python
from stiffgwpy import fast_sgwb
fast_sgwb.set_threads(32)
fast_sgwb.set_col_step(4)
```

## Reproduce the benchmark / validation

```bash
python scripts/bench_fast.py         # original vs fast wall-clock (12 cases, ~2-4 min)
python scripts/validate_fast.py      # 12-case field-by-field precision validation (~5 min)
python scripts/check_random_freq.py  # random 10-frequency spot check + plot
```

## Repository layout

```
stiffgwpy/            pip package (import stiffgwpy)
  stiff_SGWB.py       main model class LCDM_SG
  fast_sgwb.py        accelerated solver
  functions.py        FD integrals, LSODA solver
  global_param.py     constants + thermal-history splines (th.dat)
  LCDM_stiff_Neff.py  base cosmology class
  cobaya/             optional Cobaya theory/likelihood interfaces
stiff_SGWB.py         legacy top-level shims (old channel)
functions.py
global_param.py
LCDM_stiff_Neff.py
scripts/              benchmark / validation scripts
docs/                 full report + raw validation data
base_param.yml        example parameter file
pyproject.toml        PEP 621 build config
```

## License

GPL-3.0 (see `LICENSE.md`).
