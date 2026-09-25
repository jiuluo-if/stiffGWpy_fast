# `stiffgwpy_fast` package

[中文说明](README_zh.md) | English

This directory contains the installed Python package. The repository overview
and user quick start are in the [root README](../README.md).

## Engine roles

- `fast` is the sole formal user-facing profile. It uses the goal frequency
  grid, exact reheating-kink split, Numba tensor propagation, and PCHIP
  bolometric integration.
- `reference` is the independent continuous-`sigma(N)` DOP853 precision
  pipeline. It is intended for scoped comparisons and certification runs.
- `lsoda` retains the original adaptive integration path for regression and
  explicit fallback.

Historical names such as `production` and `ultra-fast` remain for validation
compatibility. The high-level API maps deprecated fast-mode names to `fast`;
they are not additional production profiles.

## Module map

| Module | Responsibility |
|---|---|
| `stiff_SGWB.py` | Public `LCDM_SG.SGWB_iter` entry point and engine dispatch |
| `LCDM_stiff_Neff.py` | Cosmological model parameters and derived quantities |
| `fast_sgwb.py` | Fast presets, Numba kernels, iteration, frequency integration, local error budget |
| `exact_background.py` | Continuous-background values and reheating-kink-aware primitives for fast preparation |
| `freq_adaptive.py` | Adaptive and goal-oriented frequency grids |
| `reference.py` | Independent continuous-background DOP853 reference |
| `config.py` | Immutable, validated per-call fast configuration |
| `_resources.py` | Context-managed access to packaged data files |
| `_metrics.py` | Shared spectrum and relative-error metrics |
| `global_param.py`, `functions.py` | Physical constants and legacy numerical helpers |
| `cobaya/` | Cobaya theory adapter and likelihood modules |

## Configuration boundary

`LCDM_SG.SGWB_iter()` selects the formal `fast` profile by default. Named
configuration is resolved per call; `FastSolverConfig` is immutable. Its bare
constructor defaults preserve legacy settings and are not the formal preset;
use `fast_sgwb.resolve_config('fast')` to obtain the named profile. Pass
`accuracy_mode=None` only when deliberately using legacy module settings.
Explicit parameter overrides apply only to the current call.

```python
from stiffgwpy_fast import LCDM_SG

model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
model.SGWB_iter()  # Uses the formal fast profile
```

See [numerical methods](../docs/numerical_method.md),
[accuracy boundaries](../docs/accuracy.md), and
[reproducibility guidance](../docs/reproducibility.md) before interpreting
dated numerical results.
