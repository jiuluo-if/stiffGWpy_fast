# stiffgwpy_fast

`stiffgwpy_fast` is a Python package for calculating stochastic gravitational-wave backgrounds in a cosmology with stiff matter. The `fast` engine is the default; the package also provides independent `reference` and `lsoda` engines for comparisons and fallback.

## Install

```bash
python -m pip install .
```

For the Cobaya adapter, install `.[cobaya]`. For development and tests, install `.[dev]`.

## Quick start

```python
from stiffgwpy_fast import LCDM_SG

model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
model.SGWB_iter()  # default fast engine
```

The model exposes the computed spectrum and derived cosmological quantities. See the [Cobaya adapter guide](docs/cobaya.md) when using the optional integration.

## Numerical method and validation

The production `fast` profile uses the package's Numba propagation kernels, reheating-kink-aware background preparation, and adaptive frequency integration. The `reference` engine uses an independent continuous-background DOP853 solve; `lsoda` is retained for regression comparisons and fallback. See [the physical model](docs/physics.md), [numerical method](docs/numerical_method.md), and [accuracy boundaries](docs/accuracy.md).

Core regression tests include fixed reference samples in `tests/data/`. Those samples are test fixtures; bulk scans, experiment scripts, reports, and generated plots are maintained locally outside the release repository.

## Development

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
ruff check stiffgwpy_fast scripts tests
python -m mypy
python -m build --sdist --wheel
```

Packaging details are in [PACKAGING.md](PACKAGING.md), contribution requirements in [CONTRIBUTING.md](CONTRIBUTING.md), and reproducibility guidance in [docs/reproducibility.md](docs/reproducibility.md).
