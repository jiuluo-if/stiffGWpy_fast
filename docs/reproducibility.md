# Reproducibility

Status: current
Date: 2026-09-03
Code version: see manifest `commit`

Validations are driven by committed scripts; every driver records git-commit +
environment metadata (`python`/`numpy`/`numba`/`scipy`/platform).

```bash
python scripts/validate_two_modes.py --phase convergence          # fast-vs-fast convergence
python scripts/validate_two_modes.py --phase param_sweep --n 400  # LHS screen
python scripts/build_two_mode_manifest.py                          # manifest (read-only replay)
python scripts/validate_fast_vs_reference.py --phase reference     # z8 matched oracle suite
python scripts/validate_plain_grid_vs_reference.py --phase plain   # plain-grid oracle suite
python scripts/validate_edges_vs_reference.py --phase reference     # axis-edge suite
python scripts/importance_posterior.py --phase all                 # Layer C posterior
python -m pytest                                                   # 99 unit tests
python -m pytest -m cobaya                                          # Cobaya adapter gate
python -m build --wheel                                             # wheel build gate
```

The regression suite is 99 passed (6 slow LSODA gates deselected by default;
opt in with `-m slow`).
