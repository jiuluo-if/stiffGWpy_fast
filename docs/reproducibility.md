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
python -m pytest                                                   # unit/regression tests
python -m pytest -m cobaya                                          # Cobaya adapter gate
python -m build --wheel                                             # wheel build gate
python scripts/smoke_installed_wheel.py dist/stiffgwpy-*.whl         # installed-resource smoke
```

The current regression suite passes 102 tests (6 slow LSODA gates are
deselected by default; opt in with `-m slow`). The installed-wheel smoke test
also checks that core, Cobaya, LIGO, and PTA package resources are available
outside the source checkout.

For performance reproduction, keep compilation separate from execution:

```bash
set FAST_THREADS=4
python scripts/bench_fast.py --reps 15 --cases 0 1 --json docs/benchmark_candidate.json
python scripts/profile_fast_breakdown.py --help
```

The canonical before/after numbers and the 1/2/4/8/16-thread results are in
`docs/performance_comparison_20260903.md`. The maintained Ruff surface and the
narrow configuration type gate are CI-enforced; legacy core lint debt remains
tracked in `docs/engineering_audit.md` and is not represented as PASS.
