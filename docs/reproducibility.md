# Reproducibility

Status: current verification procedure; evidence dates are listed per artifact
Date: 2026-09-24
Code version: record the tested HEAD and environment with each run

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
python scripts/validate_manifest.py                                # committed artifact schema
python -m build --wheel                                             # wheel build gate
python scripts/smoke_installed_wheel.py dist/stiffgwpy_fast-*.whl         # installed-resource smoke
```

Test totals vary by commit and marker selection. One recorded full GitHub CI
matrix at commit `e3123a2` passed all nine jobs. Later commits require their
own CI result; do not infer exact-HEAD status from this record. The
installed-wheel smoke test checks that core, Cobaya, LIGO, and PTA package
resources are available outside the source checkout.

For performance runs, separate compilation from timed execution. Set
`FAST_THREADS` in the shell before running the benchmark:

```powershell
$env:FAST_THREADS = "4"
python scripts/bench_fast.py --reps 15 --cases 0 1 --json docs/benchmark_candidate.json
python scripts/profile_fast_breakdown.py --help
```

The 2026-09-23 fixed-resource 25-repeat full-path profile is summarized in
`docs/benchmarks.md`; the detailed 2026-09-03 before/after numbers and its
1/2/4/8/16-thread results are historical evidence in
`docs/performance_comparison_20260903.md`. The maintained Ruff surface and the
narrow configuration type gate are CI-enforced; legacy core lint debt remains
tracked in `docs/engineering_audit.md` and is not represented as PASS.
