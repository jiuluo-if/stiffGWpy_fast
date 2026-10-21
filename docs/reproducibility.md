# Reproducibility

Status: current verification procedure; evidence dates are listed per artifact
Date: 2026-10-21
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

Test totals vary by commit and marker selection. The 2026-10-21 CI run
`39619165214` for code commit `32dbabb` passed all nine jobs. Later commits
require their own CI result; do not infer exact-HEAD status from this record.
The installed-wheel smoke test checks that core, Cobaya, LIGO, and PTA package
resources are available outside the source checkout.

For performance runs, separate compilation from timed execution and record
the thread layer, thread budget, BLAS budget, CPU affinity, and cold/warm split.
The current full-path and concurrency records link their exact commands in
`docs/benchmarks.md`. To reproduce the 2026-10-21 two-call TBB throughput A/B:

```powershell
$env:NUMBA_NUM_THREADS = "2"
$env:NUMBA_THREADING_LAYER = "threadsafe"
$env:OMP_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
python scripts/benchmark_backend_concurrency.py --output docs/benchmark_backend_concurrency_replay.json --repeats 80 --threads 2 --cases default,highT --affinity 0,1,2,3
```

The script records the selected layer and environment. It rejects a backend
that Numba does not classify as thread-safe; the timing is concurrent-call
throughput, not single-evaluation latency.

The 2026-09-23 fixed-resource 25-repeat full-path profile is summarized in
`docs/benchmarks.md`; the detailed 2026-09-03 before/after numbers and its
1/2/4/8/16-thread results are historical evidence in
`docs/performance_comparison_20260903.md`. The maintained Ruff surface and the
narrow configuration type gate are CI-enforced; legacy core lint debt remains
tracked in `docs/engineering_audit.md` and is not represented as PASS.
