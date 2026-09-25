# Test suite

[中文说明](README_zh.md) | English

Tests protect the public API, numerical invariants, solver modes, packaged
resources, and selected scientific comparison gates. Run commands from the
repository root.

## Main areas

- `test_fast_sgwb.py`, `test_engine.py`, and `test_freq_adaptive.py`: fast
  solver outputs, mode resolution, guards, grids, and evaluation nodes.
- `test_reference.py`, `test_numerics.py`, and `test_physics_limits.py`:
  reference helpers and physical limits.
- `test_cobaya_adapter.py`: adapter options, derived parameters, fallback, and
  telemetry.
- `test_compatibility_smoke.py` and `test_resource_budget.py`: supported Python
  imports and runtime resource defaults.
- `test_*_spike.py` and `test_*_round*.py`: isolated candidate or historical
  diagnostic contracts; a passing spike does not promote a solver candidate.
- `test_validation_manifest.py` and `test_posterior_tools.py`: artifact
  structure and analysis helpers.

## Run

```bash
python -m pytest -q                 # Default suite; slow tests are deselected
python -m pytest -m slow -q         # Long-running LSODA/reference checks
python -m pytest -m cobaya -q       # Cobaya integration checks
```

See [`docs/test_coverage_matrix.md`](../docs/test_coverage_matrix.md) for the
risk/cost matrix and [`docs/test_duplication_audit.md`](../docs/test_duplication_audit.md)
for how the CI jobs avoid repeating expensive solves. Test outcomes are bound
to the tested commit and environment; they do not replace artifact-backed
accuracy or parameter-space evidence.
