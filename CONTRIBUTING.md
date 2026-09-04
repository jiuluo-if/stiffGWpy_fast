# Contributing

Install the development dependencies with `pip install .[dev]`. The normal
local gate is:

```bash
python scripts/validate_manifest.py
ruff check stiffgwpy/__init__.py stiffgwpy/cobaya/__init__.py stiffgwpy/cobaya/stiffGW.py stiffgwpy/_metrics.py stiffgwpy/_resources.py stiffgwpy/config.py stiffgwpy/global_param.py stiffgwpy/exact_background.py stiffgwpy/freq_adaptive.py stiffgwpy/reference.py scripts tests
python -m mypy
python -m pytest -q
python -m pytest -q -m slow
python -m build --sdist --wheel
python scripts/verify_distribution.py dist
python scripts/smoke_installed_wheel.py dist/stiffgwpy-*.whl
```

The independent continuous-sigma reference is the precision anchor. LSODA is
for regression and runtime comparison only. Do not turn a `NOT VERIFIED` or
`FAIL` artifact into a passing claim by changing documentation alone; update
the executable validation and its provenance together.

Generated build/cache files are ignored. The pre-existing
`docs/paramsweep_oracle240/` and `scripts/validate_sobol_oracle.py` artifacts
are intentionally outside the staged changes for this audit. This checkout's
progress is pushed to the `fast` remote; do not push unrelated changes to
`origin` without an explicit request.
