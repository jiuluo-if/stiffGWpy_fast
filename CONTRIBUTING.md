# Contributing

Install the development dependencies with `pip install .[dev]`. The normal
local gate is:

```bash
python scripts/validate_manifest.py
ruff check stiffgwpy_fast/__init__.py stiffgwpy_fast/cobaya/__init__.py stiffgwpy_fast/cobaya/stiffGW.py stiffgwpy_fast/_metrics.py stiffgwpy_fast/_resources.py stiffgwpy_fast/config.py stiffgwpy_fast/global_param.py stiffgwpy_fast/exact_background.py stiffgwpy_fast/freq_adaptive.py stiffgwpy_fast/reference.py scripts tests
python -m mypy
python -m pytest -q
python -m pytest -q -m slow
python -m build --sdist --wheel
python scripts/verify_distribution.py dist
python scripts/smoke_installed_wheel.py dist/stiffgwpy_fast-*.whl
```

The independent continuous-sigma reference is the precision anchor. LSODA is
for regression and runtime comparison only. Do not turn a `NOT VERIFIED` or
`FAIL` artifact into a passing claim by changing documentation alone; update
the executable validation and its provenance together.

Generated build/cache files, local credentials, and unfinished research
outputs are ignored. Keep PyPI tokens outside version control and inspect
`git status --short` before every commit. This checkout's progress is pushed to
the `fast` remote; do not push unrelated changes to `origin` without an
explicit request.
