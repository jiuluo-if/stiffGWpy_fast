# Contributing

Install the development dependencies with `pip install '.[dev]'`, then run the
standard local checks:

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

The independent continuous-sigma reference is the precision anchor; LSODA is
for regression checks and runtime comparisons. A `NOT VERIFIED` or `FAIL`
artifact must not be presented as passing based on a documentation change.
Update and rerun the executable validation, and record its provenance.

Build and cache files, local credentials, and unfinished research outputs are
ignored. Keep PyPI tokens out of version control and check `git status --short`
before each commit. Push this checkout to the `fast` remote; do not push
unrelated changes to `origin` without an explicit request.
