# Reproducibility

Install the project with its development extra in a clean Python environment:

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
ruff check stiffgwpy_fast scripts tests
python -m mypy
python -m build --sdist --wheel
python scripts/verify_distribution.py dist
python scripts/smoke_installed_wheel.py dist/stiffgwpy_fast-*.whl
```

The GitHub Actions workflow also checks supported Python versions and the optional Cobaya adapter. Its thread limits are declared in `.github/workflows/ci.yml`.

For numerical comparisons, record the engine, full solver configuration, input parameters, source commit, frequency grid, and independent reference convention. The committed samples in `tests/data/` support narrowly scoped regression checks. Keep larger scans and raw experiment records in the local research archive; do not treat a small set of tested points as certification over the full parameter space.
