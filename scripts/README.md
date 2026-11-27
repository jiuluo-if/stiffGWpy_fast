# Repository scripts

- `verify_distribution.py` checks that wheel and source archives exclude repository-only files.
- `smoke_installed_wheel.py` installs a wheel in an isolated target and checks imports and packaged resources.
- `_resource_budget.py` and the numerical reference helpers imported by `test_fast_sgwb.py`, `test_reference.py`, and `test_wkb3_fast_transfer.py` are retained as test support. They do not enter the installed distribution.

Bulk scans, candidate prototypes, plotting, posterior workflows, and standalone benchmark drivers remain in the local research archive.
