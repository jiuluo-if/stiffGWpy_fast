# Test coverage

The default suite protects the installed package: public API and configuration, Fast solver behavior, independent reference behavior, physics limits, numerical integration, compatibility, and fixed accuracy regressions.

The independent numerical samples required by two regression tests are stored in `tests/data/`. Small reference helpers used by core solver tests are isolated in `tests/support/`; experiment drivers and bulk evidence remain in the local research archive.

CI runs the compatibility smoke across supported Python versions, the canonical suite once, static checks, source/wheel builds, an installed-wheel smoke, and the optional Cobaya integration job.
