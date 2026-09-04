# Changelog

## 0.2.0 — 2026-09-04

- Added immutable per-call fast-solver configuration with `fast` plain-grid as
  the default high-level profile; `production` remains an explicit precision mode.
- Added package-resource loading and an installed-wheel smoke test covering core,
  Cobaya, LIGO, and PTA data.
- Added Python-version CI coverage, a scheduled/on-demand slow numerical job,
  coverage reporting, and validation-manifest schema checks.
- Fixed optional-Cobaya test collection and capped adapter thread requests to
  the runner's available Numba budget.
- Fixed the installed-wheel smoke probe so optional Cobaya resources are
  checked only when Cobaya is installed.
- Kept the independent continuous-sigma reference as the precision anchor;
  existing `NOT VERIFIED` and `FAIL` validation results remain unchanged.

Verification for this release line is documented in
`docs/engineering_audit.md` and `docs/reproducibility.md`.
