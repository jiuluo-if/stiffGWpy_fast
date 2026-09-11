# Changelog

## Unreleased — 2026-09-11

- Switched the single formal fast mode's default `frequency_quadrature` to
  shape-preserving PCHIP.  Against the Oracle C WKB anchor the four named
  points move from `4.31e-4 / 1.24e-2 / 7.50e-4 / 1.29e-3` (Simpson) to
  `1.07e-5 / 1.66e-4 / 9.88e-6 / 1.15e-5` (PCHIP); `simpson` stays selectable
  and remains the embedded reference estimator.  The PCHIP hot path degrades a
  non-finite integrand to the established `nonfinite` abort/restore guard.
- Regenerated the HEAD validation evidence (runtime matrix, stability screen,
  same-grid reference, Oracle C anchors, error-budget probes, estimator
  coverage) for the new default and refreshed the validation manifest.
- Re-calibrated the `fast`/`ultra-fast` frequency-quadrature error budget from
  `1.0e-3` to `2.0e-4` (the measured PCHIP residual against the Oracle C WKB
  anchor is `1.07e-5`–`1.66e-4`), and refreshed `README.md`, `README_zh.md` and
  `docs/fast_v02_audit_report.md`.  The same-grid reference residual is
  `2.93e-4`; it is reported as a combined tail/transfer and quadrature term, so
  the branch stays `PARTIALLY VERIFIED`.

## Unreleased — 2026-09-10

- Separated native `eval_freqs` output nodes from bolometric integration support
  nodes, preserving the self-consistent `DN_gw` invariant across six probe
  regimes.
- Added Simpson/PCHIP embedded quadrature telemetry and current-HEAD evidence
  for frequency quadrature, node-count convergence, runtime, and fast-only
  stability. PCHIP and DN-driven refinement remain opt-in/rejected where the
  oracle evidence is not sufficient.
- Refreshed the release evidence for outer-reuse safety and recorded the
  rejected 89-node/PCHIP candidate; the manifest remains explicitly
  `PARTIALLY VERIFIED` because the DN and stable-runtime gates are not met.

## 0.2.1 — 2026-09-04

- Renamed the distribution metadata and Python import package to
  `stiffgwpy_fast`; PyPI normalizes the distribution display name to
  `stiffgwpy-fast`.
- Updated all Cobaya qualified names, resource paths, scripts, tests, CI
  commands, and bilingual installation instructions to use the new name.

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
