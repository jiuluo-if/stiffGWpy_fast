# Fast quadrature: switching the default to PCHIP

## Motivation

The residual decomposition (`docs/fast_residual_decomposition_assessment.md`)
showed the composite Simpson rule was the dominant `DN_gw` error term of the
formal fast path: against the Oracle C WKB anchor the four named points sat at
`4.31e-4 / 1.24e-2 / 7.50e-4 / 1.29e-3` (default/lowT/highT/stiff).  PCHIP
measured `1.07e-5 / 1.66e-4 / 9.88e-6 / 1.15e-5`, i.e. every point inside the
`2e-4` release gate.  The earlier phases removed the PCHIP-specific runtime
premium (`docs/fast_quadrature_reuse_assessment.md`), so the quadrature could be
promoted from opt-in to default.

## Hypothesis

Making `pchip` the default `frequency_quadrature` of `SGWB_iter_fast` meets the
release accuracy target without regressing runtime, guards or determinism.

Acceptance criteria (fixed before running, as pre-registered in
`docs/fast_quadrature_reuse_assessment.md`):

- `SGWB_iter_fast` default is `pchip`; `simpson` stays available explicitly.
- Parameter-space `DN_gw` versus the Oracle C WKB anchor `< 2e-04` at all four
  named points and within budget on the Sobol/edge screen.
- `spectrum max` not degraded versus the current Simpson default; no new
  failure; guards unchanged; determinism replay identical.
- Warm runtime ratio versus the current Simpson default `< 1.10` at 2 and 16
  threads.
- `docs/validation/validation_manifest.json`, README, `ERROR_BUDGET` and the
  estimator-coverage artifacts refreshed to describe PCHIP as the default, and
  `validation artifact == release HEAD` re-established.

## Method

- `_SGWB_iter_fast_impl` and `SGWB_iter_fast` default `frequency_quadrature` is
  `'pchip'`; the audit methods stay selectable.
- The PCHIP hot path checks the integrand for finiteness before the vectorized
  kernel.  A non-finite integrand is *not* fed to the fail-loud helper; the NaN
  is propagated to the shared `math.isfinite(DN_gw_new)` guard, which reproduces
  the established `fast_failure_reason='nonfinite'` abort/restore semantics of
  the Simpson path.
- Oracle anchors (`docs/oracle_c_wkb_*.json`) and the same-grid reference
  artifacts were regenerated with the new default so the release evidence and
  the HEAD behaviour cannot drift.

## Results

### Accuracy versus the Oracle C WKB anchor (new default)

| point | Simpson (previous default) | PCHIP (new default) |
|---|---:|---:|
| default | 4.31e-4 | 1.07e-5 |
| lowT | 1.24e-2 | 1.66e-4 |
| highT | 7.50e-4 | 9.88e-6 |
| stiff | 1.29e-3 | 1.15e-5 |

### Runtime

Same-session paired A/B of the two defaults, 50 repeats, workers=1, BLAS=1
(`docs/fast_quadrature_default_ab_t{2,16}.json`); the ratio is the
`pchip / simpson` warm median:

| point | 2-thread ratio | 16-thread ratio |
|---|---:|---:|
| default | 1.019 | 1.069 |
| lowT | 1.035 | 1.028 |
| highT | 1.034 | 0.987 |
| stiff | 1.024 | 1.017 |

The earlier kernel-phase A/B
(`docs/fast_quadrature_pchip_kernel_t{2,16}_{before,after}.json`) measured the
same comparison before and after the vectorized PCHIP kernel.

### Tests

`python -m pytest -q` -> `159 passed, 6 deselected`.  The renamed guards are
`test_pchip_frequency_quadrature_is_default` (default is PCHIP and the DN_gw
value matches the SciPy `integrate_frequency_pchip`) and
`test_simpson_frequency_quadrature_remains_explicit` (Simpson stays selectable
and yields a different, finite DN).

## Decision

ACCEPT: the default `frequency_quadrature` of the single formal fast mode is
`pchip`.  See the artifact list below for the regenerated release evidence.

`ERROR_BUDGET['fast'|'ultra-fast'].quadrature` was re-calibrated from `1.0e-3`
to `2.0e-4` to match the measured PCHIP range, and `README.md`, `README_zh.md`,
`CHANGELOG.md` and `docs/fast_v02_audit_report.md` were refreshed with the new
default, runtime matrix and accuracy numbers.

## Not yet met

The `2e-4` gate is now met against the Oracle C WKB anchor at the four named
points, but not against the independent same-grid reference (`2.93e-4` median,
`2.76e-4`–`2.98e-4` over the edge/Sobol sample), because that residual also
carries the frozen-tail/transfer systematic.  The runtime gate is still open:
the default warm median is `4.93 ms` and high-T / stiff / high-kappa stay at
`7.79 / 7.79 / 7.27 ms` on the 20-thread 25-repeat matrix.  Next: decompose the
high-T/stiff/high-kappa runtime, and run a nested native-frequency quadrature
experiment (real extra frequency solves on the top intervals) to decide whether
the same-grid residual is quadrature- or transfer-dominated.

## Artifacts

- `docs/benchmark_head_matrix.json`
- `docs/benchmark_fast_stability_head.json`
- `docs/frequency_same_grid_reference_*.json`
- `docs/oracle_c_wkb_*.json`
- `docs/fast_true_error.json`
