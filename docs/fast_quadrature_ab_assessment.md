# Fast frequency-quadrature A/B: Simpson versus PCHIP

## Motivation

`docs/fast_residual_decomposition_assessment.md` showed that fast's
propagation kernel is accurate to `~1e-5` while the default composite-Simpson
frequency quadrature contributes `4.2e-4..1.3e-2`.  This A/B quantifies both
sides of switching the default to PCHIP: the true DN accuracy it buys and the
warm runtime it costs.

## Protocol

`scripts/benchmark_fast_quadrature_ab.py` re-solves `SGWB_iter_fast` at the
current HEAD for each named point and quadrature, on the same native grid and
accuracy mode.  Timings alternate between the two quadratures inside every
repeat, flipping the order on odd repeats, so machine drift cannot favour one
method.  Resources: workers=1, Numba=2, BLAS=1; warmup=3, repeats=25.

## Results

| point | Simpson DN_gw | Simpson warm median | PCHIP DN_gw | PCHIP warm median | ratio | Simpson vs WKB | PCHIP vs WKB |
|---|---:|---:|---:|---:|---:|---:|---:|
| default | 2.2626969518e-03 | 6.76 ms | 2.2636474151e-03 | 8.83 ms | 1.307 | 4.31e-04 | **1.07e-05** |
| lowT | 5.2446047710e-08 | 6.86 ms | 5.3114602834e-08 | 8.92 ms | 1.300 | 1.24e-02 | **1.66e-04** |
| highT | 5.6391339588e-02 | 10.18 ms | 5.6433084271e-02 | 12.26 ms | 1.204 | 7.50e-04 | **9.88e-06** |
| stiff | 1.4965149678e-02 | 10.49 ms | 1.4984335009e-02 | 13.41 ms | 1.278 | 1.29e-03 | **1.15e-05** |

## Findings

- Accuracy: PCHIP cuts the true DN error against the Oracle C WKB anchor by one
  to two orders of magnitude — `4.31e-04 -> 1.07e-05` (default),
  `7.50e-04 -> 9.88e-06` (highT), `1.29e-03 -> 1.15e-05` (stiff),
  `1.24e-02 -> 1.66e-04` (lowT).  Every point lands inside the `2e-04` release
  gate.
- Runtime: the scipy-backed PCHIP path costs `+20%..+31%`
  (`8.83/6.76 = 1.307` default), i.e. it *fails* the pre-registered `< 10%`
  budget for an accuracy win.
- Micro-benchmarks explain the cost: `PchipInterpolator(...).integrate()` is
  `0.141 ms` per call (called ~3 times per solve: two outer-loop `g2_last`
  plus `g2c[-1]`) and `estimate_frequency_quadrature_local(..., 'pchip')` is
  `1.39 ms` per call (its per-interval `spline.integrate` loop).  Together that
  reproduces the observed `+2.1 ms`.
- The cost is entirely implementation, not mathematics: PCHIP integration is a
  *linear* functional of the node values, so the global integral and all
  per-interval integrals can be precomputed once per grid as a weight vector
  and a weight matrix, after which each solve only needs `w @ y`.

## Decision

**ACCEPTED as a measurement.**  PCHIP is the accuracy fix; the scipy
implementation is not yet fast enough to be the default.  No solver default is
changed by this phase.

## Next experiment

Hypothesis: precomputing the PCHIP integration weights (global vector plus
per-interval matrix) and using them for both the DN integral and the local
quadrature estimator removes the scipy overhead, so PCHIP can become the
default within the runtime budget.

Acceptance criteria (fixed before running):

- PCHIP warm runtime median within `10%` of the current Simpson default;
- identical DN_gw to the current scipy PCHIP path (bitwise or `< 1e-12` rel);
- true DN rel vs the Oracle C WKB anchor `< 2e-04`;
- `spectrum max` unchanged or better; no new failure; determinism replay identical.

## Artifacts

- `docs/fast_quadrature_ab.json`
- `scripts/benchmark_fast_quadrature_ab.py`
