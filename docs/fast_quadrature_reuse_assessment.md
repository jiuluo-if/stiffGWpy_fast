# Fast PCHIP quadrature: single-fit sharing and vectorized Simpson baseline

## Motivation

`docs/fast_quadrature_ab_assessment.md` established that PCHIP is the accuracy
fix (true `DN_gw` error against the Oracle C WKB anchor `1.07e-05..1.66e-04`,
all four named points inside the `2e-04` release gate) but that the scipy-backed
implementation cost `+20%..+31%` warm runtime, *failing* the pre-registered
`< 10%` budget for an accuracy win.  This phase attacks that cost.

## Hypothesis (and a falsified predecessor)

The predecessor hypothesis — "PCHIP integration is a linear functional of the
node values, so it can be precomputed once per grid as a global weight vector
plus a per-interval weight matrix" — is **FALSIFIED**.  The Fritsch-Carlson
slopes are a *nonlinear* function of the node values.  A unit-matrix probe
(`PchipInterpolator(x, np.eye(n), axis=0)`) produces weights that agree with
column-wise interval integration, yet `weights @ y` disagrees with
`PchipInterpolator(x, y).integrate(...)` by ~20x.  No fixed weight vector
exists.  The corrected hypothesis is **single-fit sharing**: build the
interpolant once and reuse it.

## Method

- `pchip_integral_breakdown(freqs, integrand)` builds one
  `PchipInterpolator` plus one `antiderivative()` and returns both the
  full-range integral and every per-interval integral of that same fit.
- `estimate_frequency_quadrature_local(..., candidate_intervals=...)` accepts
  those per-interval integrals instead of rebuilding the spline and looping
  `spline.integrate(left, right)` over every native interval.
- The non-overlapping 3-point Simpson baseline in the local estimator is
  vectorized with the same non-uniform-spacing formula scipy uses; the
  replacement is bit-identical (measured max abs difference `0.0`).
- `_SGWB_iter_fast_impl` performs exactly one breakdown per outer iteration and
  reuses it for the convergence test `g2_last`, the post-loop `g2c[-1]` and the
  local estimator, removing the two duplicated scipy fits.

## Protocol

`scripts/benchmark_fast_quadrature_ab.py --repeats 50`, workers=1, Numba=2,
BLAS=1, warmup=3.  To remove machine drift the two code states were measured
back-to-back in the same session: first the parent-commit `fast_sgwb.py`
(restored with `git checkout --`), then the working-tree change.  Timings
alternate the two quadratures inside each repeat and flip the order on odd
repeats.

## Results (same-session paired A/B)

| point | Simpson median before | Simpson median after | PCHIP median before | PCHIP median after | PCHIP extra before | PCHIP extra after | ratio before | ratio after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| default | 6.02 ms | 5.92 ms | 8.22 ms | 6.90 ms | 2.20 ms | 0.98 ms | 1.366 | 1.165 |
| lowT | 6.27 ms | 6.21 ms | 8.19 ms | 6.80 ms | 1.92 ms | 0.59 ms | 1.307 | 1.095 |
| highT | 8.87 ms | 9.84 ms | 10.88 ms | 10.56 ms | 2.01 ms | 0.72 ms | 1.228 | 1.073 |
| stiff | 9.54 ms | 9.37 ms | 11.83 ms | 10.21 ms | 2.29 ms | 0.84 ms | 1.240 | 1.089 |

`DN_gw` is **bit-identical** before and after at every point (relative
difference `0.0`): default `2.2626969518100204e-03` (Simpson) /
`2.2636474150855798e-03` (PCHIP), lowT `5.2446047710354755e-08` /
`5.3114602833937853e-08`, highT `5.6391339588221887e-02` /
`5.6433084271417577e-02`, stiff `1.4965149677851021e-02` /
`1.4984335009443728e-02`.  The `vs WKB` columns are unchanged
(`4.31e-04 / 1.24e-02 / 7.50e-04 / 1.29e-03` for Simpson and
`1.07e-05 / 1.66e-04 / 9.88e-06 / 1.15e-05` for PCHIP).

## Micro-cost decomposition (median of 300 calls, 76-node native grid)

| operation | before | after |
|---|---:|---:|
| `pchip_integral_breakdown` (1 fit + antiderivative + interval integrals) | - | 0.150 ms |
| PCHIP fit in the outer loop (2 per solve) | 2 x 0.143 ms | 2 x 0.150 ms (shared) |
| post-loop `g2c[-1]` fit | 0.143 ms | 0 (reused) |
| `estimate_frequency_quadrature_local(..., 'pchip')` | 1.390 ms | 0.290 ms |
| -- of which the Simpson panel baseline | 0.674 ms (37 scipy calls) | 0.027 ms (vectorized, bit-identical) |

## Findings

- The optimization is a strict Pareto improvement for the PCHIP path:
  bit-identical `DN_gw`, `2.0-2.3 ms -> 0.6-1.0 ms` of PCHIP-specific overhead
  (55-70% less), and the PCHIP/Simpson warm median ratio drops from
  `1.23-1.37` to `1.07-1.17`.
- The pre-registered `< 1.10` budget for switching the default is nonetheless
  **not met**: the default point measured `1.165` in the paired session
  (`1.093`/`1.123` in two earlier sessions), i.e. it straddles the budget.
- A 16-thread probe (`NUMBA_NUM_THREADS=16`, the production thread scale) shows
  the same serial overhead is relatively *more* expensive when the parallel
  part shrinks: ratio `1.154/1.134/1.137/1.101` for default/lowT/highT/stiff.
- The remaining PCHIP-specific cost is the scipy `PchipInterpolator`
  construction itself (`0.15 ms` x 2 per solve) plus the estimator's
  per-panel allocation loops (~`0.26 ms`).  Both are implementation costs, not
  mathematics.

## Decision

**ACCEPTED: the PCHIP single-fit sharing and vectorized Simpson baseline.**
Identical observables, ~60% less PCHIP-specific overhead; the opt-in PCHIP path
strictly dominates its previous self, so this is a speed win under principle 12.

**REJECTED for this phase: switching the default `frequency_quadrature` to
`pchip`.**  The `< 10%` runtime budget is not robustly satisfied.  The default
stays `simpson`; PCHIP remains the explicit accuracy path.

## Next experiment

Hypothesis: replacing the scipy PCHIP fit with a NumPy vectorized
slope+segment-integral kernel (validated against scipy) and vectorizing the
estimator's allocation loops removes the last duplicated implementation cost,
so PCHIP lands within the runtime budget at both the 2-thread and the 16-thread
production scale.

Acceptance criteria (fixed before running):

- every vectorized integral agrees with the scipy PCHIP reference to `< 1e-12`
  relative (and the existing `test_pchip_frequency_quadrature_is_opt_in`
  invariant, which pins `m.DN_gw[-1]` to `integrate_frequency_pchip`, still
  holds);
- PCHIP/Simpson warm median ratio `< 1.10` at 2 threads **and** at 16 threads;
- `DN_gw` unchanged versus the current scipy PCHIP path (`< 1e-12` rel);
- DN vs WKB `< 2e-04`; `spectrum max` not degraded; no new failure; determinism
  replay identical.

## Artifacts

- `docs/fast_quadrature_reuse_ab.json` (after, 50 repeats)
- `docs/fast_quadrature_reuse_ab_before.json` (parent-commit code, 50 repeats, same session)
- `docs/fast_quadrature_ab.json` (original 25-repeat measurement, previous session)
- `scripts/benchmark_fast_quadrature_ab.py`
- `tests/test_pchip_integral_breakdown.py`
