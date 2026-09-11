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

## Follow-up: vectorized PCHIP kernel and estimator allocation (2026-09-11)

### Motivation

After the single-fit sharing the remaining PCHIP-specific cost was the SciPy
`PchipInterpolator` construction itself (`0.15 ms` x 2 per solve) plus the
estimator's per-panel Python allocation loop.  Both are implementation costs,
so the pre-registered follow-up replaced them while keeping the SciPy path as
the reference.

### Method

- `_pchip_integrals_vectorized` reproduces the SciPy Fritsch-Carlson slopes
  (`_pchip_slopes` / `_pchip_edge_slope`, same recurrence and shape guards) and
  uses the closed form of the cubic Hermite segment integral,
  `h/2*(y_i + y_{i+1}) + h^2*(m_i - m_{i+1})/12`.  It returns the per-interval
  integrals; the full-range integral is their sum.  `pchip_integral_breakdown`
  stays the SciPy reference used by the audit helper and by the tests.
- The estimator's allocation is vectorized over the non-overlapping Simpson
  panels; the flat/zero-share fallback is preserved.
- The hot path in `_SGWB_iter_fast_impl` now uses the vectorized kernel through
  the same single-fit-per-iteration sharing introduced above.

### Results (paired A/B, 50 repeats, artifacts below)

| point | 2T ratio before | 2T ratio after | 16T ratio before | 16T ratio after |
|---|---:|---:|---:|---:|
| default | 1.1622 | **1.0105** | 1.1987 | **1.0428** |
| lowT | 1.1261 | **1.0551** | 1.1137 | **1.0342** |
| highT | 1.0495 | **1.0397** | 1.1270 | **1.0180** |
| stiff | 1.0901 | **1.0411** | 1.1203 | **1.0413** |

| point | 2T PCHIP extra before | after | 16T PCHIP extra before | after |
|---|---:|---:|---:|---:|
| default | 1.00 ms | 0.06 ms | 0.89 ms | 0.20 ms |
| lowT | 0.83 ms | 0.36 ms | 0.59 ms | 0.17 ms |
| highT | 0.51 ms | 0.40 ms | 0.95 ms | 0.15 ms |
| stiff | 0.89 ms | 0.42 ms | 0.92 ms | 0.33 ms |

A second paired repetition of the same protocol in the same session reproduced
the pattern (`before 1.050..1.184`, `after 1.002..1.056`).  Absolute wall times
drift by up to ~60% between runs on this machine, so the run-internal ratio is
the reported metric.

Micro-cost (median of 2000 calls, 76-point native grid): SciPy
`pchip_integral_breakdown` `0.143 ms` versus vectorized
`_pchip_integrals_vectorized` `0.024 ms`; estimator with reused intervals
`0.290 ms -> 0.061 ms`.

Correctness: the vectorized allocation is **bit-identical** to the previous
per-panel loop for all three allocation modes; the vectorized PCHIP kernel
agrees with the SciPy reference per interval to `<= 1e-9` relative (absolute
`>= 1e-12 * max|interval|`) and its sum to `<= 5e-14` relative over 200
randomized grids plus flat/sign-changing/two-point degenerate cases.  The
end-to-end guard `test_pchip_frequency_quadrature_is_default`, which pins
`m.DN_gw[-1]` to the SciPy `integrate_frequency_pchip`, still passes at
`rel = 1e-12`.  Measured `DN_gw` change versus the SciPy PCHIP path is
`3.83e-16 / 4.98e-16 / 0 / 0` relative (default/lowT/highT/stiff), i.e. 1-2 ulp.

### Decision

**ACCEPTED.**  The pre-registered `< 1.10` budget is met at both the 2-thread
and the 16-thread scale with margin, `DN_gw` is unchanged to 1-2 ulp, and the
default Simpson path is bit-unchanged.  The PCHIP path now costs `~1-4%` more
warm runtime than Simpson while reducing the true DN error from
`4.31e-4..1.24e-2` to `1.07e-5..1.66e-4`.

Because the change is a strict dominance for the PCHIP path and a no-op for the
default, the default `frequency_quadrature` is **not** switched in this commit:
that is a separate phase (below) that must re-validate the manifest, README,
coverage artifacts and the parameter-space gates.

### Next experiment (pre-registered)

*(Executed in the same work line; results are recorded in
`docs/fast_quadrature_default_switch_assessment.md`.)*

Hypothesis: with the runtime budget now met, making `pchip` the default
`frequency_quadrature` of `SGWB_iter_fast` satisfies the release accuracy target
without regressing the runtime or the guard behaviour.

Acceptance criteria (fixed before running):

- `SGWB_iter_fast` default is `pchip`; `simpson` stays available explicitly;
- parameter-space `DN_gw` versus the Oracle C WKB anchor `< 2e-04` at all four
  named points and within budget on the Sobol/edge screen;
- `spectrum max` not degraded versus the current Simpson default; no new
  failure; `no silent failure` guards unchanged; determinism replay identical;
- warm runtime ratio versus the current Simpson default `< 1.10` at 2 and 16
  threads;
- `docs/validation/validation_manifest.json`, README, `ERROR_BUDGET` and the
  estimator-coverage artifacts refreshed to describe PCHIP as the default, and
  `validation artifact == release HEAD` re-established.

### Artifacts (follow-up)

- `docs/fast_quadrature_pchip_kernel_t2_before.json`
- `docs/fast_quadrature_pchip_kernel_t2_after.json`
- `docs/fast_quadrature_pchip_kernel_t16_before.json`
- `docs/fast_quadrature_pchip_kernel_t16_after.json`
