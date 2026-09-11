# `fast_v0.2` audit report

Evidence generation HEAD: `a7b0ea9dae45c6a04c9a5760a2192cc21dc4a930`

Date: 2026-09-12

Authoritative manifest: [`validation_manifest.json`](validation/validation_manifest.json)

## Executive result

The formal fast path is operational and stable on the audited sample, but is
not certified against the requested release gates. The current production
configuration remains `h=0.005`, `col_step=8`, `z_tail=5`, `phase_max=0.25`,
exact kink split, `freq_grid=goal`, and (as of this HEAD) the default
`frequency_quadrature=pchip`.

## Accuracy and runtime

| Measure | Fresh result | Gate | Decision |
|---|---:|---:|---|
| Fast vs Oracle C WKB, default PCHIP, four named points (median / max) | `1.1e-5 / 1.66e-4` | `<2e-4` | pass |
| Same-grid DN rel vs independent reference at its own `z_tail=8`, PCHIP (default), six named points (median / max) | `2.93e-4 / 2.96e-4` | `<2e-4` | not met, oracle-limited (see `z_tail` attribution) |
| Same-grid DN rel vs Oracle A deepened to `z_tail=10`, PCHIP, default point | `4.96e-6` | `<2e-4` | pass |
| Nested native-frequency refinement, PCHIP, six named points, `N=4/8/12` (absolute / relative) | `1.7e-14-6.3e-10` / `7.4e-12-2.8e-9` | native grid converged | pass |
| Same-grid DN rel, legacy Simpson option, six named points (median / max) | `1.35e-3 / 1.24e-2` | `<2e-4` | not met |
| Formal default warm median / p95 | `4.77 / 5.58 ms` | stable `<4 ms` | not met |
| Formal six-point stability | `0` numerical failures, `3` explicit guards / 24 points | no silent fallback | pass |
| `eval_freqs` invariant | DN delta `0` on six probes | unchanged DN | pass |

The additional `cr0`, tilt, and Sobol oracle sample has no numerical failure;
the same-grid PCHIP DN relative error there is `2.76e-4–2.98e-4`, so the
parameter-space gate is also `NOT VERIFIED`.

The oracle-choice check is material and is now resolved by moving the *oracle*,
not fast: on the same default native grid the independent reference gives
`DN=0.0022718753` at `z_tail=5`, `0.0022643136` at `z_tail=8`, and
`0.0022636586` at `z_tail=10`, while fast (`z_tail=5`) gives `0.0022636474`.
The reference's own frozen-handoff defect `|1.5*sigma-1|/exp(z)` is `3.4e-4`
at `z=8` and `4.5e-5` at `z=10`, the same order as its distance from fast
(`2.944e-4 -> 4.956e-6`). At `z=10` the Oracle C anchor (Prüfer
amplitude-phase, `2.2636592187e-3`) and Oracle A agree to `2.4e-7`, and fast
sits `4.96e-6`/`5.21e-6` from them. The six-point `2.93e-4` row above therefore
measures the audited oracle's frozen-tail convention rather than fast's
frequency integration; the defensible default-point `DN_gw` residual against
tail-converged oracles is `~5e-6`.  Deepening the oracle is a one-point
attribution tool, not a per-point validation method: `z_tail=8` costs `94 s`,
`z_tail=10` costs `579 s`, and `z_tail=14` did not converge within `20 min`
(more modes stop triggering the analytic handoff and are integrated to today).

## Accepted changes

- exact reheating kink split and the gamma-corrected deep-tail matching remain
  enabled;
- native evaluation nodes are separated from bolometric integration support
  nodes;
- the single formal fast default `frequency_quadrature` moved from Simpson to
  shape-preserving PCHIP: against the Oracle C WKB anchor the four named points
  go from `4.31e-4 / 1.24e-2 / 7.50e-4 / 1.29e-3` to
  `1.07e-5 / 1.66e-4 / 9.88e-6 / 1.15e-5`, and the paired warm-runtime ratio
  stays `<1.10` at 2 and 16 threads (`docs/fast_quadrature_default_ab_t{2,16}.json`);
  `ERROR_BUDGET['fast'].quadrature` was re-calibrated from `1.0e-3` to `2.0e-4`
  to match that measurement;
- the cheap Simpson/trapezoid telemetry optimization was accepted after a
  fixed 50-repeat profiler showed `4.976 -> 4.544 ms` (`8.7%`) with unchanged
  output digests;
- the fast Python preparation layer was de-duplicated (goal-grid
  sort/unique collapse, per-`N` invariant hoisting, redundant `Ogw`/`Oj`/`Opgw`
  zero-fills, single `derived_param` binding): all six probe regimes keep
  bit-identical output digests (648 checked comparisons) and the paired
  A,B,B,A warm runtime improves `5.4%` (default, 20 threads) and `8.2%`
  (default, 2 threads), with no regime regressing beyond the `2%` noise margin
  (`docs/fast_pyoverhead_ab_symmetric.json`);
- outer reuse is retained for the current profile: four-case A/B has
  `false_safe_count=0` under DN `2e-4` and spectrum-max `1e-3`, although the
  default spectrum max delta is near the latter budget at `9.03e-4`.

## Rejected experiments

- global PCHIP promotion (2026-09-10 attempt): consistently better than
  Simpson, but the scipy PCHIP path cost `1.20–1.31x` warm runtime, so it was
  not promoted then. Retried and accepted on 2026-09-11 after the vectorized
  PCHIP kernel plus the shared single-fit estimator brought the paired ratio
  back inside `<1.10`;
- Chebyshev, fixed cubic/log-log alternatives, and Gauss-over-PCHIP: no
  uniform parameter-space advantage;
- DN-driven midpoint refinement and blind 76/80/90/110 seed increases:
  convergence is non-monotonic;
- 89-node (`seed_n=78`) plus PCHIP: default independent oracle remains
  `2.952e-4` DN relative error;
- phase caps `0.35` and `0.5`: no stable >5% speed gain across six points;
- Numba literal specialization of `assemble`: warm kernel regressed to about
  `150 ms`, so the source was fully reverted.

## Remaining dominant errors and hotspots

With PCHIP as the default, the recoverable frequency-quadrature error is now at
or below `1.66e-4` against the Oracle C WKB anchor, and the nested-frequency
refinement bounds the native-grid discretisation term at `<=6.3e-10` absolute
over six regimes.  The earlier headline same-grid residual (`2.93e-4`) is
dominated by the audited oracle's own `z_tail=8` frozen-handoff defect; against
tail-converged oracles the default-point DN residual is `~5e-6`.  The remaining
fast-side terms are the `z_tail=5` frozen-tail convention (shared with every
oracle) and the sparse-spectrum transfer representation.
The fast-vs-reference spectrum p95/max (`~7.8e-4/~3.6e-3` at default) shows
that solver transfer error also contributes. The fresh 20-thread 25-repeat
matrix gives a default warm median of `4.77 ms`; high-T / stiff / high-kappa
remain the slowest regimes (`6.91 / 6.81 / 6.97 ms` warm median) and are the
next runtime target. The 20-thread stage breakdown at the default point is led
by the expansion background (`1.04 ms`), the tensor solve kernel (`0.81 ms`),
and the split phase/S2 propagation (`0.63 ms`). The accepted change on this
HEAD is confined to the Python preparation layer (goal-grid de-duplication,
hoisting the per-`N` invariants, skipping redundant zero-fills, and binding
repeated `derived_param` evaluations), so it removes overhead without touching
the kernel.

## Reproducible evidence

- `benchmark_head_matrix.json`
- `benchmark_fast_stability_head.json`
- `frequency_same_grid_reference_*.json`
- `oracle_c_wkb_*.json`
- `fast_true_error.json`
- `fast_residual_decomposition.json`
- `benchmark_node_counts_head.json`
- `benchmark_eval_invariant_head.json`
- `benchmark_phase_candidate_head.json`
- `error_budget_probe_*_current.json`
- `quadrature_estimator_coverage_head.json`
- `fast_pyoverhead_ab_symmetric.json`
- `nested_frequency_quadrature_head.json`
- `nested_frequency_quadrature_simpson_control.json`
- `oracle_a_same_grid_z10_default.json`

The branch must remain `PARTIALLY VERIFIED` until the parameter-space same-grid
DN gate and the stable `<4 ms` runtime gate are both met.

## Honest limits of this evidence

The legacy `fast_transition_refine` / production artifacts
(`docs/paramsweep_z8*`, `docs/paramsweep_ref/`) were generated at their recorded
commits, before the `eval_freqs` separation and before PCHIP became the default
frequency quadrature.  The manifest now discloses them as conservative legacy
evidence rather than current-HEAD measurements.  Their spectrum bands are set by
the ODE and frequency grid, not by the bolometric quadrature, so only the
`DN_gw` aggregates are expected to move; re-running that suite is a
pre-registered follow-up, not a silent gap.

The `z_tail` attribution above is a single-point (default) measurement made
with `DN_eff` frozen at fast's self-consistent value.  The other five named
regimes still gate against the `z_tail=8` oracle, so the six-point `not met`
row is retained as an oracle-limited status instead of being reinterpreted per
point.  The nested-refinement check is likewise a sensitivity bound: it shows
the native grid is not the limiter, but it cannot certify an absolute error by
itself because `E_nested` is false-safe by 4-7 orders against the oracle
residual.
