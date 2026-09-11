# `fast_v0.2` audit report

Evidence generation HEAD: `f87e969c96bff54bcdca8146f49d8b7a86ea501f`

Date: 2026-09-11

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
| Same-grid DN rel vs independent reference, PCHIP (default), six named points (median / max) | `2.93e-4 / 2.96e-4` | `<2e-4` | not met |
| Same-grid DN rel, legacy Simpson option, six named points (median / max) | `1.35e-3 / 1.24e-2` | `<2e-4` | not met |
| Formal default warm median / p95 | `4.93 / 5.47 ms` | stable `<4 ms` | not met |
| Formal six-point stability | `0` numerical failures, `3` explicit guards / 24 points | no silent fallback | pass |
| `eval_freqs` invariant | DN delta `0` on six probes | unchanged DN | pass |

The additional `cr0`, tilt, and Sobol oracle sample has no numerical failure;
the same-grid PCHIP DN relative error there is `2.76e-4–2.98e-4`, so the
parameter-space gate is also `NOT VERIFIED`.

The oracle-choice check is material: on the same default native grid, the
independent reference gives `DN=0.0022718753` at `z_tail=5` versus
`0.0022643136` at `z_tail=8` (relative delta `3.34e-3`). The formal fast
same-grid residual must therefore be interpreted as a combined tail/transfer
and frequency-integration residual; it is not valid to attribute all of it to
the frequency quadrature alone. That is why the Oracle C WKB anchor, which
holds the shared tail convention fixed, is reported alongside it.

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
or below `1.66e-4` against the Oracle C WKB anchor, so the largest remaining
fast-scale DN term is the same-grid reference residual (`2.93e-4`), which mixes
the frozen-tail/transfer systematic with the sparse-spectrum representation.
The fast-vs-reference spectrum p95/max (`~7.8e-4/~3.6e-3` at default) shows
that solver transfer error also contributes. The fresh 20-thread 25-repeat
matrix gives a default warm median of `4.93 ms`; high-T / stiff / high-kappa
remain the slowest regimes (`7.79 / 7.79 / 7.27 ms` warm median) and are the
next runtime target. The profiler's largest single components remain expansion
background (`0.80 ms`), tensor kernel (`0.98 ms`), and phase propagation
(`0.42 ms`), with no accepted further micro-optimization.

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
