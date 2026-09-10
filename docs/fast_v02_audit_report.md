# `fast_v0.2` audit report

Evidence generation HEAD: `5ccd02e3c3a1c4426f44e957bbe90426388bb544`  
Date: 2026-09-10  
Authoritative manifest: [`validation_manifest.json`](validation/validation_manifest.json)

## Executive result

The formal fast path is operational and stable on the audited sample, but is
not certified against the requested release gates. The current production
configuration remains `h=0.005`, `col_step=8`, `z_tail=5`, `phase_max=0.25`,
exact kink split, and `freq_grid=goal`.

## Accuracy and runtime

| Measure | Fresh result | Gate | Decision |
|---|---:|---:|---|
| Same-grid DN rel, Simpson, six named points (median / max) | `1.35e-3 / 1.24e-2` | `<2e-4` | not met |
| Same-grid DN rel, PCHIP, six named points (median / max) | `2.94e-4 / 2.96e-4` | `<2e-4` | not met; opt-in |
| Formal default warm median / p95 | `4.26 / 5.29 ms` | stable `<4 ms` | not met |
| Formal six-point stability | `0` numerical failures, `3` explicit guards / 24 points | no silent fallback | pass |
| `eval_freqs` invariant | DN delta `0` on six probes | unchanged DN | pass |

The additional `cr0`, tilt, and Sobol oracle sample has no numerical failure;
PCHIP DN relative error is `2.76e-4–2.98e-4`, so the parameter-space gate is
also `NOT VERIFIED`.

## Accepted changes

- exact reheating kink split and the gamma-corrected deep-tail matching remain
  enabled;
- native evaluation nodes are separated from bolometric integration support
  nodes;
- the cheap Simpson/trapezoid telemetry optimization was accepted after a
  fixed 50-repeat profiler showed `4.976 -> 4.544 ms` (`8.7%`) with unchanged
  output digests;
- outer reuse is retained for the current profile: four-case A/B has
  `false_safe_count=0` under DN `2e-4` and spectrum-max `1e-3`, although the
  default spectrum max delta is near the latter budget at `9.03e-4`.

## Rejected experiments

- global PCHIP promotion: consistently better than Simpson, but remains above
  the DN gate and has a low-T same-spectrum discrepancy of order `1e-2`;
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

The independent same-spectrum comparisons establish frequency integration and
sparse-spectrum representation as the leading recoverable DN error source,
but the fast-vs-reference spectrum p95/max (`~7.8e-4/~3.6e-3` at default)
shows that solver transfer error also contributes. The current profiler's
warm total is about `4.54 ms`; largest measured components are expansion
background (`0.80 ms`), tensor kernel (`0.98 ms`), and phase propagation
(`0.42 ms`), with no accepted further micro-optimization.

## Reproducible evidence

- `benchmark_head_matrix.json`
- `benchmark_fast_stability_head.json`
- `frequency_same_grid_reference_*.json`
- `benchmark_node_counts_head.json`
- `benchmark_eval_invariant_head.json`
- `benchmark_phase_candidate_head.json`
- `error_budget_probe_*_current.json`

The branch must remain `PARTIALLY VERIFIED` until a new oracle-backed
quadrature/solver candidate satisfies both the DN and runtime gates.
