# Numerical method (two fast profiles)

Status: current
Date: 2026-09-03
Code version: see manifest `commit`

The fast solver (`stiffgwpy.fast_sgwb`) uses a fixed-step Magnus-type
integration of the tensor-mode equation on a grid over `N`, then assembles the
spectrum and the bolometric integral.  Two user-facing profiles differ only in
how the hard features are treated.

## fast plain-grid

* fixed step `h` (`0.02`) as the expansion-grid spacing,
* a plain `construct` frequency grid at `freq_res = 1.0`,
* **no** kink-aware refinement (`transition_refine = False`),
* **no** phase-aware horizon-crossing sub-stepping (`phase_max = 0.0`),
* the deep-subhorizon hand-off to the analytic tail at `z_tail = 5.0`.

Cost: it skips the transition refine and the adaptive grid. After the
execution-layer JIT fix, the default point is `6.879 ms/point` warm median
(`7.642 ms` p95 at 4 threads; cold JIT measured separately). Its accuracy is
**not** certified (see `accuracy.md`); the speedup does not alter this status.

## fast transition-refine (production)

* kink-aware re-meshing: `exact_background.build_kink_refined_grid` puts the
  reheating kink inside a refined sub-step so `sigma(N)` never crosses a
  discontinuity on a spline/grid,
* `phase_max = 0.5` caps the per-sub-step phase increment `dTheta = e^z dh`
  around horizon crossing (adaptive Magnus sub-stepping with the `z~0` crossing
  band entered analytically),
* curvature-adaptive frequency grid (`freq_adaptive`) seeded from the
  grid-independent grid and refined where `|y''| h^2 / 8` of `log10 Omega_GW`
  exceeds the target dex,
* a `z_tail = 8.0` analytic frozen-tail (WKB/adiabatic) hand-off, with the
  per-mode adiabaticity defect returned as a verifiable local error estimate,
* point-local a-posteriori error budget (`estimate_local_error`) with 11
  categories, each tagged `local-measured` or `calibrated-at-fiducial`, plus a
  `certification_status` (`uncertified` when the solve carries no telemetry).

## Convergence (fast-vs-fast, default point)

Measured by `scripts/validate_two_modes.py --phase convergence`:

| profile | knob | `DN_gw` change |
|---|---|---|
| plain-grid | `h` 0.04 → 0.005 | ~1.4e-2 (non-monotonic; grid-kink bias dominates) |
| plain-grid | `z_tail` 5 → 10 | ~4e-3 |
| production | `h` 0.02 → 0.005 | ~2.5e-4 |
| production | `phase_max` 0.5 → 0.125 | ~2.4e-4 |
| production | `z_tail` 7 → 10 | ~3e-4 |
| production | `freq_res` 1 → 4 | ~7.7e-6 (converged) |

`production` is internally converged to ~2-3e-4; the residual is bound by the
`z_tail` frozen-tail term (shared with the oracle), not by the ODE/grid.
