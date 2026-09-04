# Accuracy

Status: current (honestly reported limits)
Date: 2026-09-03
Code version: see manifest `commit`

The precision anchor is the independent continuous-sigma reference
(`stiffgwpy_fast.reference`).  Accuracy is layered:

1. **Level 1 — spectrum.** signal/transition-region `Omega_GW` relative error
   and dex error.
2. **Level 2 — integrated physics.** `Delta_Neff` (the bolometric `DN_gw` integral).
3. **Level 3 — likelihood.** `Delta logL` from the spectrum on the likelihood bins.
4. **Level 4 — inference.** posterior parameter shift / sigma.

The two fast profiles vs the oracle (matched z8, 9 points):

| | plain-grid | transition-refine |
|---|---|---|
| spectrum rel (signal) median | 1.9e-2 | ~2.4e-4 |
| spectrum rel (signal) max | 7.0e-2 | 7.1e-4 |
| spectrum dex (signal) max | 8.2e-3 | 3.1e-4 |
| integrated `DN_gw` rel median | 9.1e-3 | 4.3e-4 |
| `Delta logL` (posterior bulk) | — | max 7.3e-3 |
| posterior `log10 r` shift | — | -0.0011 sigma |

## The honest limit on `Delta_Neff < 1e-4`

The integrated `DN_gw` relative error does **not** reach 1e-4 (median 4.3e-4).
The cause is not a tuning knob: the reference itself carries a ~3e-4
`z_tail`-frozen-tail sensitivity.  Measured at the default point
(`reference.oracle_variants`, rtol=1e-8, signal-band subset):

| oracle choice | `DN_gw` relative change |
|---|---|
| `z_tail` 7 → 8 | 4.2e-4 |
| `z_tail` 8 → 10 | 3.0e-4 |
| `z_tail` 14 (deep / no-tail) | *infeasible* (deep-subhorizon stiff) |

So production's residual is at the level of the oracle's own model choice, not a
solver defect.  This is reported as an honest bound, never gated away.

## Local error budget honesty

`estimate_local_error` distinguishes:

* `local-measured` — computed from this solve's telemetry
  (WKB handoff, frequency-grid error, quadrature Richardson, cancellation,
  self-consistency bracket);
* `calibrated-at-fiducial` — a measured default-point anchor scaled to the
  solve's settings (background model, transition, ODE phase, interpolation, tail);
* `uncertified` — the model carries no solve telemetry, so the returned budget is
  a conservative default, not a measurement.

The combined budget never claims to be a universal per-point error estimate when
its dominant systematic terms are fiducial-calibrated; it reports
`certification_status = certified-fiducial-calibrated` in that case.
