# Numerical method

Status: current user-facing solver contract

Date: 2026-09-24

Code version: current settings are defined by
`stiffgwpy_fast.fast_sgwb.ACCURACY_MODES['fast']`; run evidence remains bound to
the commit recorded in each artifact.

## User-facing engines

`LCDM_SG.SGWB_iter()` defaults to `engine='fast'` and the single formal
goal-kink-hybrid profile. The high-level engine choices are:

| Engine | Role | Numerical path |
|---|---|---|
| `fast` | Default practical solver and MCMC hot path | Fixed-step Magnus-style tensor propagation, Numba kernels, exact reheating-kink split, goal frequency grid, PCHIP bolometric integration |
| `reference` | Independent precision comparison | Continuous `sigma(N)`, adaptive DOP853 per frequency, error-estimated frequency quadrature, analytic deep-tail handoff |
| `lsoda` | Legacy regression path and explicit fallback | Original adaptive SciPy LSODA per-frequency integration |

`production`, `transition_refine`, `plain-grid`, and other names may still
appear in old reports or internal validation code. They are historical or
validation-only labels, not additional user-facing profiles. Only `fast` is
listed in `USER_FAST_PROFILES`.

## Formal `fast` preset

The current preset uses `h=0.005`, `col_step=8`, `z_tail=5`, `freq_res=1`,
outer tolerance `1e-6`, `phase_max=0.25`, exact kink splitting, and the
goal-oriented frequency grid. Frequency quadrature defaults to shape-preserving
PCHIP; Simpson remains explicitly selectable for audits. The preset's thread
budget is constrained by available runtime resources and may be overridden by
an explicit caller setting.

The goal grid reserves nodes around relevant spectral features and preserves
`eval_freqs` as native solve nodes. Evaluation nodes are separate from the
support nodes used for the bolometric integral, so adding likelihood bins does
not change `DN_gw` through the integration grid.

During tensor propagation, `phase_max` caps the phase advance near horizon
crossing. At `z_tail`, the deep-subhorizon evolution is handed to an analytic
WKB envelope. The point-local error budget exposes measured, fiducial-calibrated,
and uncertified components through `estimate_local_error`; it is not a
universal parameter-space certificate.

## Independent reference and LSODA

The `reference` engine evolves a continuous background with the reheating kink
as an exact breakpoint and uses adaptive DOP853 for tensor modes. Its own
frozen-tail convention contributes a measurable systematic, so comparisons
must name the reference `z_tail`, grid, tolerances, and sampled parameter
points. A deep/no-tail solve can become computationally infeasible in the
stiff subhorizon regime.

LSODA remains available through `SGWB_iter(engine='lsoda')`. The fast engine
does not silently switch to LSODA; fallback occurs only when the caller enables
it, and a deterministic physical guard is not retried as a numerical failure.

## Dated profiles

Older plain-grid and transition-refine comparisons are retained under
`docs/` as historical experiments. Their convergence tables describe those
recorded configurations only. See [`experiment_catalog.md`](experiment_catalog.md)
to locate each dated result and [`benchmarks.md`](benchmarks.md) for the latest
fixed-resource full-path timing profile.
