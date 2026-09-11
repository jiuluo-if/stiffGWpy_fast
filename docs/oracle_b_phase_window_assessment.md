# Oracle B finite phase-window prototype

## Hypothesis

The phase-averaged envelope at the existing `z=5` handoff may retain a
handoff-phase systematic.  Integrating each mode through a finite oscillatory
window and applying the phase-averaged tail only at `z=5.5`, `6`, or `7` could
provide a separately auditable Oracle-B candidate.

This is a diagnostic prototype only.  It does not change the formal fast path,
the reference default, or any production error estimate.

## Protocol

- fixed fast `DN_eff` and the native goal grid;
- DOP853, `rtol=1e-10`, `workers=1`, Numba=2, BLAS=1;
- default, low-T, high-T, and stiff points;
- four representative native frequencies per point;
- initial handoff at `z=5`, followed by windows `Δz=0`, `0.5`, `1`, and `2`;
- phase-averaged today `Opgw` computed only after the finite window.

The raw rows are stored in the four `oracle_b_phase_window_*.json` artifacts.
One low-frequency mode per point did not reach `z=5` before today and is kept as
an explicit unsupported case rather than silently extrapolated.

## Result

| point | modes reaching `z=5` | max `|B(z=5)-B(z=7)|/B(z=7)` | decision |
|---|---:|---:|---|
| default | 3 | `5.321e-3` | diagnostic only |
| low-T | 3 | `5.552e-3` | diagnostic only |
| high-T | 3 | `3.927e-3` | diagnostic only |
| stiff | 3 | `6.897e-3` | diagnostic only |

The finite-window observable changes at the `1e-3` to `1e-2` level and is
frequency- and parameter-dependent.  It therefore confirms that the existing
`z=5` phase-averaged handoff carries a measurable tail systematic, but it does
not establish that the finite-window result is an independent truth anchor.

## Decision

**Not promoted.**  The prototype is useful evidence for the oracle uncertainty
budget, but it reuses the same DOP853 tensor equations and a leading-order
analytic tail.  A promotion candidate still needs either an independently
derived adiabatic invariant or a higher-order WKB transfer, followed by a
full-grid A/B/C comparison and uncertainty accounting.

## Rollback / scope

The implementation is isolated in
`scripts/benchmark_phase_averaged_oracle.py`; deleting that diagnostic and its
focused test removes the prototype without changing the solver.
