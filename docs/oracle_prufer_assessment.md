# Standalone Prüfer amplitude-phase prototype

## Hypothesis

The Cartesian tensor state becomes expensive because it explicitly resolves two
rapidly oscillating components.  A Prüfer representation

`x = A sin(theta)`, `y = A cos(theta)`

may preserve the same physics while evolving one logarithmic amplitude and one
phase variable.  It is evaluated only as a standalone oracle prototype; the
formal fast kernel is unchanged.

## Protocol

- same model, initial condition, continuous background, and DOP853 solver;
- Cartesian reference versus independently written Prüfer RHS;
- default, low-T, high-T, and stiff points;
- eight representative native frequencies per point;
- `z_tail=5` and `z_tail=7`, workers=1, Numba=2, BLAS=1;
- compare handoff amplitude, phase-averaged power, phase, function evaluations,
  and wall time.

The raw artifacts are `docs/oracle_prufer_default.json`,
`docs/oracle_prufer_lowT.json`, `docs/oracle_prufer_highT.json`, and
`docs/oracle_prufer_stiff.json`.

## Results

| point | comparable modes | max amplitude error | max power error | max phase error | mean runtime ratio |
|---|---:|---:|---:|---:|---:|
| default | 14 | `9.17e-8` | `1.83e-7` | `2.42e-5 rad` | `0.537` |
| low-T | 14 | `1.19e-7` | `2.38e-7` | `1.63e-5 rad` | `0.544` |
| high-T | 14 | `9.43e-9` | `1.89e-8` | `3.17e-6 rad` | `0.545` |
| stiff | 14 | `6.05e-8` | `1.21e-7` | `2.51e-5 rad` | `0.529` |

One low-frequency native mode per point did not reach the requested tail and is
reported as unsupported, not extrapolated.  The Prüfer solver used about 65%
of the Cartesian function evaluations in the default artifact.

## Decision

**Accepted as a standalone prototype; rejected for formal-kernel promotion.**

The prototype satisfies the standalone screening signal (>30% local speedup
with sub-`1e-6` amplitude/power differences) across the four representative
points.  It has not yet compared the complete today-spectrum and `DN_gw`
observable through the formal outer loop, nor the reheating-neighborhood edge
suite.  Those are mandatory before any production integration.

## Next experiment

Run the Prüfer state through the full independent spectrum observable on the
reheating-edge and fixed Sobol points, then compare `DN_gw`, guard behavior,
determinism, and cold/warm channel cost.  Keep it reference-only until that
comparison is complete.
