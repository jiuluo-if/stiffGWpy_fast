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

The fixed-frequency full-output replay also compared today `Ogw`, `Oj`,
`Opgw`, and integrated `DN_gw` on eight native frequencies. Across the four
points and both handoff depths, the largest `DN_gw` relative difference was
`2.26e-9`; the largest spectrum-component difference was `1.67e-7`. Repeating
the default-point run produced bitwise-identical numerical fields for the
recorded rows and spectrum summary; timing fields were excluded.

The same eight-frequency subset was then run through the reference outer
self-consistency algorithm, including its `DN_eff` guard and convergence test.
All four points converged in one or two iterations at both handoff depths, with
the Prüfer and Cartesian iteration counts identical. The largest
outer-loop `DN_gw` relative difference was `3.44e-9`, and no point crossed the
`DN_eff` guard.

The boundary follow-up used two reheating-axis points (`T_re=12.6` and
`T_re=7.94e5`) plus fixed Sobol points `sobol_000`, `sobol_002`, and
`sobol_006`, again with eight native frequencies and both handoff depths. The
low reheating edge and all three Sobol points converged with matching outer
iteration counts; their largest outer `DN_gw` relative differences were
`8.33e-10`, `2.12e-9`, `2.06e-9`, and `2.50e-9`, respectively. The high
reheating edge produced an explicit `shared_Neff_guard` on both depths in both
implementations, so it is classified as a physical rejection rather than a
numerical failure. Across the four accepted points, the largest handoff power
error was `2.37e-7` and the largest observed Prüfer/Cartesian runtime ratio was
`0.913`.

The complete ten-point parameter-axis edge suite was then replayed, together
with fixed Sobol points `sobol_000`, `sobol_002`, `sobol_006`, `sobol_010`, and
`sobol_015`. This produced 26 accepted outer comparisons and four physical
guard records (`edge_tre_hi` and `edge_nt_blue`, both handoff depths), with no
numerical failures. On accepted comparisons, outer iteration counts matched
in every case; the maximum `DN_gw` relative difference was `6.44e-9`, the
maximum handoff power error was `3.26e-7`, and the maximum runtime ratio was
`0.929`. These results expand the standalone screening evidence but do not
certify the full production frequency grid.

## Decision

**Accepted as a standalone prototype; rejected for formal-kernel promotion.**

The prototype satisfies the standalone screening signal on the accepted
representative, boundary, and Sobol points, with sub-`1e-6` amplitude/power
differences and matching physical guard behavior. It remains a reference-only
prototype: the sampled frequencies are still sparse, the full production
frequency grid has not been certified, and no formal-kernel integration is
proposed.

## Next experiment

Add a full-grid `DN_gw` comparison and a deterministic replay for the accepted
edge/Sobol cases. Keep it reference-only until the full-grid comparison and
its oracle uncertainty are documented.
