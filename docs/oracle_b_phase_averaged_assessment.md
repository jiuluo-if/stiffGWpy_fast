# Oracle B assessment: phase-averaged energy observable

Source diagnostic: `e95185e1c50752f62cb89cbda9f7a14532855443`  
Input artifact: [`tail_diagnostics_f939.json`](tail_diagnostics_f939.json)

## Hypothesis

Replacing the oscillatory hand-off state by a phase-averaged energy observable
could reduce the frozen-tail sensitivity and provide an independent Oracle B.

## Algebraic check

The current DOP853 reference already computes, at the hand-off,

```text
C = sqrt((x_f^2 + y_f^2) / 2)
T_hf(N) = C exp(-z_f + N_event - N)
```

and then constructs the post-handoff energy-density terms from that envelope.
The phase `atan2(x_f, y_f)` is not used in the current tail propagation. Thus a
phase-averaged energy observable evaluated from the same hand-off state is
algebraically the same observable as the existing frozen-tail branch. It is not
an independent oracle and cannot be used to reduce the observed systematic
bound.

## Evidence

The diagnostic solved four named points, eight native representative
frequencies per point, and `z_tail=5/6/7/8/10` with workers=1 and Numba=2:

| point | max relative `Omega_GW(z_tail)` change vs z=10 | max `|omega'/omega^2|` at z=5 | max at z=10 |
| --- | ---: | ---: | ---: |
| default | `6.602e-3` | `1.348e-2` | `9.080e-5` |
| low-T | `5.854e-3` | `1.347e-2` | `8.188e-5` |
| high-T | `6.719e-3` | `1.348e-2` | `9.080e-5` |
| stiff | `5.891e-3` | `1.348e-2` | `9.080e-5` |

The per-frequency phase, amplitude, `omega`, first/second adiabaticity
indicators, and DN weights are retained in the JSON artifact. The evidence
supports a hand-off/truncation systematic, but does not distinguish a new
phase-averaging rule from the existing envelope rule.

## Decision

**Rejected for promotion.** No production code or uncertainty reduction is
claimed. Oracle B remains unimplemented as an independent observable.

The next valid Oracle B prototype must average over a finite post-handoff
phase window or integrate an independently derived adiabatic invariant, then
compare against Oracle A and a deeper-tail reference. A relabeling of the
current envelope is insufficient.
