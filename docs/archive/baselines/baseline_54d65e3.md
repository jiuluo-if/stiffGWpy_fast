# Current `fast_v0.2` baseline

Evidence HEAD: `54d65e321727180d0f2ee7156fea07e378f38803`  
Branch: `fast_v0.2`  
Protocol: Windows, Python 3.11.9, NumPy 2.4.4, SciPy 1.17.1, Numba 0.67.0,
Numba `workqueue`, affinity `[0, 1]`, Numba/FAST threads `2`, BLAS caps `1`,
workers `1`, concurrent processes `1`.

## Solver configuration

| item | current value |
| --- | --- |
| user mode | `fast` |
| `h` | `0.005` |
| `col_step` | `8` |
| `z_tail` | `5.0` |
| `phase_max` | `0.25` |
| frequency grid | `goal`, 76 native nodes at default point |
| reheating kink | exact kink split enabled |
| outer reuse | stability-gated first full-solve reuse enabled |
| frequency quadrature | Simpson production path; PCHIP/Gauss remain diagnostic only |

## Fresh controlled runtime and stability

The default point (`bench_fast.py --cases 0 --reps 30`) measured warm median
`9.036 ms`, p95 `9.999 ms`, minimum `8.199 ms`, with zero fallback. This is a
2-thread correctness/validation baseline; it is not comparable with older
20-thread performance artifacts.

The 24-point fast-only stability screen recorded `21 ok`, `3 explicit
shared_Neff_guard`, and `0 numerical failures`. The default point returned
`DN_gw=0.0022626969518100204`; the embedded quadrature estimate was relative
`3.4134e-4`.

## Accuracy and oracle status

The latest full-grid tail study remains the provisional scientific reference:
default `z_tail=5/6/7/8/10` gave observed systematic relative uncertainty
`3.6295e-3` and a non-monotone sequence. Stage C minimal low-T/high-T/stiff
studies gave `4.9874e-4 / 4.1532e-3 / 1.7799e-3`, with high-T still
non-monotone. Those artifacts were generated before this documentation-only
HEAD and must be refreshed before release certification; they are not treated
as a new HEAD-bound release claim.

## Dominant unresolved items

1. Oracle tail uncertainty is larger than the desired `1e-3` DN budget at the
   default/high-T cases.
2. DN quadrature/grid estimates are not yet monotone or zero-false-safe across
   the certified parameter space.
3. The controlled 2-thread warm runtime is above the `<4 ms` target; no speed
   change is accepted without a fresh profiler and accuracy A/B.
4. Full parameter-space, likelihood-bin, and 76/80/90/110-node certification
   remains incomplete.

## Decision

No production-path algorithm change is accepted from this baseline. The next
experiment must isolate per-frequency tail/adiabaticity behavior and report
`reference central ± oracle systematic` before quadrature or speed promotion.
