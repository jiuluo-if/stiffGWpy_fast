# CR0 high-density Fast region: `Omega_bh2=13M–13.984375M`

## Fast coverage and the 14M Reference boundary

The full exploratory box `[13M,14M] × n_t[-.5,0]` has 19 same-input lattice/interior/Sobol samples: pinned old Fast `872a70d` succeeds `0/19`; current Fast succeeds `19/19`, all old failures `nonfinite_background_primitive`. At the exact 14M low-tilt corner, Fast succeeds, but both unchanged z8 Reference backends (Cartesian and Prüfer, `rtol=1e-9`) fail tensor integration with `Required step size is less than spacing between numbers`. The failure is a Reference numerical integration limit; it does not establish physical invalidity or a Fast failure, and it cannot count as a precision pass. At 14M high tilt, the unchanged Cartesian Reference passes.

Boundary probes at 13.875M, 13.9375M, 13.96875M and 13.984375M for low tilt pass z8 Reference. The 13.984375M high-tilt endpoint also passes. We therefore define the accepted rectangle as `[13e6,13.984375e6] × n_t[-.5,0]`, without including the failing 14M low-tilt point.

A dedicated screen of this accepted rectangle uses a 3×3 lattice, two interior points and an independent scrambled eight-point Sobol set (seed `20261142`). Old/current Fast success is `0/19→19/19`; all current points have finite `N_inf=47.7723413`, positive finite sampled H² (minimum `0.4532139935`), 75 frequencies, no warnings and pass the physical guard. CR0 tilt remains an independent input.

## Current-source z8 Reference

The exact upper-edge tilt endpoints and an independent Sobol holdout use the same Fast root and frequency grid, true input `DN_eff` origin, unchanged Reference, `z_tail=8`, `rtol=1e-9`. A 13.25M interior sample is also retained as an independent support check.

| Sample | Relative root error | Maximum same-root spectrum error | Map closure residual |
|---|---:|---:|---:|
| 13.984375M, `n_t=-.5` | `4.5617e-5` | `6.4984e-4 dex` | `1.8183e-7` |
| 13.984375M, `n_t=0` | `1.4377e-5` | `6.1562e-4 dex` | `1.1602e-7` |
| Independent Sobol `(13.0154M, -.4371)` | `2.9877e-5` | `6.7808e-4 dex` | `9.0697e-8` |
| 13.25M, `n_t=-.375` interior | `2.9473e-5` | `5.7523e-4 dex` | `1.0973e-7` |

All meet the `2e-4` root and `1e-3 dex` spectrum limits. The exact 14M low-tilt face remains excluded from the accuracy claim.

## Complete Fast runtime A/B

The A/B uses source baseline `59a23ae` and source-identical current Fast in isolated persistent processes, separate Numba caches, TBB, 2 Numba threads, CPU affinity `[0,1]`, BLAS 1, goal grid, physical guard, 3 warmups and 500 alternating complete solves.

| Case | P50 ratio | P95 ratio | Baseline/current P50 | Baseline/current P95 |
|---|---:|---:|---:|---:|
| 13.984375M, `n_t=-.5` | `0.9992` | `1.0051` | `34.14 / 34.11 ms` | `35.93 / 36.11 ms` |
| Independent Sobol | `0.9971` | `1.0080` | `34.23 / 34.14 ms` | `35.76 / 36.05 ms` |
| Trusted CR0 | `0.9957` | `1.0078` | `5.44 / 5.41 ms` | `6.11 / 6.16 ms` |

All P50/P95 ratios remain below 1.05. Roots, frequency digests and spectra match bitwise in every pair.

## Limits

No production Fast or Reference code changed. This is finite-sample support for `[13M,13.984375M]`; it does not mathematically certify the continuum. The 14M low-tilt Reference failure, both backend messages, full-box exploratory screen, narrowed-box screen, Reference arrays, timing samples, hashes, scripts and figure are included in `cr0_density_13m_to13_984375m_manifest.json`.
