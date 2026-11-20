# CR0 high-density Fast region: `Omega_bh2=13M–13.75M`

## Why the accepted rectangle is narrower than the initial screen

The broader exploratory box `[13e6,14e6] × n_t[-.5,0]` produced 19/19 current-Fast successes against old Fast 0/19; all old failures were `nonfinite_background_primitive`. At its 14M low-tilt corner, however, both unchanged z8 Reference backends (Cartesian and Prüfer, `rtol=1e-9`) failed the tensor integration with `Required step size is less than spacing between numbers`. That is a Reference validation failure, not a Fast failure and not evidence that the background is mathematically invalid. Therefore 14M is excluded from the accepted accuracy rectangle.

The narrower `[13e6,13.75e6] × n_t[-.5,0]` box was separately screened on 3×3 lattice + two interior points and an independent scrambled 8-point Sobol set (seed `20261141`). Pinned old Fast `872a70d` solves `0/19`; current Fast solves `19/19`. Every current sample has finite `N_inf=47.7723413`, positive finite H² (minimum `0.4532139935`), 75 output frequencies, zero warnings and passes the physical guard. CR0 tilt remains an independent input.

## Current-source z8 Reference for the accepted box

The 13.75M tilt endpoints and an independent Sobol holdout use the Fast root and its exact frequency grid; the fixed-point origin is the true input `DN_eff`. A further 13.25M interior lattice point was checked independently as a support probe.

| Sample | Relative root error | Maximum same-root spectrum error | Map closure residual |
|---|---:|---:|---:|
| 13.75M, `n_t=-.5` | `4.3419e-5` | `5.7501e-4 dex` | `2.5683e-7` |
| 13.75M, `n_t=0` | `9.4652e-6` | `5.4140e-4 dex` | `1.4309e-7` |
| Sobol `(13.0267e6, -.3566)` | `2.4868e-5` | `6.7529e-4 dex` | `9.3832e-8` |
| Interior 13.25M, `n_t=-.375` | `2.9473e-5` | `5.7523e-4 dex` | `1.0973e-7` |

All meet the `2e-4` root and `1e-3 dex` spectrum requirements. The 14M boundary remains unvalidated by Reference and is explicitly excluded.

## Full-solver runtime A/B

The performance gate uses source baseline `59a23ae` and source-identical current production code, isolated persistent processes, separate Numba caches, TBB, 2 threads, affinity `[0,1]`, BLAS 1, goal grid, `physical_guard=True`, 3 warmups and 500 alternating complete solves.

| Case | P50 ratio | P95 ratio | Baseline/current P50 | Baseline/current P95 |
|---|---:|---:|---:|---:|
| 13.75M, `n_t=-.5` | `1.0035` | `0.9995` | `34.09 / 34.21 ms` | `36.23 / 36.21 ms` |
| Independent Sobol | `1.0015` | `1.0013` | `34.14 / 34.19 ms` | `35.56 / 35.61 ms` |
| Trusted CR0 | `1.0075` | `1.0191` | `5.42 / 5.46 ms` | `6.19 / 6.31 ms` |

All remain below the 5% runtime budget. Paired roots, frequency digests and spectra match bitwise.

## Limits and interpretation

No production Fast or Reference code changed. The accepted finite-sample expansion is `[13M,13.75M]` only; no continuum certification is claimed. The 14M fast-only screen is retained as a diagnostic, but excluded from accuracy success because both independent Reference backends fail at that corner. Scripts, raw inputs/outcomes, the Reference failure diagnostic, source hashes and figure are captured in `cr0_density_13m_to13_75m_manifest.json`.
