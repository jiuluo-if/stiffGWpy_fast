# CR0 high-density Fast extension: `Omega_bh2=10M–11M`

## Paired old/current Fast screen

Fix `Omega_ch2=0.1201075`, `r=.01`, `T_re=2000`, `kappa10=.01`, and `cr=0`; `n_t` remains an independent input. The sampled rectangle `Omega_bh2∈[10e6,11e6] × n_t∈[-.5,0]` contains a 3×3 lattice, two interior points and a scrambled independent eight-point Sobol set (seed `20261137`).

Pinned old Fast `872a70d` succeeds on `0/19`, with all failures classified as `nonfinite_background_primitive`. Current Fast succeeds on `19/19`; each has finite `N_inf=47.7723413`, positive finite sampled H² (minimum `0.4532139935`), full 75-frequency output, no warnings and accepted physical guard.

## Current-source z8 Reference

The three selected points use the Fast root and exactly its frequency grid. Reference remains unchanged (`z_tail=8`, `rtol=1e-9`); fixed-point closure uses the true input `DN_eff` origin.

| Sample | Relative root error | Maximum same-root spectrum error | Map closure residual |
|---|---:|---:|---:|
| 11M, `n_t=-.5` | `4.5421e-5` | `6.0209e-4 dex` | `9.6058e-8` |
| 11M, `n_t=0` | `2.2955e-5` | `5.7270e-4 dex` | `2.6102e-8` |
| Independent Sobol `(10.7715e6, -.4083)` | `4.2975e-5` | `6.8906e-4 dex` | `3.1447e-8` |

All satisfy the `2e-4` root and `1e-3 dex` spectrum limits.

## Full-solver runtime A/B

The benchmark uses isolated persistent processes with source baseline `59a23ae` and source-identical current Fast, separate Numba caches, TBB, two threads, CPU affinity `[0,1]`, BLAS 1, goal grid, physical guard, 3 warmups and 500 alternating complete solves.

| Case | P50 ratio | P95 ratio | Baseline/current P50 | Baseline/current P95 |
|---|---:|---:|---:|---:|
| 11M, `n_t=-.5` | `1.0000` | `1.0067` | `32.37 / 32.37 ms` | `37.03 / 37.28 ms` |
| Independent Sobol | `0.9978` | `1.0022` | `32.08 / 32.01 ms` | `33.59 / 33.67 ms` |
| Trusted CR0 | `1.0027` | `0.9968` | `5.48 / 5.49 ms` | `6.16 / 6.14 ms` |

All ratios remain within the 5% allowance, and roots, frequency digests and spectra match bitwise in every pair.

## Limits

No Fast or Reference production code changed. The result supports sampled coverage for an adjacent CR0 density×tilt box, not every point in its continuum. Raw inputs/outcomes, Reference arrays, timings, source hashes, scripts and figure are preserved by `cr0_density_10m_to11m_manifest.json`.
