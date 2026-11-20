# CR0 high-density Fast extension: `Omega_bh2=8M–9M`

## Paired Fast screen

Fix `Omega_ch2=0.1201075`, `r=.01`, `T_re=2000`, `kappa10=.01`, and `cr=0`; CR0 `n_t` is an independent input. The sampled box `Omega_bh2∈[8e6,9e6] × n_t∈[-.5,0]` contains a 3×3 lattice, two interior points and an independent scrambled 8-point Sobol sample (seed `20261135`).

Pinned old Fast `872a70d` succeeds on `0/19`; all old failures are `nonfinite_background_primitive`. Current Fast solves `19/19`, with finite `N_inf=47.7723413`, positive finite sampled H² (minimum `0.4532139935`), full 75-frequency output, no warnings and physical guard acceptance.

## Current-source z8 Reference

Three cases use the current Fast root and exact Fast grid, with unchanged Reference (`z_tail=8`, `rtol=1e-9`) and the correct input `DN_eff` fixed-point origin.

| Sample | Relative root error | Maximum same-root spectrum error | Map closure residual |
|---|---:|---:|---:|
| 9M, `n_t=-.5` | `5.3721e-5` | `6.2684e-4 dex` | `2.6491e-8` |
| 9M, `n_t=0` | `2.5019e-5` | `6.0106e-4 dex` | `4.8786e-8` |
| Independent Sobol `(8.9459e6, -.0542)` | `2.9187e-5` | `6.1662e-4 dex` | `1.0789e-7` |

Every selected case meets the root `2e-4` and spectrum `1e-3 dex` gates.

## Full Fast runtime A/B

Isolated persistent processes compare source baseline `59a23ae` with source-identical current Fast, using separate Numba caches, TBB, CPU affinity `[0,1]`, 2 Numba threads, BLAS 1, goal grid and `physical_guard=True`. Each case includes 3 warmups and 500 alternating complete solves.

| Case | P50 ratio | P95 ratio | Baseline/current P50 | Baseline/current P95 |
|---|---:|---:|---:|---:|
| 9M, `n_t=-.5` | `0.9900` | `1.0167` | `36.81 / 36.44 ms` | `41.73 / 42.42 ms` |
| Independent Sobol | `0.9934` | `0.9953` | `34.45 / 34.23 ms` | `40.11 / 39.93 ms` |
| Trusted CR0 | `0.9941` | `0.9974` | `5.97 / 5.93 ms` | `7.37 / 7.35 ms` |

All P50/P95 ratios remain within the 5% gate. All 500 paired roots, frequency digests and spectra match bitwise.

## Limits

No production Fast or Reference code changed in this tranche. The result extends sampled evidence for the unified CR0 density×tilt mechanism through an adjacent box; it does not certify every point in the continuum. Raw screens, Reference arrays, timings, source hashes and the figure are preserved in `cr0_density_8m_to9m_manifest.json`.
