# CR0 high-density Fast extension: `Omega_bh2=12M–13M`

## Paired Fast screen

Fix `Omega_ch2=0.1201075`, `r=.01`, `T_re=2000`, `kappa10=.01`, and `cr=0`; retain `n_t` as an independent input. The sampled box `Omega_bh2∈[12e6,13e6] × n_t∈[-.5,0]` includes a 3×3 lattice, two interior points and an independent scrambled Sobol sample of eight points (seed `20261139`).

Pinned old Fast `872a70d` succeeds on `0/19`, with every failure `nonfinite_background_primitive`. Current Fast succeeds on `19/19`; sampled backgrounds have finite `N_inf=47.7723413`, positive finite H² (minimum `0.4532139935`), full 75-frequency outputs, zero warnings and accepted `physical_guard`.

## Current-source z8 Reference

Three cases compare the Fast fixed-point root and its exact frequency grid with unchanged Reference (`z_tail=8`, `rtol=1e-9`), iterating from the actual input `DN_eff` origin.

| Sample | Relative root error | Maximum same-root spectrum error | Map closure residual |
|---|---:|---:|---:|
| 13M, `n_t=-.5` | `3.3581e-5` | `6.7689e-4 dex` | `8.0921e-9` |
| 13M, `n_t=0` | `3.8469e-6` | `6.4345e-4 dex` | `8.6884e-8` |
| Independent Sobol `(12.8247e6, -.1076)` | `3.0760e-6` | `5.6246e-4 dex` | `5.4801e-8` |

All samples satisfy the root `2e-4` and spectrum `1e-3 dex` limits.

## Full Fast performance A/B

Two isolated persistent processes compare source baseline `59a23ae` with current source-identical Fast using separate Numba caches, TBB, 2 threads, CPU affinity `[0,1]`, BLAS 1, goal grid and `physical_guard=True`. Each case has three warmups and 500 alternating full solves.

| Case | P50 ratio | P95 ratio | Baseline/current P50 | Baseline/current P95 |
|---|---:|---:|---:|---:|
| 13M, `n_t=-.5` | `1.0029` | `1.0170` | `34.20 / 34.30 ms` | `35.66 / 36.26 ms` |
| Independent Sobol | `1.0030` | `0.9905` | `34.16 / 34.27 ms` | `36.06 / 35.72 ms` |
| Trusted CR0 | `1.0016` | `0.9956` | `5.45 / 5.46 ms` | `6.16 / 6.14 ms` |

All ratios are within 5%; the 500 paired roots, frequency digests and spectra are bitwise identical.

## Limits

This validation tranche changes no production Fast or Reference code. It supplies finite-sample evidence for another adjacent CR0 density×tilt region; it does not certify the whole rectangle. Raw screen/Reference/timing data, source hashes, figure and runners are captured in `cr0_density_12m_to13m_manifest.json`.
