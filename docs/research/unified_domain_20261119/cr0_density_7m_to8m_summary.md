# CR0 high-density Fast extension: `Omega_bh2=7M–8M`

## Paired Fast screen

Fix `Omega_ch2=0.1201075`, `r=.01`, `T_re=2000`, `kappa10=.01`, and `cr=0`; retain `n_t` as an independent input. The tested box `Omega_bh2∈[7e6,8e6] × n_t∈[-.5,0]` uses a 3×3 lattice, two interior points and an independent scrambled 8-point Sobol sample (seed `20261134`).

Pinned old Fast `872a70d` solves `0/19`; all failures are classified as `nonfinite_background_primitive`. Current Fast solves `19/19` with finite `N_inf=47.7723413`, positive finite H² over every generated grid (minimum `0.4532139935`), 75 frequencies, no warnings and `physical_guard` acceptance.

## Current-source z8 Reference

Three selected cases use the current Fast root and the same frequency grid. Reference remains unchanged (`z_tail=8`, `rtol=1e-9`), and fixed-point closure uses the actual input `DN_eff` origin.

| Sample | Relative root error | Maximum same-root spectrum error | Map closure residual |
|---|---:|---:|---:|
| 8M, `n_t=-.5` | `6.4442e-5` | `6.4205e-4 dex` | `9.4208e-9` |
| 8M, `n_t=0` | `3.6862e-5` | `6.1775e-4 dex` | `1.2519e-7` |
| Independent Sobol `(7.0702e6, -.2434)` | `7.5815e-5` | `6.6046e-4 dex` | `2.7791e-8` |

All points pass the root `2e-4` and spectrum `1e-3 dex` gates.

## Full-solver A/B

Two isolated persistent processes use production source baseline `59a23ae`, separate Numba caches, TBB, 2 Numba threads, CPU affinity `[0,1]`, BLAS 1, goal grid and `physical_guard=True`. The production Fast code is source-identical on both sides. Each case has 3 warmups and 500 alternating pairs.

| Case | P50 ratio | P95 ratio | Baseline/current P50 | Baseline/current P95 |
|---|---:|---:|---:|---:|
| 8M, `n_t=-.5` | `0.9983` | `1.0266` | `34.77 / 34.71 ms` | `41.25 / 42.34 ms` |
| Independent Sobol | `1.0072` | `1.0103` | `33.77 / 34.01 ms` | `36.63 / 37.01 ms` |
| Trusted CR0 | `1.0141` | `1.0418` | `5.62 / 5.70 ms` | `6.60 / 6.87 ms` |

All ratios remain within the 5% budget. The 500 paired roots, frequency digests and spectra match bitwise.

## Limits

No production Fast or Reference code changed. The samples support another adjacent CR0 density×tilt improvement through the unified numerical path, but do not prove the entire rectangle. Raw screens, Reference arrays, timing pairs, figure and source hashes are preserved in `cr0_density_7m_to8m_manifest.json`.
