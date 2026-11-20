# CR0 high-density Fast extension: `Omega_bh2=6M–7M`

## Paired Fast screen

Hold `Omega_ch2=0.1201075`, `r=.01`, `T_re=2000`, `kappa10=.01`, and `cr=0` fixed; `n_t` remains an independent input. The sampled rectangle is `Omega_bh2∈[6e6,7e6] × n_t∈[-.5,0]`, with a 3×3 lattice, two interior points and an independently scrambled 8-point Sobol set (seed `20261133`).

Against pinned old Fast `872a70d`, success is `0/19`; all old failures are `nonfinite_background_primitive`. Current Fast succeeds on all `19/19`, with finite cutoff `N_inf=47.7723413`, positive finite H² over the Fast grid (minimum `0.4532139935`), 75 frequencies, no warnings and physical guard acceptance. The supplied CR0 tilt is preserved.

## Current-source Reference checks

Three current-source z8 Reference cases use the Fast root and identical Fast frequency grid, `rtol=1e-9`, with fixed-point origin taken from the input `DN_eff`:

| Sample | Relative root error | Maximum same-root spectrum error | Map closure residual |
|---|---:|---:|---:|
| 7M, `n_t=-.5` | `8.5756e-5` | `6.9349e-4 dex` | `6.2642e-8` |
| 7M, `n_t=0` | `6.4879e-5` | `6.7064e-4 dex` | `2.2945e-8` |
| Independent Sobol `(6.5986e6, -.0751)` | `9.7423e-5` | `6.5242e-4 dex` | `1.2114e-7` |

All pass the `2e-4` root and `1e-3 dex` spectrum criteria. Reference code is unchanged.

## Full-solver performance

A/B compares isolated persistent Fast processes from source baseline `59a23ae` and current source-identical code, separate Numba caches, TBB, two threads, CPU affinity `[0,1]`, BLAS one thread, goal grid, and `physical_guard=True`. Each case uses 3 warmups and 500 alternating pairs.

| Case | P50 ratio | P95 ratio | Baseline/current P50 | Baseline/current P95 |
|---|---:|---:|---:|---:|
| 7M, `n_t=-.5` | `0.9954` | `1.0005` | `35.59 / 35.42 ms` | `37.76 / 37.78 ms` |
| Independent Sobol | `0.9954` | `0.9953` | `35.56 / 35.39 ms` | `37.62 / 37.44 ms` |
| Trusted CR0 | `0.9957` | `0.9849` | `5.96 / 5.93 ms` | `6.82 / 6.72 ms` |

Every pair has identical roots, frequency hashes and spectra. All P50/P95 ratios stay within 5%.

## Limits

No Fast or Reference production code changed. This is additional finite-sample evidence for the unified background path across another adjacent CR0 density×tilt slice; it does not prove every point of the box succeeds. Earlier 3M–4M, 4M–5M and 5M–6M slices and the distinct CR1 temperature×`kappa10` slice provide the multi-region context. Raw data, source hashes, figure and reproduction scripts are included in `cr0_density_6m_to7m_manifest.json`.
