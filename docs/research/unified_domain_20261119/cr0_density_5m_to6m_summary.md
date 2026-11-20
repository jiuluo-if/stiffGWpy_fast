# CR0 high-density Fast extension: `Omega_bh2=5M–6M`

## Domain and old/current Fast screen

This CR0 slice fixes `Omega_ch2=0.1201075`, `r=.01`, `T_re=2000`, `kappa10=.01`, and `cr=0`; `n_t` is an independent input. The box `Omega_bh2∈[5e6,6e6] × n_t∈[-.5,0]` was sampled with a 3×3 lattice, two interior points and a scrambled independent eight-point Sobol set (seed `20261132`).

On identical inputs, pinned old Fast `872a70d` succeeds on `0/19`; all failures are `nonfinite_background_primitive`. Current Fast succeeds on `19/19`, each with finite `N_inf=47.7723413`, positive finite sampled H² (minimum `0.4532139935`), a full 75-frequency spectrum, no warnings and accepted physical guard. `cr=0` keeps the supplied tilt independent.

## Current-source z8 Reference checks

Three representative cases use the same Fast fixed-point root and its exact frequency grid. Reference remains unchanged (`z_tail=8`, `rtol=1e-9`), and its fixed-point closure uses the true input origin `DN_eff`.

| Sample | Relative root error | Maximum same-root spectrum error | Map closure residual |
|---|---:|---:|---:|
| 6M, `n_t=-.5` | `1.2161e-4` | `6.7416e-4 dex` | `3.4958e-8` |
| 6M, `n_t=0` | `1.4323e-4` | `6.5340e-4 dex` | `6.6083e-8` |
| Independent Sobol `(5.8403e6, -.2778)` | `1.3738e-4` | `7.3053e-4 dex` | `4.7165e-8` |

All three pass the `2e-4` root and `1e-3 dex` spectrum gates.

## Full Fast performance A/B

The benchmark compared isolated persistent processes using source baseline `59a23ae`, separate Numba caches, TBB, CPU affinity `[0,1]`, two Numba threads, BLAS one thread, goal frequency grid, and `physical_guard=True`. Both production source trees are identical; the comparison checks cost parity at newly solvable inputs. The acceptance table uses 500 alternating pairs and 3 warmups.

| Case | P50 ratio | P95 ratio | Baseline/current P50 | Baseline/current P95 |
|---|---:|---:|---:|---:|
| 6M, `n_t=-.5` | `0.9977` | `1.0028` | `34.13 / 34.06 ms` | `35.61 / 35.70 ms` |
| Independent Sobol | `0.9969` | `1.0009` | `35.18 / 35.07 ms` | `38.35 / 38.39 ms` |
| Trusted CR0 | `1.0003` | `1.0070` | `5.76 / 5.76 ms` | `6.47 / 6.51 ms` |

All 500 pairs have bitwise identical roots, frequency digests and spectra. A 100-pair pilot showed a noisy trusted-control P95 ratio of `1.0525`; the predeclared 500-pair recheck brought the side-wise P95 ratio to `1.0070`, so the resolved measurement stays within the 5% gate. The 100-pair raw data is retained in the manifest.

## Interpretation and limits

No Fast production code, Reference code, equation, dispatch threshold, or parameter meaning changed. The sampled area extends the same unified numerical mechanism across an additional adjacent CR0 density×tilt slice. Together with the validated 3M–4M and 4M–5M CR0 slices and the distinct CR1 temperature×`kappa10` slice, results support improvements in multiple two-dimensional sampled regions. They do not mathematically certify every point in the rectangles. Raw paired results, the pilot and 500-pair timings, Reference arrays, figure, and source hashes are captured in `cr0_density_5m_to6m_manifest.json`.
