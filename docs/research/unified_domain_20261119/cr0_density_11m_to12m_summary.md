# CR0 high-density Fast extension: `Omega_bh2=11M–12M`

## Paired old/current Fast screen

Fix `Omega_ch2=0.1201075`, `r=.01`, `T_re=2000`, `kappa10=.01`, and `cr=0`. The input tensor tilt remains independent. The sampled box `Omega_bh2∈[11e6,12e6] × n_t∈[-.5,0]` uses a 3×3 lattice, two interior points and an independent scrambled 8-point Sobol design (seed `20261138`).

Pinned old Fast `872a70d` succeeds on `0/19`; all old failures remain `nonfinite_background_primitive`. Current Fast solves `19/19` with finite cutoff `N_inf=47.7723413`, positive finite H² (minimum sampled `0.4532139935`), 75 frequencies, no warnings and accepted `physical_guard`.

## Current-source z8 Reference checks

Three selected samples use the Fast root and identical frequency grid. Reference is unchanged (`z_tail=8`, `rtol=1e-9`) and the fixed-point mapping starts from the input `DN_eff` origin.

| Sample | Relative root error | Maximum same-root spectrum error | Map closure residual |
|---|---:|---:|---:|
| 12M, `n_t=-.5` | `3.7307e-5` | `5.9156e-4 dex` | `1.1203e-7` |
| 12M, `n_t=0` | `3.4971e-6` | `5.6071e-4 dex` | `2.4526e-8` |
| Independent Sobol `(11.1456e6, -.4472)` | `4.2548e-5` | `6.1427e-4 dex` | `6.0530e-8` |

All satisfy the root `2e-4` and spectrum `1e-3 dex` thresholds.

## Complete-solve runtime A/B

Use isolated persistent Fast processes from source baseline `59a23ae` and source-identical current code, separate Numba caches, TBB, two Numba threads, CPU affinity `[0,1]`, BLAS 1, goal grid, physical guard, three warmups and 500 alternating full solves.

| Case | P50 ratio | P95 ratio | Baseline/current P50 | Baseline/current P95 |
|---|---:|---:|---:|---:|
| 12M, `n_t=-.5` | `0.9977` | `0.9994` | `33.90 / 33.83 ms` | `38.00 / 37.98 ms` |
| Independent Sobol | `0.9973` | `0.9934` | `31.88 / 31.79 ms` | `33.57 / 33.35 ms` |
| Trusted CR0 | `1.0051` | `1.0221` | `5.41 / 5.43 ms` | `6.00 / 6.13 ms` |

All ratios remain below 1.05. Paired roots, frequency digests and spectra match bitwise.

## Limits

No production Fast/Reference source changed. This is finite-sample evidence for another adjacent CR0 density×tilt slice; it does not certify every point in the box. Raw screens, Reference arrays, timing pairs, source hashes, scripts and the figure are in `cr0_density_11m_to12m_manifest.json`.
