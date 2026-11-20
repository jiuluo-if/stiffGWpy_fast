# CR0 high-density Fast extension: `Omega_bh2=9M–10M`

## Paired Fast screen

Fix `Omega_ch2=0.1201075`, `r=.01`, `T_re=2000`, `kappa10=.01`, and `cr=0`; the tensor tilt remains an independent input. The sampled rectangle `Omega_bh2∈[9e6,10e6] × n_t∈[-.5,0]` uses a 3×3 lattice, two interior points and an independent scrambled eight-point Sobol set (seed `20261136`).

Pinned old Fast `872a70d` succeeds on `0/19`; every old failure is `nonfinite_background_primitive`. Current Fast succeeds on `19/19`, with finite cutoff `N_inf=47.7723413`, positive finite sampled H² (minimum `0.4532139935`), full 75-frequency outputs, no warnings and accepted physical guard.

## Current-source z8 Reference checks

Three points use the current Fast root and identical frequency grid. Reference remains unchanged (`z_tail=8`, `rtol=1e-9`) and fixed-point closure starts from the true input `DN_eff`.

| Sample | Relative root error | Maximum same-root spectrum error | Map closure residual |
|---|---:|---:|---:|
| 10M, `n_t=-.5` | `5.1285e-5` | `6.2001e-4 dex` | `3.5154e-8` |
| 10M, `n_t=0` | `3.0930e-5` | `5.9117e-4 dex` | `4.0239e-8` |
| Independent Sobol `(9.7344e6, -.2248)` | `4.0288e-5` | `6.0586e-4 dex` | `3.5936e-8` |

All meet the `2e-4` root and `1e-3 dex` spectrum criteria.

## Complete Fast runtime A/B

The A/B uses isolated persistent processes, source baseline `59a23ae` and source-identical current Fast, separate Numba caches, TBB, CPU `[0,1]`, 2 Numba threads, BLAS 1, goal frequency grid, `physical_guard=True`, 3 warmups and 500 alternating pairs.

| Case | P50 ratio | P95 ratio | Baseline/current P50 | Baseline/current P95 |
|---|---:|---:|---:|---:|
| 10M, `n_t=-.5` | `0.9977` | `0.9896` | `32.32 / 32.25 ms` | `34.76 / 34.40 ms` |
| Independent Sobol | `0.9973` | `1.0001` | `32.46 / 32.37 ms` | `37.18 / 37.19 ms` |
| Trusted CR0 | `1.0020` | `1.0047` | `5.52 / 5.53 ms` | `6.15 / 6.18 ms` |

All ratios are within 5%; every paired root, frequency digest and spectrum is bitwise identical.

## Limits

No production Fast or Reference implementation changed. The result is finite-sample evidence for one more adjacent CR0 density×tilt box, not a continuum proof. Source hashes, raw screen and Reference data, timing pairs, reproduction scripts and figure are in `cr0_density_9m_to10m_manifest.json`.
