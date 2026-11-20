# CR0 high-density Fast extension: `Omega_bh2=4M–5M`

## Domain and old/current Fast screen

This CR0 slice fixes `Omega_ch2=0.1201075`, `r=.01`, `T_re=2000`, `kappa10=.01`, and `cr=0`; `n_t` remains an independent input. The box is `Omega_bh2∈[4e6,5e6] × n_t∈[-.5,0]`. Paired inputs include a 3×3 density/tilt lattice, two interior points, and an independent scrambled eight-point Sobol set (seed `20261131`).

Against pinned old Fast `872a70d`, success was `0/19`; all 19 failures were `nonfinite_background_primitive`. Current Fast at source commit `59a23ae` succeeds on `19/19`; each has finite cutoff `N_inf=47.7723413`, positive finite sampled H² (minimum `0.4532139925`), 75 frequencies, no warnings, and passes `physical_guard`. The derived tilt equals the independently supplied CR0 tilt. This is a finite sample over an adjacent rectangle, not proof for every point in it.

## Current-source z8 Reference checks

Three representative points (both 5M tilt corners and a separate Sobol holdout) were solved on the same Fast root and frequency grid. The Reference source is unchanged, `z_tail=8`, `rtol=1e-9`; the fixed-point origin uses the actual input `DN_eff` and closes from the Fast final total.

| Sample | Relative root error | Maximum same-root spectrum error | Final map closure |
|---|---:|---:|---:|
| 5M, `n_t=-.5` | `1.2257e-4` | `6.9341e-4 dex` | `1.3632e-7` |
| 5M, `n_t=0` | `1.4912e-4` | `6.7403e-4 dex` | `1.6677e-7` |
| Independent Sobol `(4.2279e6, -.0779)` | `9.6515e-5` | `7.2695e-4 dex` | `1.2746e-9` |

All meet the root `2e-4` and spectrum `1e-3 dex` gates. Reference is used only for independent validation; its implementation and thresholds were not changed.

## Complete Fast runtime A/B

The timing experiment used two isolated persistent processes, separate Numba caches, TBB, two Numba threads, CPU affinity `[0,1]`, BLAS one thread, goal grid, and `physical_guard=True`. Each case used 3 warmups and 100 alternating full solves. The pinned source baseline is `59a23ae`, which has the same production solver source as the candidate; this establishes runtime parity for the new solve points, not a speed gain over a slower old solver that cannot solve them.

| Case | P50 ratio | P95 ratio | P50 baseline/current | P95 baseline/current |
|---|---:|---:|---:|---:|
| 5M, `n_t=-.5` | `0.9928` | `1.0010` | `34.85 / 34.60 ms` | `36.97 / 37.01 ms` |
| Independent Sobol | `0.9970` | `0.9918` | `34.96 / 34.85 ms` | `37.51 / 37.20 ms` |
| Trusted CR0 control | `1.0068` | `1.0076` | `6.09 / 6.13 ms` | `6.90 / 6.95 ms` |

Roots, 75-frequency digests and spectra match bitwise in every pair. The trusted control stays within the 5% P50/P95 budget.

## Interpretation and limits

No production Fast or Reference code changed in this tranche. Success uses the already deployed unified background representation and WKB conditioning path; there are no density-specific thresholds or parameter exceptions. Combined with the separately sampled 3M–4M CR0 slice and the CR1 temperature×`kappa10` slice, this supports adjacent finite-sample improvements in more than one two-dimensional parameter slice. The result is still sampled coverage, not a mathematical continuum certificate. Raw inputs/outcomes, current-source Reference arrays, timings, figure, and hashes are captured in `cr0_density_4m_to5m_manifest.json`.
