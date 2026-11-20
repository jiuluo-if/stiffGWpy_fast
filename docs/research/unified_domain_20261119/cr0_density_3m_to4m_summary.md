# CR0 high-density Fast extension: `Omega_bh2=3M–4M`

## Parameter domain and paired Fast evidence

This CR0 slice fixes `Omega_ch2=0.1201075`, `r=.01`, `T_re=2000`, `kappa10=.01`, and `cr=0`; the tensor tilt is an independent input. The tested box is `Omega_bh2∈[3e6,4e6] × n_t∈[-.5,0]`. The paired sample contains a 3×3 lattice, two interior points, and an independently scrambled eight-point Sobol set (seed `20261130`).

Against old Fast `872a70d`, success is `0/19`; all failures are `nonfinite_background_primitive`. Current Fast `59a23ae` succeeds on `19/19`. Each candidate has finite cutoff `N_inf`, positive finite H² over its Fast grid (minimum sampled value `0.4532139925`), a complete 75-frequency output, no warnings, and passes `physical_guard`. CR0’s input `n_t` is not replaced by any CR1 relation.

## Correct-origin Reference checks

Current-source direct z8 checks (`rtol=1e-9`, `z_tail=8`) use each Fast root and its exact frequency grid. The fixed-point origin comes from `model.DN_eff_orig`, with the invariant `final DN_eff = origin + Fast root` checked to roundoff. Three maps were evaluated at each selected sample.

| Sample | Fast root relative error | Maximum same-root spectrum error | Map closure residual |
|---|---:|---:|---:|
| `Omega_bh2=4e6, n_t=-.5` | `8.2454e-5` | `7.1759e-4 dex` | `2.1580e-7` |
| `Omega_bh2=4e6, n_t=0` | `8.5923e-5` | `7.0057e-4 dex` | `9.2587e-9` |
| Independent Sobol `(3.9593e6, -.3625)` | `8.2954e-5` | `7.4337e-4 dex` | `1.2740e-8` |

The accepted 3M boundary was replayed with current Fast: roots, frequency arrays, and Fast spectra match the saved Reference input arrays bitwise. Its previous validator’s optional closure field used the post-solve total `DN_eff` as the input origin, so that field is withdrawn. Correct-origin closure on the three saved 3M maps gives maximum root error `3.4836e-5` and maximum map residual `2.5578e-7`; the same-root spectrum errors remain unchanged (maximum `7.7100e-4 dex`). Reference production code remains unchanged.

## Full Fast cost

A 100-pair alternating full-solver A/B used two isolated persistent processes and Numba caches, CPU affinity `[0,1]`, TBB, 2 Numba threads, BLAS 1 thread, goal frequency grid and `physical_guard=True`. Both sides used production source `59a23ae`; this is a throughput check of the existing Fast path, not a new solver-code performance claim.

- 4M, `n_t=-.5`: candidate/base P50/P95 `0.9910 / 0.9901`; warm candidate P50/P95 `35.1 / 37.2 ms`.
- Independent Sobol: `1.0063 / 1.0152`; warm candidate P50/P95 `34.9 / 40.6 ms`.
- Trusted CR0 control: `0.9920 / 0.9919`; root, frequency digest and spectrum are bitwise identical in every case.

The first isolated call, including JIT compile, took `23.002 / 23.023 s` for baseline/candidate; later pairs are warm. The existing 500-pair trusted-control A/B from the same production source hashes records P50/P95 `0.9963 / 1.0110`, below the 5% gate.

## Method and limits

No production equation, user parameter meaning, or dispatch rule changed in this tranche. The high-density samples use the already deployed unified background scaling and WKB-basis conditioning guard; no density-specific exception was added. The checked box is finite sampled evidence, not a proof that every point in the continuum is solvable. Source hashes, raw paired inputs/outcomes, current-source Reference arrays, corrected 3M closure, full timing samples and the figure are in `cr0_density_3m_to4m_manifest.json`.
