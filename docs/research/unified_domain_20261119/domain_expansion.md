# 2026-11-19 unified Fast evaluation-domain candidate

## Mathematical change

For `N < N_re_abs = N_inf - N_re`, the existing physical background is the analytic matter-dominated extension. The exact horizon coordinate in this branch is

`f_hor(N) = f_hor(N_re_abs) - (N - N_re_abs) / (2 ln 10)`.

The prior `grid_independent_freqs.f_hor_cont` first called `H2_vec(N)` and only then replaced the result with this analytic expression. That `H2` was unused. At `N_inf=740`, evaluating `exp(N_inf-N)` at `N=0` is outside float64 range, even though the requested matter-era horizon value and all subsequent physical background quantities are finite. The Fast-side frequency builder now applies the matter-era expression before evaluating `H2`. The exact transition-grid builders likewise evaluate `H2` only on and after the physical reheating breakpoint, then extend the earlier branch analytically. No physical transition is moved or smoothed.

The propagation-scale representation uses `P(N)=Phi(N)-N`, so `Psi=2P`, `S2=exp(2P)`, and `S2inv=exp(-P)`. For each mode starting at `j0`, the candidate stores the transfer state multiplied by the constant `exp(P(j0))`. The linear transfer equation preserves this constant factor. When assembling an observable at node `k`, the state is multiplied by `exp(P(k)-P(j0))`, giving `exp(P(k))*state_original`, exactly the factor previously supplied by `sqrt(S2(k))`. When the stored state is renormalized by a factor `q`, `log(q)` is accumulated and restored at assembly. The WKB handoff uses the same reconstructed state and the identity `exp(N_k-N_slot-z_k) * exp((f_slot-f_k) ln 10)` for its terminal tensor amplitude.

## Reproduction and region evidence

Baseline commit: `3eb9589952833301b197b650a90d0709a5b229cc` (`fast_v0.2`). The 3x3 `cr=0` rectangle spans `T_re={10, 2000, 100000} GeV` and `kappa10={0.001, 0.01, 0.1}`. For every point, `DN_re` is derived as `740-N_re`, giving `N_inf=740`; `r=0.01`, `n_t=1`, and `A_s` places the physical blue-tilt cutoff at the same pivot-adjacent location.

- Baseline: 9/9 points raise `OverflowError: math range error` while constructing the frequency support.
- Candidate: 9/9 points return finite spectra. Root values are `4.81640928069` to the shown precision. The transition location moves from `log10(f_re/Hz)=-5.02` to `4.05` across the rectangle.
- Independent 27-node Reference maps at the two rectangle corners differ from the Fast root by `8.36e-7` relative.
- The `cr=0`, `DN_re x n_t` cutoff-admissible slice has 12/12 baseline exceptions and 12/12 candidate successes. It is sampled finite-region evidence, not a continuous-domain proof.
- A seeded 32-point 3D Sobol training set over `T_re`, `kappa10`, and `n_t` has 32/32 baseline `OverflowError` failures and 32/32 candidate successes. A separate seed gives a 16-point holdout with the same 16/16 versus 16/16 result. Direct 27-node Reference checks on every holdout point give maximum root relative error `8.365e-7` and maximum spectral error `3.834e-7 dex`; full frequencies and spectra are archived per point.
- The cutoff-to-reheating-feature distance was screened at 2, 4, 6, 8, and 10 log-frequency units. Reference checks at margins 2, 4, 6, 8, and 10 fail either the spectrum or the self-consistent-root gate. At margin 10, the 27-node Reference fixed-point root closes at `2.7611924`, while the 14-node Fast result is `2.8224694`. The main error source is under-resolved frequency integration around the cutoff. Raising adaptive frequency resolution recovers the Reference root but requires hundreds to over a thousand modes and costs about 0.9–2.0 s per call in this environment, so it is not an accepted Fast path.
- A 3x2 wide-spectrum slice that brings the transition into the sampled support fails the accuracy gate. At two points the 127-node Reference maps differ from the default-grid Fast roots by 31% and 60%. An independently closed 395-node Reference root also differs from Fast by `9.03e-4` on another transition-in-support point. These results are retained as rejected regimes.
- The `cr>0` 3x3 control grid is valid under the derived consistency tilt and cutoff constraints, but baseline also succeeds there. It is a regression control, not an expansion result.

Artifacts in this directory include the exact baseline and candidate grid rows, Sobol train and independent holdout rows, frequency-refinement records, Reference root closure, per-mode error table, timing screen, and baseline/candidate source hashes. `domain_coverage_20261119.png` visualizes the sampled grids.

## Performance and verification

An interleaved trusted-point benchmark used eight alternating fresh processes, each with one warmup and 12 timed calls (96 samples per version). Baseline/candidate P50/P95 were 10.34/12.33 ms and 10.16/11.58 ms. The candidate shows no measured regression on this workload. The separate cold-to-warm process startup medians were 262/258 ms. New-region calls on the 3x3 narrow-cutoff slice took 46–47 ms warm after compilation. The broad feature-support cases are rejected on accuracy.

The added regression test checks the 3x3 rectangle. The focused file passes 11 tests. The current full suite passes 422 tests with 9 deselected. Ruff, Mypy, source distribution and wheel builds also pass. GitHub Actions run `44462717007` for commit `b6140947c909fe2891d562b2e3cef05acf409b4b` completed successfully across canonical regression, static checks, package build, Cobaya smoke, and Python 3.9–3.13 compatibility.

## 2026-11-19 candidate screen: exact matter-era Bessel transfer (rejected)

For the Fast tensor equation `T'' + (3 - 1.5 sigma) T' + exp(2z) T = 0`, a constant matter-era branch has `sigma=1` and `z'=1/2`. With `x=2 exp(z)`, it reduces to `x^2 T_xx + 4x T_x + x^2 T = 0`; its independent solutions are proportional to `(sin x - x cos x)/x^3` and `(cos x + x sin x)/x^3`. This is the exact order-3/2 Bessel solution. The state conversion follows from the Fast matrix system and the existing `S2` gauge; the matching point remains the actual unsmoothed reheating event.

The source check read Turner, White & Lidsey, “Tensor Perturbations in Inflationary Models as a Probe of Cosmology,” §2.1 and Eq. (13), which states the matter-era `3 j1(k tau)/(k tau)` solution, and Greene et al., “Observable CMB Tensor Modes from Cosmological Phase Transitions,” §II, Eqs. (15)–(21), which derives the source-free tensor equation and the radiation/matter Bessel solutions. Those papers assume exact constant-equation-of-state eras and a homogeneous source-free mode. They establish the analytic basis, not the Fast state conversion or accuracy of a discretized reheating match; those were checked separately against DOP853 and an exact fractional-event unit case.

The unit screen first failed on the old frozen-z propagation, then passed after the Bessel transfer was implemented. A paired end-to-end A/B with 22 eligible modes (`cr=0`, `T_re=2000 GeV`, `kappa10=.01`, `n_t=1`, `A_s=1e-15`, `DN_re=740.003`) changed the self-consistent root by `19.30` on a root near `1.50e7`, and the interpolated same-grid spectrum by at most `1.13e-6 dex`. The A/B used the same source and disabled only the new shortcut while retaining the exact event metadata; it is a local branch-isolation probe, not an old-binary comparison. Four paired warm runs gave P50/P95 `74.94/80.77 ms` with the transfer and `76.30/81.49 ms` without it. This is too little numerical change to justify the added production path, so the candidate and its unit tests were removed.

The `395`-frequency feature-band point had previously appeared to fail accuracy against Reference `z_tail=5`: Fast root `0.1613390744` versus the closed shallow-tail root `0.1614848239` (`9.03e-4` relative), with maximum spectrum difference `0.004664 dex`. A same-input, same-grid Reference check at `z_tail=8` gives map `0.1613459923` at the Fast root and maximum/P95 spectrum differences `0.0009420/0.0006924 dex`. One updated-root map gives `0.1613462882` with relative residual `1.83e-6`; the Fast root is `4.47e-5` relative from this closed deep-tail estimate. At the shared node `log10(f/Hz)=-16.186935755`, the Fast error is `-0.004664 dex` relative to Reference `z_tail=5`, and `-3.45e-5 dex` relative to `z_tail=8`; the Reference spectra themselves shift by `-0.004629 dex`. The z8 solve does not use a tail for that mode. Independent shared-node trajectory checks up to `z=5` for two modes show maximum Fast/Reference state differences `7.11e-5` and `7.95e-5`, with maximum `z` differences below `2.7e-7`. Thus the earlier feature-band rejection at this point was driven mainly by the shallow Reference tail convention; it was not evidence that an exact matter-era Fast transfer was needed. The point was already a Fast success, so this revises its accuracy classification but does not count as a new Fast capability.

The wide-support `T_re=1e-3 GeV, kappa10=1e-3` probe remains unresolved: its reheating event landed on a background grid node, and prior Reference checks had up to `59.6%` map discrepancy. The feature-band z8 result does not validate that separate region. Its previous z5 measurements must be reconsidered only if an independently sourced shallow-tail discrepancy is present there. The new per-mode trajectory and root-closure artifacts are `mode_trajectory_shared_nodes.json`, `feature_reference_root_z8_next_map.json`, and `exact_matter_transfer_feature_probe.json`; scripts and hashes accompany them.

## Literature and scope

The accepted background rewrite follows directly from the project’s piecewise matter-dominated horizon equation and its breakpoint condition. The state rescaling is a gauge change of a linear two-component initial-value problem, with the mode-specific factor fixed by its initial node. No new numerical-analysis method or paper-derived approximation was introduced. `Reference` source and acceptance thresholds remain unchanged.

## Current decision

The long-background overflow fix passes the sampled narrow-cutoff `cr=0` rectangle and Reference root checks. One previously rejected `395`-frequency point now passes the root and spectrum gates against a converged `z_tail=8` Reference; this was an accuracy-reference correction for an existing Fast success, not new domain coverage. The wide-support grid, `cr>0` expansion, and complete delivery gates remain unresolved. The exact matter-era Bessel candidate was rejected. Do not treat the remaining failed wide-spectrum regimes as accepted.

## 2026-11-19 continuation — stable H2 primitives and WKB screening

### Mathematical and implementation change

On each smooth post-reheating branch, the background definitions imply

`sigma = 1 - (d ln H2 / dN)/3` and `F(N)-F(N_re) = (N-N_re) - (ln H2(N)-ln H2(N_re))/3`.

The candidate evaluates this exact identity from positive finite `H2`, anchored at the unmodified physical reheating event. Before reheating it retains the exact `sigma=1` matter branch. No transition is smoothed or moved. If `Psi` or either required exponential scale is outside float64 representability, the implementation uses a common relative scale and the same H2 primitive; it reports a nonfinite or nonpositive H2 when the mathematical background is invalid. It does not clamp invalid values.

The high-sigma panel trigger `h*|sigma_right-sigma_left| > 0.05` is a local resolution indicator. It is independent of parameter names and also protects finite-but-unresolved terminal panels. The WKB handoff screen uses the local small parameter `epsilon=|1.5 sigma-1| exp(-z)`. Its fused Numba implementation uses a suffix maximum and monotone-`Phi` binary search; if sampled `Phi` is not monotone it uses the general scan. Before building the suffix array, a sufficient global bound `max|1.5 sigma-1| * exp(-(min_tail - max_delta_Phi)) <= 0.02` skips the detailed scan for smooth trusted backgrounds. This bound includes the predecessor node used by the transfer kernel. Scalar-equivalence and terminal-unsafe regressions pass.

This change follows from the project’s own `H2_vec`/`sigma_vec` equations and the linear tensor system. No external approximate solver or paper-derived method was adopted; the identity is an algebraic consequence of the background definitions. The assumptions are a smooth branch with positive finite H2 and the exact unchanged reheating breakpoint.

### Paired domain evidence

All screens used the same input objects, Fast configuration (`h=0.005`, `z_tail=5`, goal frequency grid, two threads, `tol=1e-6`, `physical_guard=True`), and baseline source `872a70d`:

- `cr=0`, `Omega_bh2={1e5,3e5,1e6}` by `T_re={1000,2000,4000}` plus two independent holdouts: baseline `1/11`, candidate `11/11`. Baseline failures were nine `nonfinite_background_primitive` and one 15-second resource timeout; the one prior success remains successful.
- `cr=0`, `Omega_bh2={1e5,3e5,1e6}` by `Omega_ch2={1e5,3e5,1e6}` at fixed `T_re=2000` plus two holdouts: baseline `0/11`, candidate `11/11`; all baseline failures were `nonfinite_background_primitive` with positive finite H2.
- A separate `cr=1` consistency-relation density rectangle used the model-derived tilt (`n_t=-0.00125`) and required a finite derived `N_inf=68.81384` plus physical-guard acceptance. Baseline `0/11`, candidate `11/11`, all old failures the same primitive classification. This is mathematical-computation coverage only; the high-density corner fails the accuracy gate described below and is not counted as accepted expansion.

The two `cr=0` slices have direct same-grid `z_tail=8`, `rtol=1e-9` Reference checks at both corners and an independent holdout in each slice, plus an old-success control and a former baseline timeout point. The maximum selected root-relative proxy is `1.5121e-4`; maximum full-spectrum error is `8.7021e-4 dex`. Current-source replay after fusing the WKB screen reproduced all eight stored Fast roots, frequency arrays, and spectra exactly. Finite samples support local region improvement; they do not prove a continuous parameter domain.

For `cr=1`, the low-density corner and independent holdout pass the root and spectrum gates, but the high-density corner does not: `z_tail=8` root-relative proxy `4.3667e-4`, spectrum max `9.5425e-4 dex`. One same-grid `z_tail=10` check leaves the root proxy at `4.3525e-4` and gives max spectrum error `1.0938e-3 dex`; extending the tail does not resolve it. A 254-node adaptive Fast grid also leaves the root proxy at `4.2471e-4`. Therefore the `cr=1` successes are not credited as accurate domain expansion; high-density `cr>0` remains a known boundary.

### Performance and artifacts

On the same Windows host, two Numba threads (`workqueue`), 64 alternating warm pairs on each of four trusted full Fast inputs, roots, grids, warnings, and spectra were exact. Aggregate candidate/baseline ratios were P50 `0.8850` and P95 `0.8603`; all four individual cases also stayed below the 5% regression limit. Cold first-solver P50 was `0.5935/0.5745 s` (candidate/base), and process-to-first-solve P50 was `1.644/1.599 s`. In the newly covered rectangles, candidate Fast solve times were `0.094–0.193 s`; old Fast failed on those cases. The independent Reference calls remain slow and were used only for selected precision checks.

Paired screens, raw same-grid Reference spectra, the exact-replay artifact, benchmark samples, z8/z10 CR1 checks, runner hashes, and [`h2_identity_domain_coverage.png`](h2_identity_domain_coverage.png) are in this directory. [`plot_h2_identity_domain.py`](plot_h2_identity_domain.py) regenerates the failure-map and precision figure. The source hashes are stored with each data artifact; the supplementary replay artifact links earlier CR0 Reference arrays to the current candidate by exact output equality.

### Decision

Accept the stable H2 identity and fused WKB screen for the sampled `cr=0` expanded regions only. The final trusted samples show exact output with aggregate P50/P95 ratios `0.8850/0.8603`; all four individual cases are below 5%. The `cr=1` high-density corner remains outside accepted accuracy coverage because it misses the root gate and its z10 spectrum misses the spectrum gate. Final local gates passed: `428 passed, 9 deselected`, Ruff, Mypy, compileall, sdist/wheel build, and installed-wheel smoke. Commit `5895d29ddbb685f5ed3a04574bf87565cd91e249` was pushed to `fast_v0.2`; GitHub Actions run `44541293060` passed all nine jobs.

## 2026-11-19 final-spectrum representability classification

The reconstructed wide-support probe (`cr=0`, `T_re=1e-3`, `kappa10=1e-3`, `n_t=0`, `r=0.01`, `A_s=2.100549e-9`, `DN_re=717.4979895910536`) had returned a converged model while one final linear spectral bin was zero and `log10OmegaGW` was `-inf` with a divide-by-zero warning. This is a nonzero tensor input whose final linear output is not representable; it is not an exact-zero limit. The Fast output assembly now fails closed with `nonrepresentable_output/final_spectrum_underflow` for a zero final signal in the linear-amplitude path. Nonfinite and other nonpositive signals retain separate failure reasons; the explicit exact-zero path and log-amplitude recovery are unchanged. The same regression was red before the change and passes afterward.

This is an output-validity correction only, not an expansion result: the wide-support parameter region remains unresolved and is not counted as Fast-success coverage. Verification: `423 passed, 9 deselected` (two existing deprecation warnings), plus the focused affected tests and Cobaya marker; Ruff, Mypy, compileall, diff check, sdist/wheel build, distribution boundary check, and installed-wheel smoke pass. GitHub Actions run `44474144299` passed all 9 jobs, including canonical regression/manifest, static checks, packaging, Cobaya, and Python 3.9–3.13 compatibility. The changed `fast_sgwb.py` SHA256 is `29a398be7bcfa112df0705c8c56d3b14b372ebe5f2fdb274ca378b2edd9db1fc`. Prior expansion hashes in `source_hashes_20261119.json` still refer to their original benchmark source state.



## 2026-11-19 continuation — precision-qualified cr=1 density rectangle

### Accepted finite-sample coverage

The earlier `cr=1` 3×3 rectangle extending to `Omega_bh2=Omega_ch2=1e6` remains outside the accepted precision domain because its upper corner misses the root gate. The precision-qualified square now reaches `Omega_bh2, Omega_ch2 ∈ [1e5,6e5]`, with `cr=1`, consistency-derived tilt `n_t=-0.00125`, finite `N_inf=68.81384`, and `physical_guard=True`:

- The paired 3×3 lattice on `{1e5,5e5,6e5}` has old Fast `0/9` and current Fast `9/9`; all old failures are `nonfinite_background_primitive`. It combines earlier same-input rows with the new 600k cross points. Base SHA is `872a70dc5ab3bb2ff9d3c0d78f05850341b42a45`, candidate source SHA `b236c5632fa06cbd113f11627efadee92947c243`.
- Independent paired holdouts `(2e5,2.5e5)` and `(5.5e5,5.5e5)` are both old Fast `0/1`, current Fast `1/1`, with the same old failure. The square plus these holdouts samples `0/11 → 11/11`.
- Same-grid z8 Reference checks cover each distinct square total `200k,600k,700k,1M,1.1M,1.2M`. The overlapping asymmetric rectangle adds totals `400k,800k`, and the independent lower holdout adds `450k`. Across nine direct Reference checks plus two exact-sum holdout reuses, maximum root-relative proxy is `1.4193e-4` and maximum spectrum error is `9.0781e-4 dex`, below `2e-4` and `1e-3 dex`.
- `Omega_bh2` and `Omega_ch2` enter downstream only through `Omega_mh2=Omega_bh2+Omega_ch2`. Swapped coordinates with the same total density produce identical modeled backgrounds; each distinct total used in these samples has direct Reference coverage. This source-level symmetry does not certify unsampled points.
- New 600k square points run in roughly `0.09–0.15 s`; the 1.1M and 1.2M direct Reference calls take about `71 s` each. The `7e5 × 7e5` probe also succeeds in Fast, but two same-grid z8 Reference maps close at residual `8.81e-8` and leave the Fast root `2.0355e-4` from the closed map, just beyond the `2e-4` gate; its spectrum error `8.696e-4 dex` passes. That point is excluded, as is the previously screened `1e6 × 1e6` corner. This brackets the current square's accuracy boundary with a paired Fast success that is not counted as a valid solution. No interpolation is used to claim unsampled accuracy.

### Additional asymmetric CR1 rectangle

A second overlapping sampled rectangle extends the baryon-density axis while keeping the total density within the Reference-validated range: `Omega_bh2 ∈ [1e5,7e5]`, `Omega_ch2 ∈ [1e5,3e5]`. On the `{1e5,3e5,5e5,7e5} × {1e5,3e5}` lattice, combining prior paired rows and new exact baseline screens gives old Fast `0/8`, current Fast `8/8`; every old failure is `nonfinite_background_primitive`. The independent holdout `(6.5e5,1.5e5)` is old `0/1`, current `1/1`.

The lattice totals are `200k,400k,600k,800k,1M`, already covered by the direct z8 Reference checks above. For the independent holdout, the current Fast frequency nodes match the saved 800k Reference nodes bitwise; its root proxy is `3.4397e-5` and spectrum max error is `8.7012e-4 dex`. This reuse is valid because both density inputs enter only through their exact sum, and the source hash matches. The holdout Fast run took `0.35 s`; new rectangle corners took `0.15–0.16 s`.

The square and asymmetric rectangles overlap but are distinct sampled regions. They are finite coverage evidence, not a guarantee for unsampled points or a continuous-domain proof. The updated plot outlines both regions.

### Extended CR1 rectangle through 1.3M total density

The validated square and asymmetric rectangle are contained in the larger rectangle `Omega_bh2 ∈ [1e5,7e5]`, `Omega_ch2 ∈ [1e5,6e5]`. A 3×2 lattice on `{1e5,6e5,7e5} × {1e5,6e5}` has old Fast `0/6` and current Fast `6/6`; all old failures are `nonfinite_background_primitive`. Three paired holdouts inside the rectangle `(2e5,2.5e5)`, `(6.5e5,1.5e5)`, `(6.5e5,5.5e5)` also change `0/3 → 3/3`.

Direct same-grid z8 Reference checks cover the new total densities `700k` and `1.3M`; existing raw data cover `200k`, `800k`, and `1.2M`. The maximum root proxy across accepted direct checks is `1.6726e-4`; maximum spectrum error remains `9.0781e-4 dex`. The `1.2M` holdout reuses the same-sum Reference with a bitwise-equal frequency grid. Fast corners take about `0.14–0.16 s`.

The adjacent `7e5 × 7e5` point has a paired Fast success, but two Reference maps close at residual `8.81e-8` and place the Fast root `2.0355e-4` from the closed Reference root. Its spectrum error passes (`8.696e-4 dex`), but the root gate does not, so that point remains excluded. The current rectangle is supported by finite samples; no continuous-domain claim is made.

### Rejected local transfer-curvature candidate

A separate diagnostic found that high-error CR1 mode trajectories accumulate transfer-state discrepancy near the terminal active panels. A quadratic quarter-point local transfer reduced the largest discrepancies for two selected modes, but a full high-corner Fast solve moved the root only from `0.2117557174` to `0.2117574375`; relative to the saved same-grid z10 Reference map this remains `4.43e-4`, slightly worse than the previous `4.35e-4`. The candidate production edit and its focused test were removed. This mechanism is not accepted and provides no capability or performance claim. Its fixed-root mode traces and the no-change Fast probe remain in the diagnostic artifacts.

The updated figure marks the precision-checked lower `cr=1` rectangle and includes its Reference points. Raw spectra, paired baseline/candidate rows, diagnostic scripts, and source hashes are indexed by `cr1_low_rectangle_manifest.json`.
