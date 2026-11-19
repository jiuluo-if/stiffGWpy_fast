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

The unit screen first failed on the old frozen-z propagation, then passed after the Bessel transfer was implemented. A paired end-to-end A/B with 22 eligible modes (`cr=0`, `T_re=2000 GeV`, `kappa10=.01`, `n_t=1`, `A_s=1e-15`, `DN_re=740.003`) changed the self-consistent root by `19.30` on a root near `1.50e7`, and the interpolated same-grid spectrum by at most `1.13e-6 dex`. The A/B used the same source and disabled only the new shortcut while retaining the exact event metadata; it is a local branch-isolation probe, not an old-binary comparison. Four paired warm runs gave P50/P95 `74.94/80.77 ms` with the transfer and `76.30/81.49 ms` without it. This does not repair the feature-band errors: on the prior `395`-frequency point, whose reheating event was grid-aligned and therefore did not invoke the shortcut, the spectrum remained `0.004664 dex` off Reference and the root error remained about `9.03e-4`. The wide-support probe also had a grid-aligned event. The candidate was removed from production and its unit tests; the probe and A/B artifacts are retained locally for provenance.

This result redirects work toward the spectral-error source. In the wide-support probe, the reheating event landed exactly on a background grid node, so the fractional-event shortcut was correctly ineligible; the frequency grid itself had eight nodes from reheating through the next two dex. In the separate `395`-node feature-band comparison, propagation error persisted even though the grid was dense. More work is needed to distinguish the post-reheating transfer residual from background/initial-state effects before choosing another production candidate.

## Literature and scope

The accepted background rewrite follows directly from the project’s piecewise matter-dominated horizon equation and its breakpoint condition. The state rescaling is a gauge change of a linear two-component initial-value problem, with the mode-specific factor fixed by its initial node. No new numerical-analysis method or paper-derived approximation was introduced. `Reference` source and acceptance thresholds remain unchanged.

## Current decision

The long-background overflow fix passes the sampled narrow-cutoff `cr=0` rectangle and Reference root checks. This evidence does not support broad-spectrum or transition-in-support accuracy, a `cr>0` expansion, or the complete delivery gates. The exact matter-era Bessel candidate was rejected. Do not treat the failed wide-spectrum or transition-in-support regimes as accepted.

