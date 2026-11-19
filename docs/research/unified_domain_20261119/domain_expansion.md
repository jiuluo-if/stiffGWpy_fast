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
- A 3x2 wide-spectrum slice that brings the transition into the sampled support fails the accuracy gate. At two points the 127-node Reference maps differ from the default-grid Fast roots by 31% and 60%. An independently closed 395-node Reference root also differs from Fast by `9.03e-4` on another transition-in-support point. These results are retained as rejected regimes.
- The `cr>0` 3x3 control grid is valid under the derived consistency tilt and cutoff constraints, but baseline also succeeds there. It is a regression control, not an expansion result.

Artifacts in the experiment root include the exact baseline and candidate grid rows, frequency-refinement records, Reference root closure, per-mode error table, timing screen, and baseline/candidate source hashes. `domain_coverage_20261119.png` visualizes the sampled grids.

## Performance and verification

An interleaved trusted-point benchmark used eight alternating fresh processes, each with one warmup and 12 timed calls (96 samples per version). Baseline/candidate P50/P95 were 10.34/12.33 ms and 10.16/11.58 ms. The candidate shows no measured regression on this workload. The separate cold-to-warm process startup medians were 262/258 ms. New-region calls on the 3x3 narrow-cutoff slice took 46–47 ms warm after compilation. The broad feature-support cases are rejected on accuracy.

The added regression test checks the 3x3 rectangle. The focused file passes 11 tests. The current full suite passes 422 tests with 9 deselected. Ruff, Mypy, source distribution and wheel builds also pass. GitHub CI has not yet been run for this candidate.

## Literature and scope

The accepted background rewrite follows directly from the project’s piecewise matter-dominated horizon equation and its breakpoint condition. The state rescaling is a gauge change of a linear two-component initial-value problem, with the mode-specific factor fixed by its initial node. No new numerical-analysis method or paper-derived approximation was introduced. `Reference` source and acceptance thresholds remain unchanged.

## Current decision

The long-background overflow fix passes the sampled narrow-cutoff `cr=0` rectangle and Reference root checks. This evidence does not support broad-spectrum, transition-in-support accuracy, `cr>0` expansion, cold/warm P50/P95 qualification, or the complete delivery gates. Do not treat the failed wide-spectrum or transition-in-support regimes as accepted.

