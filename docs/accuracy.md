# Accuracy and evidence boundaries

Status: current summary of scoped evidence; no full parameter-space release certification

Date: 2026-10-23

## Exact horizon-crossing step boundary (2026-10-23)

In `scaled_step`, `w² = exp(2 z_mid)` equals one at the horizon-crossing
boundary. The oscillatory formula `sin(Ωh)/Ω` is undefined there because
`Ω=0`; the continuous limit is finite (`cos=1`, `sin(Ωh)/Ω=h`). Route equality
through the branch that evaluates this limit directly. A compiled direct probe
and the `_phase_segment` path now return finite values matching the limit
exactly. After selecting the formal `fast` preset explicitly, default, low-T,
high-T, stiff, and high-kappa retain all five output digests from the preceding
profile. The direct boundary values and five-case digest comparison are
recorded in [the boundary-fix result](scaled_step_horizon_boundary_20261023.json).
Those five-case digests were recorded for the preceding algebraically equivalent
series-limit implementation. The exact boundary is now written explicitly as
`c=1`, `si=h`; the artifact records the current source hash and a fresh direct
compiled probe. The five full cases were not rerun after this expression rewrite.
This edge correction adds no new Oracle coverage or general accuracy claim.

The formal fast profile is the single `goal-kink-hybrid` mode. The independent
continuous-sigma DOP853 pipeline is the precision reference; LSODA is a
regression/runtime path, not a truth oracle. Every numerical result below is
bounded by the artifact's commit, reference convention, parameter points, and
frequency grid.

## Earlier named-point baseline (2026-09-12)

The four-regime quadrature audit below is a dated baseline from
[`fast_v02_audit_report.md`](fast_v02_audit_report.md), not the latest Oracle C
edge evidence or a full parameter-space certification. The newer two-edge
`DN_re` results are recorded in the
[matched Oracle C assessment](oracle_c_wkb_assessment.md) and the dated section
below.

| Comparison | Recorded result | Scope and interpretation |
|---|---:|---|
| PCHIP fast vs Oracle C WKB `DN_gw` | `1.07e-5`–`1.66e-4` | Four named points; below the `2e-4` gate on those points |
| Fast vs reference, same native grid, reference `z_tail=8` | median `2.93e-4`, max `2.96e-4` | Six named points; includes the reference frozen-tail convention |
| Fast vs reference, default point, reference `z_tail=10` | `4.96e-6` | One-point deep-tail attribution, not a parameter-space result |
| Nested native-frequency refinement | relative change `7.4e-12`–`2.8e-9` | Six named points; supports grid convergence in that test |

The same audit reports Simpson errors above the PCHIP results and records
PCHIP as the formal default. Its `z_tail=8` reference difference is limited by
the reference's own frozen-handoff tail sensitivity; the one-point `z_tail=10`
comparison is an attribution check, not a substitute for broader validation.

## 2026-10-22 matched `DN_re` edge probe

A full-grid comparison adds two reheating-parameter edges to the three
previously sampled Oracle C points. Its original edge anchors were generated
from solver commit `5e6e78875cf9998fb64b72c6b216c32c7848abe8`; the matched fast
comparison used `07293069a64c16d9efc7d7896590798b76d9591a`. Oracle WKB and
deep-reference `DN_gw` agree within `2.4e-6` and `5.7e-6`, respectively. The
formal fast profile is within `1.2e-5` on the earlier three points, but differs
by `2.74e-3` (`edge_dnre_lo`) and `2.50e-3` (`edge_dnre_hi`) on the
new edges. Those two results fail the scoped `2e-4` gate; this supersedes any
interpretation that the gate has been verified across representative parameter
edges. It remains a sampled probe, not a full parameter-space certification.

With Oracle anchors regenerated at the same `DN_eff` and native frequency
grid, reducing `h` to `0.0003125` lowers the two errors to `1.65e-4` and
`1.63e-4` versus WKB (both below the scoped gate). Warm median runtimes rise
from `5.69`/`4.53 ms` at `h=0.005` to `62.7`/`71.7 ms` (about `11.0x`/`15.8x`).
The matched `transition-refine`, `h=0.0025` candidate also passes these edges,
at `3.97e-6` and `1.64e-4` versus deep `z=10`; its medians are `46.5`/`73.7 ms`
(about `8.2x`/`16.3x` slower than standard fast). It is faster than uniform
fine-step only on `edge_dnre_lo`. Neither candidate is promoted; improving
latency remains open. Full grids, configuration, repeated timings, and reproduction commands are in
[`oracle_c_wkb_assessment.md`](oracle_c_wkb_assessment.md) and
[`oracle_c_wkb_corner_20261022/`](oracle_c_wkb_corner_20261022/).

Reducing the transition-refinement radius from `2.0` to `0.25` e-folds while
keeping `h=0.0025` and factor `8` lowers the matched medians to `35.6`/`51.7 ms`
and improves the two errors to `5.17e-6`/`1.62e-4` versus deep `z=10`. It is
about `1.29x`/`1.42x` faster than the radius-2.0 candidate, while still about
`6.3x`/`11.4x` slower than standard fast. This two-edge candidate remains
validation-only. Its exact configuration and raw records are in
[`oracle_c_wkb_corner_20261022/transition_h00025_r025/`](oracle_c_wkb_corner_20261022/transition_h00025_r025/).

The radius `0.125` sweep holds the local step at `0.0003125` while increasing
the outer step from `h=0.0025` to `0.02`. The `h=0.02/factor=64` candidate
passes both tested edges at `1.16e-4` and `8.25e-5` versus deep `z=10`, with
warm medians `7.63/9.49 ms`. This is the fastest passing point in the sweep,
about `1.34x/2.09x` slower than standard fast. Coarsening further to
`h=0.04/factor=64` reduces runtime but fails both edges (`4.49e-4/3.17e-4`).
These candidates remain validation-only and are documented with their matched
Oracle anchors in [`oracle_c_wkb_assessment.md`](oracle_c_wkb_assessment.md).

At `h=0.04`, widening the radius to `0.25` and increasing the factor to `128`
restores the same `0.0003125` local step, but only the high `DN_re` edge passes:
the low/high errors are `4.51e-4` and `1.38e-4`. Exact anchor checks pass, so
the low-edge miss is not a grid mismatch. This candidate is rejected; its
factor-64 comparison at the same radius fails both edges (`4.51e-4`/`2.95e-4`).
At `h=0.02`, reducing the radius from `0.125` to `0.0625` retains the
`0.0003125` local step and passes at `1.14e-4`/`8.30e-5`, with medians
`6.49/8.74 ms`. It remains slower than standard fast and is not promoted. At
radius `0.015625`, errors are `1.12e-4/8.27e-5` with medians `6.15/7.98 ms`;
at `0.0078125`, they are `1.10e-4/8.27e-5` with medians `6.31/7.94 ms`. A
fixed-affinity paired comparison found no stable timing advantage for the
smaller radius. Matched Oracle data are linked from the assessment report.

## Current full-path solver profile is separate

The corrected 2026-10-22 formal `fast` profile used solver source at commit
`6707572e73d576a0ccf294172579cedc26e5e394`, with the preset's `tol=1e-6`,
PCHIP, goal grid, and exact kink split. It covers default, high-T, stiff, and
high-kappa with 25 warm repeats per case; the first process call is reported
separately. Warm medians are `5.987`, `10.049`, `9.179`, and `9.458 ms`. It is a
runtime profile with repeatability digests, not a matched Oracle accuracy
sweep or LSODA comparison, so it adds no precision coverage and establishes no
LSODA speed ratio. Full per-case values and resource provenance are in
[`benchmarks.md`](benchmarks.md).

The earlier 2026-10-21 profile at commit
`cdbdb387479e11cd87a38d76a0e58348ca13e70f` omitted the low-level solver's
`tol` argument and therefore used `1e-4`, not the formal preset's `1e-6`. It
remains a historical stage measurement, not the formal profile baseline.

Later PTA likelihood-cache timings are adapter-level measurements, not
full-path solver profiles or accuracy evidence; see
[`benchmarks.md`](benchmarks.md).

The 2026-10-21 tensor sample-index change matched the baseline output arrays
bit-for-bit in six benchmark regimes. This is an exact differential check, not
new oracle-accuracy evidence or a full parameter-space certification.

The 2026-10-22 `H2_vec` batching A/B also matched every solver output array
exactly on two `DN_re` edges, but it measured the validation `transition-refine`
path. It improves runtime on that path and adds no Oracle accuracy points or
evidence for the formal exact-kink profile.

## Likelihood and posterior evidence

The 2026-09-03 importance-sampling validation compares the fast spectrum with
the continuous-sigma reference on 240 posterior-bulk points and 11 likelihood
bins. It reports max per-bin dex error `3.10e-4`, max absolute `Delta logL`
`7.30e-3`, ESS `4167`, and a `log10(r)` posterior shift of `-0.0011 sigma`.
These are results of the recorded mock-data/importance-reweighting procedure,
not an independent reference-engine MCMC chain; see
[`mcmc_posterior/posterior_validation.md`](mcmc_posterior/posterior_validation.md).

The separate 2026-09-24 LVK report runs four chains for each method and
scenario. It contains convergence diagnostics and posterior-spectrum samples;
SageNet+ does not meet the report's diagnostic references in two scenarios.
The low-reheating scenario has 0% LVK-band coverage for all methods and cannot
support a data-fit conclusion. See
[`mcmc_sagenet_compare/report.md`](mcmc_sagenet_compare/report.md).

## Local error budget

`estimate_local_error` distinguishes telemetry measured during the solve,
terms calibrated from a fiducial point, and uncertified defaults when no solve
telemetry is attached. A `certification_status` of
`certified-fiducial-calibrated` is limited to that calibration model; it does
not imply a uniform error guarantee over all parameter values.

The formal `kink_split=True`, `transition_refine=False` path and the
`freq_grid='goal'` grid do not have matching calibration anchors in the
historical error-budget experiments. Their existing proxy values remain
available for diagnostics, but the affected terms appear in
`uncertified_terms` and the status is `uncertified-path-configuration`.
When `auto_escalate=True`, this status forces a higher-accuracy solve before
the likelihood gate accepts the fast result. The default fast path is unchanged
when auto-escalation is disabled.

The `self_consistency` term is the last-step change in `DN_gw` divided by the
current `DN_gw` magnitude. At low signal amplitude this relative diagnostic can
be large even when the absolute change is small; it is an iteration-change
proxy, not an error bound.

正式的 `kink_split=True`、`transition_refine=False` 路径以及
`freq_grid='goal'` 尚无对应的历史误差预算校准样本。现有代理值仍供诊断，
但相关分项会列在 `uncertified_terms` 中，状态为
`uncertified-path-configuration`。启用 `auto_escalate=True` 后，求解会先升级精度，
再由似然门限决定是否接受；关闭自动升级时，默认快速路径行为不变。

`self_consistency` 表示最后两次 `DN_gw` 的差值相对于当前 `DN_gw` 幅值的比值。
信号幅值很小时，即使绝对变化较小，该相对诊断值也可能很大；它表示迭代变化，不是误差上界。

## Accepted execution optimization evidence

The dated H2 endpoint-cache A/B recorded identical spectrum and background
output digests over its 13 named/Sobol/edge points and 50 repeats. The measured
high-kappa warm median improved by 6.54%; the low-T case regressed within the
reported noise range. This supports that specific cache change and protocol,
not a universal runtime improvement. Consult the linked artifact before
reusing the result.

## 2026-10-23 transition-refinement frontier

At `h=0.0296875`, radius `0.015625`, and factor `95`, the candidate passes both
sampled `DN_re` edges. Errors versus deep `z_tail=10` are `8.59e-5` and
`9.10e-5`; the matched Oracle WKB/deep residuals are `2.88e-6` and `5.84e-6`.
The tested `h=0.03` neighbor failed the high edge at `2.49e-4`. A fixed-affinity
30-pair A/B found mixed timing: candidate/formal-fast paired median ratios are
`0.812` (95% interval `0.792–0.837`) on the low edge and `1.177`
(`1.158–1.200`) on the high edge. The sampled accuracy improves while latency
trades between edges; this is not a general speedup or a full-space guarantee.
The candidate remains validation-only. Matched records are in
[`oracle_c_wkb_corner_20261023/transition_h00296875_r0015625_f95/`](oracle_c_wkb_corner_20261023/transition_h00296875_r0015625_f95/)
and raw paired timings are in
[`benchmark_transition_frontier_ab_20261023.json`](benchmark_transition_frontier_ab_20261023.json).

## 2026-10-30 kink-path terminal step correction

`gen_fast` appends the exact present-day endpoint `N_inf`. When `kink_split`
was enabled, the propagation call still received `h_arr=None`, so its scalar
step `h` was also used for the shorter final interval. The fix passes that
single terminal width as a scalar; general transition-refined grids continue
to use their per-interval step array. It avoids a per-step array lookup in the
ordinary kink path.

The red regression records the actual final interval and verifies that the
kernel receives it. A six-point, two-thread, 15-pair A/B found the scalar path
equivalent to the general array-step path within `2.5e-14 dex` maximum spectrum
difference and `<1e-14` relative `DN_gw` difference. Against the legacy scalar
step path, a 101-pair replay on high-T, stiff, and low-r measured P50/P95
changes of `-0.16%/-1.06%`, `-0.27%/+1.78%`, and `-0.98%/-0.16%`.

For an independent accuracy check, the four 76-frequency Oracle C records use
WKB-corrected `z_tail=5` values where available and deep `z_tail=10` values for
the remaining modes. The comparison uses `log10(Ogw_today - Oj_today)`, matching
Fast's public spectrum convention. Results are bound to Oracle source commit
`1bf840d` and Fast commit `d3fff31`:

| Point | Fast spectrum P95/max error before | After | Fast `DN_gw` relative error vs WKB after |
|---|---:|---:|---:|
| Default | `4.05e-4 / 5.45e-4 dex` | `2.98e-4 / 4.16e-4 dex` | `1.07e-5` |
| Low-T | `1.03e-3 / 1.84e-3 dex` | `2.56e-4 / 3.91e-4 dex` | `6.25e-7` |
| High-T | `4.12e-4 / 9.23e-4 dex` | `2.55e-4 / 3.87e-4 dex` | `9.88e-6` |
| Stiff | `9.03e-4 / 2.97e-3 dex` | `2.53e-4 / 3.93e-4 dex` | `1.15e-5` |

All 76 positive modes at each point are below `5e-4 dex` after the correction;
the corresponding Oracle WKB-to-deep `DN_gw` differences are at most
`1.25e-5`. Raw pointwise comparisons and provenance are in
[`kink_terminal_step_fix_oracle_c_20261030.json`](kink_terminal_step_fix_oracle_c_20261030.json),
reproduced by [`benchmark_kink_terminal_step_fix.py`](../scripts/benchmark_kink_terminal_step_fix.py).
This is four named-point evidence, not a global parameter-space certification.

The same fix also recovers the retained random-corpus sample11 at
`h=.005`, `freq_grid='grid_independent'`, and `z_tail=8` on the 94 native
tail10 Reference frequencies. Corrected Fast gives `DN_gw=4.85759964`, about
`2.62e-5` relative from the separately iterated Reference root
`4.85747263`; median/P95/max spectrum error is `4.01e-5/1.43e-4/1.65e-4 dex`.
The legacy terminal step had P95/max `5.17e-4/6.33e-4 dex` and seven modes
over `5e-4`. At `z_tail=5`, a `1.70e-3 dex` worst-bin error remains because
the shallower frozen-tail convention dominates that bin; the endpoint fix
does not remove this separate tail error. The underlying sample parameters
and earlier transition-refined recovery are recorded in
[`fast_domain_pressure_20261030.md`](fast_domain_pressure_20261030.md). The
endpoint-path A/B is reproducible with
[`benchmark_sample11_terminal_step.py`](../scripts/benchmark_sample11_terminal_step.py)
using the saved 94-mode Reference spectrum
[`sample11_reference_tail10_spectrum_20261030.json`](sample11_reference_tail10_spectrum_20261030.json)
and its second-map fixed-point check
[`sample11_reference_second_map_20261030.json`](sample11_reference_second_map_20261030.json).

The separate blue-tilt res4 `h=.02, transition_refine=True` candidate's
independent fixed-background Reference and local secant-root check are now
complete; see the 2026-11-04 record below. This is a separate parameter-point
check and does not expand the four-point Oracle C scope above.

## 2026-11-04 res4 aligned-transition scalar-step path

The optional `transition_refine=True`, `h=.02`, zero-radius path uses an exact
reheating node, but its grid is otherwise uniform for this parameter point.
The kernel now detects that shape and uses scalar `h` plus the exact final
interval width. If any true interior refinement or near-coincident physical
anchor exists, it retains the full interval-width array. The optimization is
guarded by a floating-point roundoff check and does not change the formal
default grid. Code and regression tests are committed at `671aa11`.

On the blue-tilt res4 point (`r=.01`, `cr=0`, `n_t=2`, `T_re=1e4 GeV`,
`kappa10=.01`, `DN_re=10`) and 260 common native frequencies:

- The refreshed uniform `h=.005` full-root spectrum has maximum tail9 error
  `6.455e-4 dex`; 18 modes exceed `5e-4 dex` and none exceed `1e-3 dex`.
- The `h=.02`, zero-radius refined candidate has Fast root
  `DN_eff=0.13372910340092184`. Independent tail9/tail10 Reference maps at
  that exact input differ from the Fast map by `4.622e-5` and `5.647e-5`
  relative. Its spectrum maximum error is `3.4487e-4 dex` against tail10, with
  zero modes above `5e-4 dex`.
- The independent Reference fixed-point root is closed at
  `DN_eff=0.13373681807168084` (input `0.13373681817755848`, relative map
  residual `-7.92e-10`). The Fast root differs by `-5.77e-5` relative. At the
  Reference root, Fast's spectrum median/P95/maximum error is
  `6.109e-5 / 3.162e-4 / 3.4532e-4 dex`.
- Reference tail9-to-tail10 spectrum change on the same 260 frequencies is
  median/P95/maximum `9.94e-8 / 7.282e-5 / 2.0604e-4 dex`. Reference estimates
  at the closed root are `1.84e-20` quadrature error and `1.82e-10`
  interpolation error.
- A fresh fixed-background factorial decomposition at the candidate
  `DN_eff=0.13374433666555344` puts the maximum tail9 spectrum error of the
  uniform fast primitive, uniform `sigma_exact`, and transition-refined path
  at `4.1984e-4 dex`; none of the 260 modes exceeds `5e-4 dex`. The old
  53-mode `>1e-3 dex` record is stale for current source. Comparing complete
  roots instead, the default root has 18 modes above `5e-4 dex` while the
  candidate has none. This isolates the remaining default-root discrepancy
  mainly to the outer `DN_eff` map/root, rather than the continuous primitive,
  tensor propagation, or analytic-tail handoff. The 260-node Reference's own
  tail9-to-tail10 change and quadrature/interpolation estimates above bound
  the Reference contribution separately.
- Per-mode current full-root and fixed-background rows are archived in
  [`res4_mode_errors_full_root_20261104.csv`](res4_mode_errors_full_root_20261104.csv)
  and [`res4_mode_errors_fixed_background_20261104.csv`](res4_mode_errors_fixed_background_20261104.csv).
  The latter separates uniform fast primitive, uniform `sigma_exact`, and
  transition-refined propagation on the same DN and native grid. The raw
  same-background tail9 Reference is in
  [`res4_reference_transition_tail9_20261104.npz`](res4_reference_transition_tail9_20261104.npz).
- An independent DOP853 integration of the continuous-coefficient tensor
  system on the 18 former edge modes agrees with the refined Fast values to
  maximum `2.56e-6 dex`; the two DOP853 tolerance settings differ by at most
  `1.40e-6 dex`. Both methods stay within `2.5e-5 dex` of the tail9 Reference
  on those modes. They do not trigger the analytic tail handoff. Raw per-mode
  probes are [`historical53`](res4_dop853_historical53_20261104.json) and
  [`current18`](res4_dop853_current18_20261104.json); the 53-mode cohort is
  historical and is included to check those former outliers against the
  continuous-coefficient solver, not to imply they remain current failures.
- A targeted frequency-density check held `h=.005` and `z_tail=8` fixed while
  increasing `freq_res` from 4 to 8 and 16. The Fast root moved from
  `0.1337074583` to `0.1337063774`, but the same 18 modes remained above
  `5e-4 dex` and the maximum spectrum difference stayed near `6.455e-4 dex`.
  This rules out frequency-grid density as the main cause for these outliers;
  it is not a replacement for the separate background/root comparisons. Raw
  measurements are in
  [`res4_default_frequency_density_20261104.json`](res4_default_frequency_density_20261104.json).

The retained stiff random sample11 is a second independent hard point. Its
current h=.02 Fast root (`DN_eff=4.857519441699757`) was checked by a same-input
tail10 Reference map: Reference maps to `4.857473394469992`, a relative
residual of `-9.48e-6`. On the same 94 native frequencies and exact same DN,
the Reference spectrum median/P95/maximum difference is
`3.016e-5 / 6.953e-5 / 9.149e-5 dex`, with no modes above `5e-4 dex`.
Reference estimates are `6.68e-19` for quadrature and `4.42e-7` for
interpolation. Raw Fast and Reference spectra plus the map record are
[`sample11_h02_fast_20261104.npz`](sample11_h02_fast_20261104.npz),
[`sample11_h02_reference_map_20261104.npz`](sample11_h02_reference_map_20261104.npz),
and [`sample11_h02_reference_map_20261104.json`](sample11_h02_reference_map_20261104.json).
The reproducible exact-map driver is
[`benchmark_sample11_reference_map.py`](../scripts/benchmark_sample11_reference_map.py);
the 94-mode tail10 Reference run took `2158 s` with four workers on four pinned
CPUs. This expensive path is for verification only.

With 15 alternating pairs on fixed affinity, candidate/default warm P50/P95
ratios are `0.849/0.894` at 2 threads, `0.938/0.926` at 7 threads, and
`0.907/0.893` at 8 threads after the primitive-failure guard. That run measured
candidate `57.7/61.9 ms` versus default `63.7/69.4 ms`, below the 100 ms difficult-point
target. In that run, scalar steps vs the full interval array differed by
`2.78e-16` in `DN_eff` and at most `2.85e-14 dex` in spectrum, while using
`0.695/0.718` of array-path P50/P95. The corrected three-variant batch rate
was `14.52 cases/s`. A separate 10-pair resource profile
measured 14.8 cases/s and post-call RSS maxima of 159.5 MiB (default) and
157.0 MiB (candidate); these are sampled RSS values, not OS high-water marks.
The scalar-step output matches the previous full-array candidate to
`6.04e-14 dex` maximum and the root to `1.67e-15` absolute.

Reproduction uses [`benchmark_res4_scalar_step.py`](../scripts/benchmark_res4_scalar_step.py)
with the saved 260-mode Reference spectrum
[`res4_reference_secant_root_tail10_20261104.npz`](res4_reference_secant_root_tail10_20261104.npz).
The expensive independent root map can be rerun with
[`benchmark_res4_reference_root.py`](../scripts/benchmark_res4_reference_root.py);
the aggregate results and provenance are in
[`res4_scalar_step_profile_20261104.json`](res4_scalar_step_profile_20261104.json).
The raw closed-root Reference map is saved in
[`res4_reference_secant_root_tail10_20261104.json`](res4_reference_secant_root_tail10_20261104.json).
These are two difficult parameter points, not global parameter-space
certification.

## 2026-11-04 continuous-background default for transition refinement

When `sigma_exact` is omitted, `transition_refine=True` now recomputes the
background primitive from continuous piecewise `sigma`. Explicit
`sigma_exact=False` retains the previous node-based primitive. The regular
Fast path remains on the same primitive: its res4 output fingerprint is
bitwise identical to the pre-change benchmark.

On res4, this gives a full Fast root `DN_eff=0.13373208300428413`, whose
relative gap to the independently closed tail10 Reference root is
`-3.54e-5`. Across all 260 frequencies the median/P95/maximum spectrum error
is `3.914e-5 / 3.164e-4 / 3.456e-4 dex`, with zero modes above `5e-4 dex`.
The direct tail10 Reference map at this exact Fast root returns
`0.13373671803614506`, a map residual of `+3.466e-5`; the same-input spectrum
error is `3.901e-5 / 3.167e-4 / 3.458e-4 dex` median/P95/maximum. Its
quadrature/interpolation estimates are `1.84e-20` and `1.82e-10`.
The prior default h=.005 result had 18 modes above that threshold and maximum
error `6.455e-4 dex`. In 15 paired 8-thread runs, the refined continuous
candidate measured P50/P95 `66.9/69.3 ms` versus `68.5/71.5 ms` for the
default path; sampled post-call RSS maxima were `160.9/161.3 MiB`.

On sample11, the transition-refined path retains the recovery of this
historically nonfinite random case; continuous-background recomputation
further reduces its spectrum error. The candidate full root is
`DN_eff=4.857503279562753`; an independent
94-mode tail10 Reference map at that exact input is `4.8574731267571964`, a
relative residual of `-6.21e-6`. This Fast root is `6.31e-6` relative above
the independently second-map-confirmed Reference root. The same-input spectrum error is
`2.078e-5 / 5.743e-5 / 7.057e-5 dex` median/P95/maximum, with quadrature and
interpolation estimates `6.68e-19` and `4.42e-7`. The Reference run used
`rtol=1e-11` and took `2106 s` with four workers. In 15 paired 8-thread Fast
runs, transition refinement with the continuous primitive measured P50/P95
`35.7/39.5 ms`, within the `100 ms` difficult-point target.

The exact sample11 map and raw Fast/Reference spectra are
[`sample11_h02_sigma_exact_reference_map_20261104.json`](sample11_h02_sigma_exact_reference_map_20261104.json)
and [`sample11_h02_sigma_exact_reference_map_20261104.npz`](sample11_h02_sigma_exact_reference_map_20261104.npz).
The paired timing/error records are
[`sample11_sigma_exact_benchmark_20261104.json`](sample11_sigma_exact_benchmark_20261104.json)
and [`res4_sigma_exact_benchmark_20261104.json`](res4_sigma_exact_benchmark_20261104.json).
Reproduce with [`benchmark_sample11_sigma_exact.py`](../scripts/benchmark_sample11_sigma_exact.py)
and [`benchmark_res4_scalar_step.py`](../scripts/benchmark_res4_scalar_step.py).
For the exact-input sample11 comparison, pass
`--reference docs/sample11_h02_sigma_exact_reference_map_20261104.npz`
and `--reference-map docs/sample11_h02_sigma_exact_reference_map_20261104.json`.
The prior sample11 nonfinite recovery and res4 error attribution are recorded
in [`fast_domain_pressure_20261030.md`](fast_domain_pressure_20261030.md).

The direct same-input res4 map and candidate spectrum are archived as
[`res4_sigma_exact_reference_map_20261104.json`](res4_sigma_exact_reference_map_20261104.json),
[`res4_sigma_exact_reference_map_20261104.npz`](res4_sigma_exact_reference_map_20261104.npz),
[`res4_sigma_exact_same_input_20261104.json`](res4_sigma_exact_same_input_20261104.json),
and [`res4_sigma_exact_same_input_20261104.npz`](res4_sigma_exact_same_input_20261104.npz).
The 260-mode tail10 map used `rtol=1e-11`, four workers, and took `5692 s`.

## 2026-11-04 opt-in Reference Prüfer backend

`solve_reference_mode`, `spectrum_reference`, `run_reference`, and
`apply_reference_to_model` now accept `backend="prufer"`. The default remains
`backend="cartesian"`. The new backend applies the exact state change
`x=A sin(theta)`, `y=A cos(theta)` and evolves `(z, log(A), theta)`; it keeps the
same continuous background, DOP853 method, tail event, and frozen-amplitude
tail convention. Non-finite ODE results and amplitude/frequency overflow or
underflow fail explicitly. The production backend is opt-in pending validation
over additional parameter regimes.

At the independently closed res4 Reference root
`DN_eff=0.13373681817755848`, the production Prüfer backend was run on the same
260 frequencies as the archived Cartesian tail10 spectrum. Candidate settings
were `rtol=1e-9`, `atol=[1e-12,1e-12,1e-12]`, one worker, `z_tail=10`; the
Cartesian archive used `rtol=1e-11`. The candidate map is
`DN_gw=0.1337368198597089`, versus the Cartesian map
`0.13373681807168084`: relative map difference `1.337e-8`, and fixed-point
residual relative to the input root `1.258e-8`. On all 260 same-input modes the
median/P95/maximum spectrum error is `1.465e-7 / 2.602e-6 / 4.279e-5 dex`.
Candidate quadrature/interpolation estimates are `1.84e-20` and `1.83e-10`.
The single-worker candidate propagation took `657.4 s`; this is a candidate
runtime, not a paired full-grid speedup measurement.

A repeated paired screen on four representative modes (two alternating
Cartesian/Prüfer pairs per mode) used the same root background and frequency
points. The Cartesian side used `rtol=1e-11`; the candidate used `rtol=1e-9`.
Aggregate candidate speedup was `4.269x`, per-mode speedup P50/P95 was
`4.144x / 4.704x`, and the selected-mode maximum spectrum difference was
`1.401e-5 dex`. The low-frequency endpoint mode was only `2.03x` faster, so the
four-mode timing cannot be generalized to every parameter or frequency.

The root-grid result and repeat samples are archived in
[`reference_prufer_production_res4_refroot_20261104.json`](reference_prufer_production_res4_refroot_20261104.json)
and [`reference_prufer_production_ab_smoke_20261104.json`](reference_prufer_production_ab_smoke_20261104.json).
Reproduce with [`benchmark_res4_prufer_candidate.py`](../scripts/benchmark_res4_prufer_candidate.py),
using `--reference docs/res4_reference_secant_root_tail10_20261104.npz`,
`--reference-map docs/res4_reference_secant_root_tail10_20261104.json`,
`--dn-eff 0.13373681817755848`, `--rtol 1e-9`, `--candidate-only --full-grid`.
The older experimental Prüfer helper produced the initial screening artifacts;
these current figures come from the package backend.

### Cartesian DOP853 dense-output removal

The Cartesian Reference path now disables global `dense_output` and consumes
the DOP853 event state from `y_events`. This preserves the event handoff while
avoiding a dense interpolant at every accepted integration step. On two fixed
res4 modes at `rtol=1e-11`, `z_tail=10`, and identical backgrounds, all
`Ogw/Oj/Opgw` outputs matched bitwise. The no-tail/tail modes improved by
`1.275x/1.292x`; sampled peak RSS deltas fell from `81.5/94.0 MB` to
`29.4/35.2 MB`. This scoped A/B does not establish a full-grid speed ratio.
See [`reference_dense_output_benchmark_20261104.json`](reference_dense_output_benchmark_20261104.json)
and [`benchmark_reference_dense_output.py`](../scripts/benchmark_reference_dense_output.py).

The public author ARDC package `riccati==2.0.0` was also screened on one mode.
At `eps=epsh=1e-9`, its spectrum error was `3.585e-4 dex` in `0.145 s`, and
tightening the tolerance did not reduce the observed offset. The 260-mode run
did not finish within a 45-minute CPU budget; it was cancelled without a
partial spectrum or accuracy result. This is recorded as an incomplete resource
screen, not as a full-grid method rejection, in
[`reference_ardc_res4_timeout_20261104.json`](reference_ardc_res4_timeout_20261104.json).
