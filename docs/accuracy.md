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
- The two archived Cartesian tail10 maps around this root use the same
  `reference.py` bytes and bitwise-identical 260-mode frequency arrays. Their
  secant slope is `0.02113`, giving a first-order fixed-point shift estimate
  `-1.08e-10` from the second map's residual. This is an interval secant, not a
  local derivative or rigorous root enclosure. See
  [`res4_reference_map_sensitivity_20261104.json`](res4_reference_map_sensitivity_20261104.json)
  and reproduce with
  [`analyze_reference_map_sensitivity.py`](../scripts/analyze_reference_map_sensitivity.py).
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
The Reference hot path also evaluates the existing background `CubicSpline`
coefficients with a Numba Horner evaluator. It matched the four-spline SciPy
interpolants within `2.22e-16` absolute error on a 250,000-point screen. On two
hard modes, the solve was `1.61x–1.76x` faster; the four-mode archived-spectrum
check at `rtol=1e-11` had maximum error `1.93e-8 dex`. A no-JIT scalar path
avoids compiling for the measured cheap low-frequency single-mode case. The
updated Prüfer backend also passes a full 260-mode same-input tail10 check at
the independent Cartesian root: map difference `4.39e-7` relative and maximum
spectrum difference `5.25e-5 dex`. A fresh full-grid paired Cartesian rerun is
now available below. The exact local root derivative remains unmeasured; the
available sensitivity is an archived interval secant. Raw
data are in [`reference_spline_eval_benchmark_20261104.json`](reference_spline_eval_benchmark_20261104.json).
For cold requests, the first measured cheap low-frequency mode avoids JIT and
remains `17.6 ms`; a costly mode pays `0.3 s` once and then uses the compiled
interpolant. A four-mode 4-worker smoke matched serial outputs bitwise but was
slower (`34.79 s` vs `17.50 s`), so it establishes thread correctness only.
The fresh-process low/high mode records and worker smoke are
[`reference_spline_cold_low_20261104.json`](reference_spline_cold_low_20261104.json),
[`reference_spline_cold_high_20261104.json`](reference_spline_cold_high_20261104.json),
and [`reference_spline_workers_20261104.json`](reference_spline_workers_20261104.json).
After the Cartesian spline change, the four-mode Prüfer A/B is `3.137x`
aggregate with maximum spectrum error `5.22e-7 dex`; see
[`reference_prufer_production_ab_post_spline_20261104.json`](reference_prufer_production_ab_post_spline_20261104.json).
The full-grid run record is
[`reference_prufer_spline_res4_refroot_20261104.json`](reference_prufer_spline_res4_refroot_20261104.json).
An expanded 8-frequency paired A/B (two alternating pairs per mode across
indices `0,65,100,130,160,195,220,259`) measured `3.626x` aggregate speedup,
or `3.586x` including JIT warmup, with per-mode P50/P95 `2.778x/4.174x` and
maximum error `3.94e-6 dex`. This eight-frequency subset remains a sample; the
complete 260-mode paired timing is reported below. Raw data are in
[`reference_prufer_production_ab_8mode_20261104.json`](reference_prufer_production_ab_8mode_20261104.json).

The same Prüfer/Cartesian comparison was then run once on every one of the
260 archived frequencies, alternating solver order by mode, at identical
`DN_eff=0.13373681817755848` and `z_tail=10`, with Prüfer `rtol=1e-9` and
Cartesian `rtol=1e-11`. All 260 pairs completed. The direct paired spectrum
absolute-dex difference has median/P95/maximum
`1.575e-7 / 6.654e-6 / 5.246e-5 dex`. Total propagation time was
`574.13 s` Prüfer versus `1992.45 s` Cartesian (`3.470x` aggregate); including
the separately measured `0.313 s` spline JIT gives `3.468x`. Per-mode speedup
P50/P95 is `2.632x/4.332x`; single-mode runtime P50/P95 was `2.627/5.690 s`
and `7.180/23.707 s`, respectively. These are serial measurements on this
res4 point, not general parameter or thread-scaling guarantees.

Reconstructed from the paired per-mode outputs, the Prüfer integral map is
`0.1337367594174389`, while the current Cartesian map is `0.1337368185045801`.
The latter differs from the archived Cartesian map by `3.24e-9` relative. The
Prüfer fixed-point residual is `-4.394e-7` relative; current Cartesian is
`+2.45e-9`. The paired map difference is `4.42e-7` relative and gives an
archived-secant first-order root-shift estimate of `-6.04e-8` in `DN_eff`
(`-4.51e-7` relative). The slope is an interval secant, not a local derivative
or a rigorous root bracket. Both methods report quadrature/interpolation
estimates near `1.84e-20 / 1.82e-10`. The raw 260-mode rows and reconstructed
summary are in
[`reference_prufer_cartesian_paired_fullgrid_20261104.json`](reference_prufer_cartesian_paired_fullgrid_20261104.json)
and
[`reference_prufer_cartesian_paired_fullgrid_analysis_20261104.json`](reference_prufer_cartesian_paired_fullgrid_analysis_20261104.json).
Reproduce with [`benchmark_res4_prufer_candidate.py`](../scripts/benchmark_res4_prufer_candidate.py)
and [`analyze_reference_paired_fullgrid.py`](../scripts/analyze_reference_paired_fullgrid.py).
The run telemetry records one worker and the software/resource environment;
this run did not measure peak RSS or worker scaling.

### Standalone Newton–Kantorovich phase-transfer prototype

The scalar damped tensor equation is written as
`T'' + p T' + exp(2z) T = 0`, with `p=3-1.5*sigma`. Setting
`u=exp(0.5*integral(p dN))*T` gives
`u'' + Q u = 0`, where `Q=exp(2z)-p'/2-p^2/4`. On each smooth high-frequency
interval, the prototype solves the logarithmic-derivative Riccati equation
`r'+r^2+Q=0` by Chebyshev collocation and Newton steps initialized by the
Liouville–Green expansion, then transfers the state through a phase basis.
Unlike the earlier Kummer-after-fit prototype, the phase is obtained by
solving the collocated Riccati residual directly. The implementation remains
standalone and does not alter the Reference API or production path.

For `r=.01`, `cr=0`, `n_t=2`, `T_re=1e4 GeV`, `kappa10=.01`, `DN_re=10`, and
`DN_eff=.13373681817755848`, DOP853 propagated production `z'=1.5*sigma-1`
with `(T,T')` from `z=-12` to the configured phase start. The phase path splits
the background at `nu=.1` and `nu=100`, using the matching radiation, FD-spline,
and high-frequency branches. At order 24, the `log10(f/Hz)=-12.6723` mode
crosses `nu=.1` and reaches the `z=10` tail; the `-10.772906` mode reaches the
tail before `nu=.1`. Their errors against current Cartesian Reference
(`rtol=1e-11`) are `1.7461e-8` and `4.1274e-9 dex`. A no-tail mode at
`log10(f/Hz)=-14.6013`, with phase start `z=6`, crosses both background
boundaries and is propagated to today without a tail approximation; its error
is `1.1137e-7 dex`. Order-24 dense Riccati residual maxima for the three cases
are `1.32e-14`, `8.68e-15`, and `2.66e-14`.

In fresh single runs, total candidate times including the DOP853 prefix were
`0.138`, `0.137`, and `0.499 s`, against `10.935`, `6.457`, and `9.520 s` for
Cartesian Reference. These timing ratios (`79.1x`, `47.3x`, `19.1x`) are
single-run prototype observations, not P50/P95 benchmarks; spline JIT warmup
(`~0.3 s`) is recorded separately. The no-tail case requires the later phase
start to pass the high-frequency screen. These are still only three modes at
one parameter point: no full-grid integral map, self-consistency, or
parameter-domain check has been completed. Reproduce with
[`prototype_nonosc_phase_res4_mode.py`](../scripts/prototype_nonosc_phase_res4_mode.py)
using `--frequency -12.6723`, `--frequency -10.772906114757099`, and
`--frequency -14.6013 --phase-start 6`; results are
[`nonosc_phase_res4_mode_20261104.json`](nonosc_phase_res4_mode_20261104.json),
[`nonosc_phase_res4_highfreq_mode_20261104.json`](nonosc_phase_res4_highfreq_mode_20261104.json),
and
[`nonosc_phase_res4_notail_mode_20261104.json`](nonosc_phase_res4_notail_mode_20261104.json).

A fixed-background screen was then run on the complete archived 260-frequency
res4 grid at the same `DN_eff` and tail depth. The Cartesian module digest was
checked against the paired archive before propagation. The phase path completed
213/260 modes (82.0%) with orders 16/24/32; 45 modes were rejected by the
high-frequency interval screen and two had non-positive normal-form `Q`. The
batch driver begins at `z_phase=5`, moves later after a screen/residual failure
and earlier when the handoff event is not reached, in `0.5` increments within
`[0,10]`; it records every attempt.
The Newton iteration limit remains 20: when its update test is not met, the method
now accepts only if the existing scaled collocation (`1e-10`) and dense
residual (`1e-9`) gates both pass. This residual-based check covered 14 modes
and 22 segments; their maximum observed collocation/dense residuals were
`1.18e-13/1.68e-13`. On completed modes, order 32 has maximum/P95 error
`1.548e-4/5.562e-8 dex`.
Two modes near the tail threshold use different handoff semantics (phase to
today, Cartesian to tail) and contribute the two largest errors, both below
`1.55e-4 dex`. For the 213 completed modes, the single-run order-32 phase time
has P50/P95 `0.172/0.282 s`; summed phase and archived Cartesian times are
`39.16/1990.97 s` (`50.8x` on this supported subset). This is not an all-grid
speedup claim or an interleaved benchmark; peak RSS was not measured. Failed
modes are preserved as failures, not included in the accuracy statistics. The
complete per-mode record, including phase-start attempts, residuals, branch
segments, timings, and failure messages, is
[`nonosc_phase_res4_fullgrid_20261104.json`](nonosc_phase_res4_fullgrid_20261104.json),
reproduced with
[`benchmark_nonosc_phase_res4_fullgrid.py`](../scripts/benchmark_nonosc_phase_res4_fullgrid.py).
This remains one fixed parameter point and does not validate an integral map,
self-consistency root, or parameter domain.

As a separate screen-sensitivity experiment, the same full grid was rerun with
`--minimum-gamma 0`, while retaining the unchanged Newton step limit, scaled
collocation and dense residual gates, and Cartesian comparison. It completed
219/260 modes (six more than the default screen); 39 still failed the Riccati
residual gate and two had non-positive `Q`. The maximum/P95 error remained
`1.548e-4/5.525e-8 dex`. The strict three-order screen recovers six additional
modes on this one parameter point and leaves 41 failures, so the default `gamma_min=10` screen is
retained pending broader validation. The separate raw record is
[`nonosc_phase_res4_fullgrid_gamma0_20261104.json`](nonosc_phase_res4_fullgrid_gamma0_20261104.json),
reproduced with `--minimum-gamma 0`.

The 39 remaining residual failures were then replayed independently at orders
16, 24, and 32 with the same recorded phase start. Order 16 completed 14,
order 24 completed 7, and order 32 completed none. Five no-tail modes
(`indices 212, 218, 224, 237, 244`) pass both 16 and 24; their maximum
order-pair difference is `5.47e-11 dex`, and maximum Cartesian error is
`3.28e-9 dex`. This raises the set of checked modes to 224/260 if these five
two-order cases are counted separately from the 219 strict three-order
results. The other 34 residual cases and two non-positive-`Q` modes remain
unresolved. Full per-order failures, residuals, and reference errors are in
[`nonosc_phase_res4_residual_order_probe_20261104.json`](nonosc_phase_res4_residual_order_probe_20261104.json),
reproduced with
[`analyze_phase_residual_order_failures.py`](../scripts/analyze_phase_residual_order_failures.py).

A second fixed-background cross-check at `DN_eff=0.133732083004284` (same
cosmological parameters, different background normalization) also passed at
three modes: the two tail modes have errors `8.42e-9` and `4.12e-9 dex`, and
the no-tail mode has error `3.39e-9 dex`. Candidate/Cartesian single-run times
were `0.134/11.573`, `0.131/6.578`, and `0.493/9.945 s`; these are not repeated
timing distributions. This supports stability across two nearby `DN_eff`
backgrounds only; it is not another cosmological parameter point or a
self-consistency check. Raw results are
[`nonosc_phase_res4_dn133732_tail_cross_20261104.json`](nonosc_phase_res4_dn133732_tail_cross_20261104.json),
[`nonosc_phase_res4_dn133732_high_20261104.json`](nonosc_phase_res4_dn133732_high_20261104.json),
and
[`nonosc_phase_res4_dn133732_notail_20261104.json`](nonosc_phase_res4_dn133732_notail_20261104.json).

The prototype CLI now records all base cosmological inputs explicitly. A
separate pressure point from the minimized Fast failure set uses the default
`Omega_bh2`, `Omega_ch2`, `H0`, and `A_s`, with `r=0.19649188544593907`,
`n_t=1.0622501865308032`, `T_re=1761.314039059193 GeV`,
`kappa10=0.011598483151630803`, `cr=0`, `DN_re=10`, and fixed `DN_eff=0`.
The two tail modes and one no-tail mode match current Cartesian within
`9.37e-9`, `2.56e-9`, and `3.86e-8 dex`; three additional middle frequencies
`-11.722586`, `-13.621946`, and `-14.0` are tail modes and match within
`2.42e-8`, `1.24e-8`, and `7.25e-9 dex`. Dense Riccati residuals stay below
`2.72e-13`. These are six fixed-background mode checks, not a complete
spectrum or self-consistency validation of the difficult point. Reproduce
with the explicit `--r`, `--cr`, `--n-t`, `--t-re`, `--kappa10`, `--dn-re`,
and `--dn-eff` options; raw records are
[`nonosc_phase_sample2_tail_cross_20261104.json`](nonosc_phase_sample2_tail_cross_20261104.json),
[`nonosc_phase_sample2_high_20261104.json`](nonosc_phase_sample2_high_20261104.json),
and
[`nonosc_phase_sample2_notail_20261104.json`](nonosc_phase_sample2_notail_20261104.json),
[`nonosc_phase_sample2_mid1_20261104.json`](nonosc_phase_sample2_mid1_20261104.json),
[`nonosc_phase_sample2_mid2_20261104.json`](nonosc_phase_sample2_mid2_20261104.json),
and
[`nonosc_phase_sample2_mid3_20261104.json`](nonosc_phase_sample2_mid3_20261104.json).

### Sample2 full native-grid check at the Reference root

The sample2 cosmology was then compared on its full `freq_res=8` native grid at
the previously converged `z_tail=5` Reference root anchor
`DN_eff=2.291411867022825`. A fresh one-pass Cartesian map on the same 800
frequencies returned `DN_gw=2.29141186744725`, for relative fixed-point
residual `-1.85e-10`. This is a fixed-background map at the archived root, not
a new outer root solve.

Before the residual-based Newton stop fix, the `gamma_min=10` phase screen
completed 700/800 modes; `gamma_min=0` completed 709/800. On the latter strict
16/24/32 set, maximum and P95 error were `8.668e-7` and `2.413e-7 dex`. The
phase path took about `145 s` against `45 s` for Cartesian. No tail-semantic
mismatches remain after the handoff fix. The original 91 unsupported modes
comprised 88 residual-gate failures and three non-positive-`Q` failures. The
candidate grid was incomplete, so no phase integral map or candidate-root
residual was claimed.

An independent order replay of the 88 residual failures found 14 no-tail modes
that pass two adjacent orders and Cartesian checks: 13 at 16/24 and one at
24/32. The maximum order-pair difference is `2.54e-10 dex`; maximum Cartesian
error is `2.72e-8 dex`. Count these separately from the 709 strict three-order
modes: 723/800 modes have at least one validated tier, with 77 still
unresolved.

The same grid exposed a tail-event bug: four modes crossed `z_tail=5` and later
had `z_today<5`, while the phase code inferred tail usage only from today's
`z`. It now uses the tail state at the phase endpoint when `z_phase=z_tail`;
the regression in [`test_nonosc_phase_prototype.py`](../tests/test_nonosc_phase_prototype.py)
failed before and passes after the change. Root-grid reference arrays and
phase results are
[`sample2_reference_root_res8_tail5_20261104.json`](sample2_reference_root_res8_tail5_20261104.json),
[`sample2_reference_root_res8_tail5_20261104.npz`](sample2_reference_root_res8_tail5_20261104.npz),
[`sample2_phase_reference_root_res8_tail5_20261104.json`](sample2_phase_reference_root_res8_tail5_20261104.json),
[`sample2_phase_reference_root_res8_tail5_gamma0_20261105.json`](sample2_phase_reference_root_res8_tail5_gamma0_20261105.json),
and
[`sample2_phase_root_residual_order_probe_20261105.json`](sample2_phase_root_residual_order_probe_20261105.json).
Reproduce with
[`benchmark_phase_sample2_reference_root_grid.py`](../scripts/benchmark_phase_sample2_reference_root_grid.py)
and
[`analyze_sample2_root_phase_order_failures.py`](../scripts/analyze_sample2_root_phase_order_failures.py).

The positive consistency-relation branch was also sampled with `cr=0.1` and
`r=.01`, for which the model derives `n_t=-r/8=-.00125` (the supplied `n_t=2`
does not control the spectrum on this branch). At the same res4 background
normalization, a tail mode and a no-tail mode match Cartesian within
`3.54e-8` and `3.21e-9 dex`. This checks branch semantics at two frequencies,
not a full grid or root; see
[`nonosc_phase_cr_positive_tail_20261104.json`](nonosc_phase_cr_positive_tail_20261104.json)
and
[`nonosc_phase_cr_positive_notail_20261104.json`](nonosc_phase_cr_positive_notail_20261104.json).

### Sample2 full native-grid check at the Reference root

The minimized sample2 cosmology was then checked on the full `freq_res=8`
native grid at the previously converged `z_tail=5` Reference root anchor
`DN_eff=2.291411867022825`. A fresh one-pass Cartesian map over the same 800
frequencies returned `DN_gw=2.29141186744725`, for relative fixed-point
residual `-1.85e-10`. This is a fixed-background map at the archived root, not
a new outer root solve.

At `gamma_min=10`, the phase method completed 700/800 modes with maximum/P95
error `8.668e-7/2.416e-7 dex`. At `gamma_min=0`, it completed 737/800, with
maximum/P95 `8.668e-7/2.400e-7 dex`; the phase path used `132.2 s` versus
`43.0 s` for Cartesian Reference, so it is slower on this grid. It did not
produce a complete candidate integral map because 63 modes remain unsupported
or fail the residual gates. The output keeps those failures explicit.

The Newton solve stops when the existing scaled collocation residual gate
(`1e-10`) is first met and still checks the existing dense residual gate
(`1e-9`). When a full Newton step does not lower the residual, it tries up to
eight halvings; if none lowers it, it tries a minimum-norm least-squares step
on the same linearized system and backtracks that step up to eight halvings.
If neither direction lowers the residual, it retains the original full step.
The 20-step cap and both gates remain unchanged. Residual stopping added five
strict modes to the prior gamma-zero count of709; conditional damping added
four, `[685, 712, 713, 751]`; the least-squares fallback added19 more
`[656, 671, 672, 687, 700, 701, 702, 703, 714, 716, 718, 726, 730, 733,
734, 755, 756, 767, 769]`. No prior strict mode was lost and maximum spectral
error is unchanged. The latest gamma-zero run took `132.2 s`; gamma-10 remains
700/800 and took `131.3 s`. Both remain slower than Cartesian. Focused modes
`656`, `672`, and `687` pass residual, dense, and spectrum checks on the
recorded Windows environment.

The Windows grid archive records mode687 as a strict pass, but Linux CI replay
of its order-16 regression returned collocation residual `7.76e-10`, above the
unchanged `1e-10` gate, and raised `ArithmeticError`. The candidate was safely
rejected rather than used. This is a platform-sensitive boundary in the
research prototype: the `737` strict / `786` total counts above are for the
recorded Windows run, not cross-platform certification. The test now requires
either a gate-compliant, reference-accurate candidate or this explicit
residual rejection for mode687; it does not relax any threshold.

A tail-event audit found and fixed a semantic bug exposed by these modes: four
frequencies crossed `z_tail=5` during propagation but later had `z_today<5`.
The prototype had inferred tail use from `z_today`, causing up to
`4.80e-3 dex` error against Reference. When the phase prefix ends exactly at
the tail threshold, it now uses the tail state at that endpoint. The focused
regression in [`test_nonosc_phase_prototype.py`](../tests/test_nonosc_phase_prototype.py)
fails before the fix and passes after it; the corrected 800-bin run has no tail
semantic mismatches and maximum spectral error below `1e-6 dex` on supported
modes.

Independent order replays of the 57 remaining gamma-zero residual failures find
20 modes with adjacent-order agreement and Cartesian validation. Their maximum
adjacent-order difference is `1.20e-10 dex`, and maximum Cartesian error is
`1.96e-7 dex`.

A follow-up independently validated all 70 completed order candidates (15 at
order 16, 30 at order 24, 25 at order 32) across 44 residual-failure
frequencies. Tight Cartesian (`rtol=1e-11`) and Prüfer (`rtol=1e-9`) references
agree within `1.19e-7 dex` maximum; candidate differences are at most
`8.15e-9` and `1.19e-7 dex`, respectively. All 20 adjacent-order modes are in
this dual-reference set. The unique validated count is therefore 781/800,
with 13 residual failures and six non-positive-`Q` modes unresolved. The
phase spectrum remains incomplete, so no phase `DN_gw` map or root is claimed.
Reproduce the dual-reference pass with
[`validate_sample2_residual_orders.py`](../scripts/validate_sample2_residual_orders.py)
and see [`sample2_phase_residual_orders_dual_reference_20261105.json`](sample2_phase_residual_orders_dual_reference_20261105.json).

A targeted phase-start scan then refined the internal handoff position for the
13 residual modes not yet checked. It tried offsets from `-0.5` through
`+0.4921875` in `0.0078125` increments around each recorded start, kept candidates
that passed the unchanged residual gates, tail semantics, and `1e-3 dex`
archived Cartesian check, then cross-checked the best candidate with tighter
Cartesian and Prüfer references. Six modes `[757, 759, 780, 782, 783, 786]`
pass; their maximum candidate errors are `3.51e-10 dex` vs tight Cartesian
and `8.54e-8 dex` vs Prüfer. Refining the scan from 1/64 to 1/128 spacing
found no additional modes. Unique validated coverage is now787/800, with
seven residual failures and six non-positive-`Q` modes unresolved. No complete phase
map or root is claimed. Reproduce with
[`probe_sample2_phase_start_offsets.py`](../scripts/probe_sample2_phase_start_offsets.py)
and see [`sample2_phase_start_offset_probe_20261105.json`](sample2_phase_start_offset_probe_20261105.json).

A lower-order screen of the same 57 remaining modes at orders 8/12/16 completed
0/0/15 per order and found no converged pair. It adds no new unique
frequency coverage. It is recorded as a negative coverage result and is not
used to alter the candidate or its residual gates; see
[`sample2_phase_root_lower_order_probe_20261105.json`](sample2_phase_root_lower_order_probe_20261105.json)
and [`probe_sample2_lower_collocation_orders.py`](../scripts/probe_sample2_lower_collocation_orders.py).

The remaining seven residual-failure modes were also tested with each original
branch interval split into two and four propagation segments, using the current
Newton/backtracking/least-squares solve and unchanged residual gates. No order
produced a screened candidate in this `1.67 s` probe; segment splitting is not
retained as a recovery path. See
[`sample2_phase_residual_segment_split_probe_20261105.json`](sample2_phase_residual_segment_split_probe_20261105.json)
and [`probe_sample2_residual_segment_splits.py`](../scripts/probe_sample2_residual_segment_splits.py).

Raw Reference and phase arrays, exact 800-bin rows, and provenance hashes are
in [`sample2_reference_root_res8_tail5_20261104.json`](sample2_reference_root_res8_tail5_20261104.json),
[`sample2_reference_root_res8_tail5_20261104.npz`](sample2_reference_root_res8_tail5_20261104.npz),
[`sample2_phase_reference_root_res8_tail5_20261104.json`](sample2_phase_reference_root_res8_tail5_20261104.json),
and
[`sample2_phase_reference_root_res8_tail5_gamma0_20261105.json`](sample2_phase_reference_root_res8_tail5_gamma0_20261105.json).
Its isolated Reference inputs are
[`sample2_reference_root_res8_tail5_gamma0_20261105.json`](sample2_reference_root_res8_tail5_gamma0_20261105.json)
and [`sample2_reference_root_res8_tail5_gamma0_20261105.npz`](sample2_reference_root_res8_tail5_gamma0_20261105.npz).
The same-code high-screen run is
[`sample2_phase_reference_root_res8_tail5_gamma10_20261105.json`](sample2_phase_reference_root_res8_tail5_gamma10_20261105.json);
its isolated Reference JSON/NPZ inputs are archived as
[`sample2_reference_root_res8_tail5_gamma10_20261105.json`](sample2_reference_root_res8_tail5_gamma10_20261105.json)
and [`sample2_reference_root_res8_tail5_gamma10_20261105.npz`](sample2_reference_root_res8_tail5_gamma10_20261105.npz).
Reproduce with
[`benchmark_phase_sample2_reference_root_grid.py`](../scripts/benchmark_phase_sample2_reference_root_grid.py)
and the residual probe
[`analyze_sample2_root_phase_order_failures.py`](../scripts/analyze_sample2_root_phase_order_failures.py).
The multi-order dual-reference validation record is
[`sample2_phase_residual_orders_dual_reference_20261105.json`](sample2_phase_residual_orders_dual_reference_20261105.json).

Before the scalar-spline optimization, rechecking Prüfer against the
dense-output-optimized Cartesian baseline gave
`3.349x` aggregate speedup on the same four modes, with per-mode P50/P95
ratios `3.085x/3.990x` and maximum spectrum error `1.401e-5 dex`; the endpoint
mode is only `1.736x` faster. The updated paired samples are in
[`reference_prufer_production_ab_postdense_20261104.json`](reference_prufer_production_ab_postdense_20261104.json).

The public author ARDC package `riccati==2.0.0` was also screened on one mode.
At `eps=epsh=1e-9`, its spectrum error was `3.585e-4 dex` in `0.145 s`, and
tightening the tolerance did not reduce the observed offset. The 260-mode run
did not finish within a 45-minute CPU budget; it was cancelled without a
partial spectrum or accuracy result. This is recorded as an incomplete resource
screen, not as a full-grid method rejection, in
[`reference_ardc_res4_timeout_20261104.json`](reference_ardc_res4_timeout_20261104.json).

## 2026-11-05 Reference evaluation cost reduction

The continuous-background Cartesian Reference now has a Numba path for the
existing piecewise `H2`/`sigma` equations. It evaluates the same SciPy cubic
coefficients and preserves the original Python evaluator as the cold path. A
batch of at least 32 frequencies that includes a costly mode warms the compiled
kernel; a one-frequency request does not pay this additional compilation.
Exact-zero tensor amplitude also skips this warmup.

On sample2's exact `freq_res=8`, 800-node Reference-root background
(`DN_eff=2.291411867022825`, `z_tail=5`, `rtol=1e-9`, Cartesian DOP853, one
worker), three alternating full-grid pairs gave Python P50/P95 `41.392/41.549 s`
and compiled P50/P95 `27.513/27.652 s`. The compiled/Python median ratio is
`0.6647` (33.53% lower warm runtime); first-call compile/warm cost was `0.579 s`
and was measured separately. Every paired `DN_gw`, spectrum node, frequency,
and tail mask is bitwise identical, including the archived Reference output;
maximum spectrum difference is zero. This is one difficult fixed background,
not a universal Reference speed claim. See
[`reference_h2_numba_fullgrid_ab_20261105.json`](reference_h2_numba_fullgrid_ab_20261105.json)
and reproduce with
[`benchmark_reference_h2_numba_fullgrid.py`](../scripts/benchmark_reference_h2_numba_fullgrid.py).

The `32`-frequency warmup threshold was measured on stratified sample2 batches.
A cold 16-node batch left the H2 kernel uncompiled; a cold 32-node batch compiled
it in `0.280 s`. Including that compile, the first 32-node candidate run took
`1.347 s`, versus the warm Python P50 of `1.619 s` (0.272s lower); the warm
compiled P50 was `1.072 s`. At 64 nodes, warm P50 was `2.188 s` versus `3.279 s`
Python. All batch outputs match bitwise. This supports the current dispatch
threshold on this hard background, not every workload. See
[`reference_h2_numba_batch_crossover_20261105.json`](reference_h2_numba_batch_crossover_20261105.json)
and [`benchmark_reference_h2_batch_crossover.py`](../scripts/benchmark_reference_h2_batch_crossover.py).

The compiled background matches the Python evaluator bitwise on 2,017 sample2
points, including nextafter probes at the reheating, thermal-table and neutrino
branch boundaries. The targeted Reference suite passes all44 non-slow tests;
the baseline profiler measured 4.32 million DOP853 RHS calls and 19.55s
inclusive background time in a 45.16 s 800-mode propagation. The compiled profile
measured 5.71 s in `_H2_and_sigma` and 30.63 s total instrumented propagation. Profile
and scalar probe records are
[`reference_fullgrid_profile_20261105.json`](reference_fullgrid_profile_20261105.json)
and [`reference_fullgrid_profile_compiled_20261105.json`](reference_fullgrid_profile_compiled_20261105.json),
plus [`reference_h2_numba_probe_20261105.json`](reference_h2_numba_probe_20261105.json).

The tensor RHS was then fused with the compiled background evaluator. The batch
path warms this kernel only for at least32 frequencies when a costly mode is
present; a cold one-frequency call keeps the original Python RHS. On the same
sample2 800-node map, three alternating pairs match the archive and each other
bitwise for `DN_gw`, every spectrum value, the grid and tail mask. Warm P50/P95
falls from `27.258/27.361 s` to `23.851/23.959 s` (12.5% lower P50); combined
JIT/background warm cost is `0.679 s`, measured separately. The 32-node batch
screen also repaid the extra RHS compilation on its first batch (`1.036 s`
including the measured `0.106 s` compile versus current P50 `1.098 s`). This
remains one cosmology and one worker; cold process startup, multi-worker memory
and other full-grid backgrounds are not certified. Reproduce with
[`benchmark_reference_numba_rhs_fullgrid.py`](../scripts/benchmark_reference_numba_rhs_fullgrid.py)
and [`benchmark_reference_numba_rhs_batch_crossover.py`](../scripts/benchmark_reference_numba_rhs_batch_crossover.py);
raw results are
[`reference_numba_rhs_fullgrid_ab_20261105.json`](reference_numba_rhs_fullgrid_ab_20261105.json)
and [`reference_numba_rhs_batch_crossover_20261105.json`](reference_numba_rhs_batch_crossover_20261105.json).

The fused RHS matches the original derivative tuple bitwise on five parameter
backgrounds at piecewise/transition boundaries and seeded states, and preserves
overflow and singular-background handling through the original Python fallback.
The complete Reference focused suite passes 56 tests with one slow-marked test
deselected; that slow test also passes separately. This is scoped regression
evidence, not a claim of global parameter-domain correctness.

The same-day `_find_start_N` change reuses the fixed horizon normalization
anchors and stops only after evaluating the first floating-point-stagnant
midpoint. All 800 sample2 start coordinates and five complete representative
mode outputs match the old code bitwise. After the H2/sigma optimization,
three paired full maps show P50/P95 `27.853/27.895 s` before and
`27.608/27.638 s` after (0.88% lower P50 runtime). See
[`reference_start_search_fullgrid_ab_20261105.json`](reference_start_search_fullgrid_ab_20261105.json)
and [`benchmark_reference_start_search_fullgrid.py`](../scripts/benchmark_reference_start_search_fullgrid.py).
The 800-coordinate count and search-call reduction are reproduced with
[`benchmark_reference_start_search.py`](../scripts/benchmark_reference_start_search.py).

A second-background check used `cr=1`, `r=0.01`, `T_re=2000 GeV`,
`kappa10=0.01`, `DN_eff=0.133732083004284`, `rtol=1e-11`, and two frequencies
(`-14.6013`, `-12.6723`). Python and compiled `H2`/`sigma` paths returned
bitwise-identical Cartesian and Prüfer outputs, event states, and `nfev` at
both frequencies. Cartesian versus Prüfer `Omega_GW` differed by at most
`3.80e-8` relative (`1.65e-8 dex`); this is a two-mode cross-formulation check,
not a full-grid or root certification. Single paired timings ranged from
`1.29x` to `1.45x` faster with the compiled evaluator, but these observations
are not a timing distribution or a performance claim; the separate first
compile/warm cost was `0.598 s`. See
[`reference_h2_numba_cr_positive_20261105.json`](reference_h2_numba_cr_positive_20261105.json)
and reproduce with
[`validate_reference_h2_numba_cr_positive.py`](../scripts/validate_reference_h2_numba_cr_positive.py).

The Reference focused suite now reports `56 passed, 1 deselected`, including
pointwise branch/boundary checks across five backgrounds (`cr=0/1`, low/high
reheating temperatures, and strong coupling). This broadens evaluator
equivalence evidence but does not replace multi-background full-spectrum
validation.

## 2026-11-05 256-point Fast domain pilots

Fresh 256-point Sobol and 256-point Latin-hypercube pilots sampled the broad
signed/log parameter mapping with alternating `cr` branches, seed `20261105`,
two threads, and the existing `N_inf<=120` resource-only preflight. Sobol
statuses were 6 success, 4 returned nonfinite, 14 solver failed, 60 exception,
78 resource skip and 94 resource preflight reject. Latin statuses were 7
success, 20 solver failed, 66 exception, 75 resource skip and 88 resource
preflight reject. None of the 13 successful outputs has an independent
accuracy certification; the four nonfinite results are not successes. These
are deliberately broad pilots, not success-rate or parameter-boundary claims.
Raw rows and run summaries are linked from
[`fast_domain_pressure_20261030.md`](fast_domain_pressure_20261030.md#2026-11-05-256-point-sobol-and-latin-hypercube-pilots).

## 2026-11-05 subnormal Reference integral and extreme-grid convergence

The `sobol_20261105_0146` extreme sample had a Reference spectrum spanning
subnormal floating-point values. SciPy PCHIP's harmonic-mean slope calculation
emitted overflow warnings while evaluating reciprocal secant slopes. In
`integrate_spectrum`, the spectrum is now divided by its peak only when a
nonzero subnormal value is present and the peak is below one; the integral and
both error estimates are restored by the same positive scale. A regression
test makes RuntimeWarnings errors and checks scale equivalence. The focused
Reference suite passes `57 passed, 1 deselected`.

At this point's native 77-node map, Cartesian and Prüfer tail5 agree within
`2.89e-7 dex` per mode and `3.28e-10` in `DN_gw`, so propagation backend error
is much smaller than the bolometric grid error. Default Fast versus the
independent tail10 modes has maximum spectrum difference `0.00645 dex` and
`DN_gw` difference `1.47%`, outside the requested gates. The `sigma_exact` plus
tail10 diagnostic reduces those same-grid differences to `0.000346 dex` and
`5.09e-5` relative, but does not establish a resolved bolometric integral or
justify a production change. A critical tail mode's tail10-to-tail12 change is
`1.60e-5 dex`; the separate mode/tail files document that convergence check.
These are fixed-map diagnostics for one deliberately extreme point.

On the exact fixed background, independent Prüfer Reference maps used the same
`rtol=1e-9`, `z_tail=5`, one worker, and frequency resolutions 1/2/4/8
(295/590/1181/2362 nodes). Their `DN_gw` values are
`6.174897172552e-106`, `6.114755710207e-106`, `6.109445931037e-106`, and
`6.109049829317e-106`. The res=4 to res=8 change is `6.48e-5` relative, below
the `2e-4` gate; res=1 to res=8 changes by `1.066e-2`, demonstrating material
frequency-grid under-resolution at the coarser setting. The native 77-node
Fast map's `DN_gw=1.127720207096e-105` is not a certified integral for this
extreme point. After scaling, the res=1 integral differs from the archived
pre-change result by about `5e-16` relative, and res=8 reproduces the
pre-change result within `2e-15`; no PCHIP warnings were emitted in the final
res=1/2/4/8 run. This is a single fixed-background Reference convergence check,
not a self-consistent root or an expansion of the certified Fast domain.
The archived native 77-mode Prüfer tail5/tail10 solutions were also reintegrated
without rerunning the ODEs: both `DN_gw` values changed by only `1.74e-16`
relative and emitted zero warnings. This isolates the integration fix from
mode propagation; see the [re-integration record](qmc_sobol146_post_scale_reintegration_20261105.json)
and [driver](../scripts/reintegrate_qmc_sobol146_post_scale.py).

## 2026-11-05 stable Fast PCHIP slopes and adaptive-grid point check

The Fast vectorized PCHIP slope routine had the same reciprocal-slope
overflow for subnormal spectral tails. It now scales each valid same-sign
secant pair by its smaller magnitude before forming the weighted harmonic
mean. A warning-as-error regression fails on the former reciprocal formula
and passes with the scaled expression. On a normal Fast point, 30 alternating
warm pairs are bitwise identical for `DN_gw` and every spectrum node; P50/P95
is `0.004932/0.005355 s` before and `0.004823/0.005249 s` after. This sample's
Fast run also emitted zero PCHIP warnings after the change.

The adaptive frequency grid is now selectable through the backward-compatible
optional `SGWB_iter(freq_grid=..., freq_grid_target=...,
freq_grid_max_points=...)` arguments; leaving them unset preserves the current
goal-grid default. For sample146 the default 77-node goal grid's `DN_gw`
relative error against the res=8 reference is `84.6%` (`45.8%` using the
larger absolute value as denominator). Existing adaptive at target `3e-4` uses
295 nodes but still differs by `2.58%` relative to Reference; tightening to `1e-6`, with
`sigma_exact=True` and `z_tail=10`, uses 321 nodes and differs by `7.07e-5`.
An independent Prüfer tail10 map at that exact 321-node array (`DN_eff` fixed)
gives `DN_gw` relative error `5.06e-5`, maximum spectrum error `3.72e-4 dex`
and P95 `4.90e-5 dex`, meeting the pointwise and integral gates. The Reference
run took `660.6 s`, used one worker and emitted no warnings. This validates one
extreme point through the explicit adaptive path only; the general error-budget
label remains `uncertified-path-configuration`, and default Fast settings are
unchanged.

Reproduce the threshold sweep with
[`probe_qmc_sobol146_fast_adaptive_grid.py`](../scripts/probe_qmc_sobol146_fast_adaptive_grid.py)
and the independent exact-map check with
[`validate_qmc_sobol146_fast_adaptive_reference.py`](../scripts/validate_qmc_sobol146_fast_adaptive_reference.py).
Raw records are
[`qmc_sobol146_fast_adaptive_grid_20261105.json`](qmc_sobol146_fast_adaptive_grid_20261105.json),
[`qmc_sobol146_fast_adaptive_reference_20261105.json`](qmc_sobol146_fast_adaptive_reference_20261105.json),
and [`fast_pchip_stable_slopes_ab_20261105.json`](fast_pchip_stable_slopes_ab_20261105.json).
After exposing the same low-level settings through `SGWB_iter`, the high-level
candidate reproduced the exact 321-node frequency hash and bitwise-identical
`DN_gw`; the core solver hash was unchanged. The prior Reference result did not
store a separate spectrum copy for that wrapper-only continuity check. See
[`qmc_sobol146_adaptive_api_continuity_20261105.json`](qmc_sobol146_adaptive_api_continuity_20261105.json)
and [`verify_qmc_sobol146_adaptive_api_continuity.py`](../scripts/verify_qmc_sobol146_adaptive_api_continuity.py).

The full grids, spectra, frequency/source hashes, warnings, runtime and
environment are recorded in
[`qmc_sobol146_reference_grid_convergence_20261105.json`](qmc_sobol146_reference_grid_convergence_20261105.json);
rerun with
[`validate_qmc_extreme_reference_grid_convergence.py`](../scripts/validate_qmc_extreme_reference_grid_convergence.py).
On a separate 800-node, normal-scale synthetic spectrum, 24 alternating
integrator pairs returned bitwise-identical `g2` and error estimates. P50 was
`0.079655 s` before and `0.079686 s` after scaling support (0.04% higher);
P95 was `0.080338/0.080249 s`. This measures only the integrator overhead, not
a cosmological Reference map; raw timings and inputs are in
[`reference_integrate_spectrum_scale_ab_20261105.json`](reference_integrate_spectrum_scale_ab_20261105.json)
and the [benchmark driver](../scripts/benchmark_reference_integrate_spectrum_scale.py).
The same-native-grid backend, Fast-mode, and tail-depth records are
[`qmc_sobol146_reference_crosscheck_20261105.json`](qmc_sobol146_reference_crosscheck_20261105.json),
[`qmc_sobol146_full_prufer_tail_convergence_20261105.json`](qmc_sobol146_full_prufer_tail_convergence_20261105.json),
[`qmc_sobol146_mode_tail_probe_20261105.json`](qmc_sobol146_mode_tail_probe_20261105.json),
and [`qmc_sobol146_prufer_tail_depths_20261105.json`](qmc_sobol146_prufer_tail_depths_20261105.json).

## 2026-11-05 adaptive Fast positive-cr fixed-background cross-check

On the standard density/Hubble background with `r=0.01`, input `n_t=2`,
`cr=1`, `T_re=2000 GeV`, `kappa10=0.01`, and `DN_eff_input=0`, the positive-cr
consistency branch derives `n_t=-0.00125`. The opt-in Fast adaptive grid at
target `1e-4`, `sigma_exact`, and `z_tail=10` returned 256 frequencies,
`frequency_grid_error=9.2055e-5`, and a self-consistent `DN_gw=0.00226202795`
in `0.736 s`. Its fixed-background spectrum was checked on the exact same
array by both Cartesian DOP853 and the independent Prüfer Reference (`rtol=1e-9`, one worker). Fast versus Cartesian and Prüfer have maximum positive-bin
errors `5.3041e-4` and `5.3030e-4 dex`; `DN_gw` map differences are
`8.86e-6` and `9.92e-6` relative. Cartesian versus Prüfer also agrees within
`2.14e-5 dex` on a stratified 32-bin subset. The full Fast/Prüfer frequencies
and raw positive spectra, exact parameters, hashes, and timing are in
[`adaptive_cr1_prufer_fullgrid_20261105.json`](adaptive_cr1_prufer_fullgrid_20261105.json);
reproduce with
[`validate_adaptive_cr1_prufer_fullgrid.py`](../scripts/validate_adaptive_cr1_prufer_fullgrid.py).

The same Prüfer Reference on a denser native `freq_res=2` grid (484 nodes)
returns `DN_gw=0.0022622765591`, a `9.9973e-5` relative change from the
256-node adaptive-grid Prüfer map. The Fast output differs from that dense map
by `1.0989e-4` relative. Both are within the `2e-4` DN gate and provide a
two-grid convergence screen for this fixed background. The raw 484-node
frequencies and spectrum, environment, and hashes are in
[`adaptive_cr1_prufer_res2_20261105.json`](adaptive_cr1_prufer_res2_20261105.json);
reproduce with
[`validate_adaptive_cr1_prufer_res2.py`](../scripts/validate_adaptive_cr1_prufer_res2.py).

This two-grid difference is not a formal quadrature-error bound or a
self-consistent Reference root certification. The higher Prüfer `freq_res=4`
map and an alternate Cartesian `freq_res=2` map each reached the 1200-second
process CPU budget without a saved result; these are resource timeouts, not
accuracy failures. The Fast-only target sweep is non-monotone: the `1e-4` and
`1e-6` outputs differ by `1.12e-4` relative in `DN_gw`, while the `3e-5` and
`1e-6` outputs differ by `2.94e-4`. Keep this target sweep separate from the
independent Reference result. The partial Cartesian run and resource outcome
are recorded in
[`adaptive_standard_background_validation_20261105.json`](adaptive_standard_background_validation_20261105.json).

On 32 stratified modes, Reference `workers=1` took `84.36 s` and `workers=4`
took `132.56 s`; outputs were bitwise equal. Four workers are not retained as
a speed path for this workload. Timings, the five-target Fast sensitivity
screen, and sampled-memory limitations are recorded in the full-grid JSON.

## 2026-11-05 Fast validation above the former `N_inf=120` pilot cutoff

The earlier Sobol/LHS harness skipped `N_inf>120` points as a resource policy;
this is not a Fast code limit. A standard-positive-density point with
`T_re=1e-3`, `DN_re=100`, `r=0.01`, `n_t=0.5`, `cr=0`, and `N_inf=122.502`
is directly evaluable. The default 75-node grid returns `DN_gw=6.4825083e-9`;
adaptive target `1e-4` uses 242 nodes and returns `6.4111361e-9` with an
internal frequency-grid estimate `7.3651e-5`. The two Fast outputs differ by
about `1.11%`, so neither finite output nor that estimator alone certifies the
bolometric value.

On eight stratified bins at `z_tail=5`, same-background Cartesian and Prüfer
References agree within `1.3e-7 dex`, while both differ from Fast by up to
`2.89e-3 dex`. `sigma_exact=True` leaves this discrepancy essentially
unchanged; transition refinement also does not remove it. Tail sensitivity
identifies the dominant error: the selected-bin maximum falls to
`1.46e-4 dex` at `z_tail=8` and `1.61e-5 dex` at `z_tail=10`.

The tail8, target `1e-4` candidate passed full 242-bin same-grid Prüfer
comparison (`1.587e-4 dex` maximum, `4.52e-5` relative `DN_gw` map error),
but failed the dense-grid integral check: adaptive 242 versus native res2 451
Prüfer maps differed by `3.84e-4` relative. Tightening the adaptive target to
`1e-5` produced 299 bins; full same-grid Prüfer comparison then gave maximum
spectral error `1.656e-4 dex` and `2.50e-5` relative `DN_gw` error. Native res2
Prüfer changed `DN_gw` by `8.95e-5` relative from the adaptive-grid map; Fast
versus this dense map differs by `1.145e-4` relative. These pass the stated
fixed-background spectrum/map screens for this one point. They do not certify
the self-consistent Reference root or imply a global automatic tail rule.

Fast runtime was `0.667 s` for the tail8/target `1e-5` candidate; Reference
Prüfer took `97.22 s` on the 299-node exact Fast grid and `175.84 s` on the
451-node native res2 grid. Raw frequencies, spectra, source hashes, tolerance,
tail depth, and resource data are in
[`fast_ninf122_tail8_target1e5_validation_20261105.json`](fast_ninf122_tail8_target1e5_validation_20261105.json);
reproduce with
[`validate_fast_ninf122_tail8_target1e5.py`](../scripts/validate_fast_ninf122_tail8_target1e5.py).
The six-configuration tail/sigma/transition screen is in
[`fast_ninf122_tail_sensitivity_20261105.json`](fast_ninf122_tail_sensitivity_20261105.json),
with driver
[`probe_fast_ninf122_tail_sensitivity.py`](../scripts/probe_fast_ninf122_tail_sensitivity.py).

## 2026-11-05 Fast self-consistent root screen above `N_inf=120`

The tail8/target-`1e-5` Fast result at `N_inf=122.502` gives
`DN_eff_final=6.4131520531e-9`. To check whether the fixed-background map
accuracy carries through the self-consistent equation
`x = DN_eff_input + DN_gw(x)`, two additional Prüfer maps were evaluated at
the native res2 grid: one at the Fast value and one at a `+1e-3` relative
offset. Their residuals bracket a native-grid secant root at
`6.4138866031e-9`; the root is `1.1454e-4` relative from Fast. The secant
residual at that candidate is `4.55e-8` relative after correction. On the
299-node adaptive Fast grid, a separate interval secant gives
`6.4133122911e-9`; the two grid-specific estimates differ by `8.9542e-5`
relative. The adaptive-grid residual changes sign between the Fast value and
the native-grid secant candidate over a `1.1454e-4` relative interval.

This meets the `2e-4` root-error screening threshold for this one point using
two frequency grids and two interval secants. It is not a rigorous interval
enclosure, a formal quadrature bound, or a full iterative Reference solve;
the secant slopes are finite-interval estimates. The artifact also records
that the native-grid residual does not change sign over the narrower
Fast-to-native-secant interval, so the two grids do not share one common
bracket. Do not generalize this result into a global root guarantee or an
automatic tail/grid dispatch rule. The exact maps, residuals, grid hashes,
resource metadata, and derivation are in
[`fast_ninf122_reference_root_secant_20261105.json`](fast_ninf122_reference_root_secant_20261105.json),
reproduced by
[`probe_fast_ninf122_reference_root_secant.py`](../scripts/probe_fast_ninf122_reference_root_secant.py).

## 2026-11-05 Fast full-grid screen at `N_inf=222.502`

Increasing `DN_re` from `100` to `200` at the same other inputs gives a
derived `N_inf=222.502`. The tail8/adaptive-target-`1e-5` Fast solve completed
with finite positive output on the same 299-node frequency array as the
`N_inf=122.502` point. The full 299-node Prüfer comparison at the Fast fixed
background has maximum spectrum error `1.4962e-4 dex` (P95
`1.4563e-4 dex`) and same-grid `DN_gw` difference `2.50e-5` relative. Native
451-node `freq_res=2` Prüfer changes the map by `8.95e-5` relative; Fast differs
from this denser map by `1.1453e-4`. All four pointwise spectrum/map screens
pass. The Fast runtime is `0.623 s`; Prüfer took `90.39 s` on the exact Fast
grid and `160.15 s` on the native grid. On the same 299 modes, the candidate
Prüfer spectrum differs from the parent point by at most `2.03e-5 dex`,
consistent with an e-fold coordinate shift for this particular parameter
change.

This extends full fixed-background Fast validation to `N_inf=222.502`. It does
not certify the independent self-consistent root at the new point or prove
global invariance under `DN_re` shifts. The original 12-mode early screen and
the full Fast/Prüfer arrays, map convergence, source hashes, and resource data
are in [`fast_ninf222_selected_prufer_20261105.json`](fast_ninf222_selected_prufer_20261105.json)
and [`fast_ninf222_fullgrid_prufer_20261105.json`](fast_ninf222_fullgrid_prufer_20261105.json),
with drivers [`validate_fast_ninf222_selected_prufer.py`](../scripts/validate_fast_ninf222_selected_prufer.py)
and [`validate_fast_ninf222_fullgrid_prufer.py`](../scripts/validate_fast_ninf222_fullgrid_prufer.py).

## 2026-11-05 Fast low-temperature positive-cr cross-check

For the independent `cr=1` branch, the input `n_t=0.5` and `DN_re=200` are
not independent: the consistency branch derives `n_t=-0.00125` and overwrites
`DN_re` with `57.9462`. At `T_re=1e-3`, the resulting `N_inf=80.4482` Fast
solve completes on 339 adaptive nodes (`target=1e-5`, `z_tail=8`,
`sigma_exact=True`) in `0.568 s`; its local grid estimate is `9.621e-6` and
`DN_gw=3.10600045e-9`.

The full same-grid 339-node Prüfer map gives maximum spectrum error
`1.5079e-4 dex`, P95 `1.4550e-4 dex`, and `DN_gw` difference `3.03e-7`
relative. The first native `freq_res=2` grid (383 nodes) misses the integral
screen narrowly: adaptive-to-res2 difference `2.038e-4`, Fast-to-res2
`2.041e-4`. Refinement shows res2 is under-resolved: res2-to-res4 (765 nodes)
changes `DN_gw` by `2.764e-4`. The converged pair res4→res8 (1530 nodes)
differs by `8.971e-5`; Fast-to-res8 differs by `1.743e-5`, and the adaptive
map-to-res8 difference is `1.713e-5`. Thus the fixed-background spectrum and
map gates pass against the denser converged Prüfer pair even though the
intermediate res2 result does not.

This validates this fixed background only. It does not solve the independent
self-consistent Reference root, and the `cr>0` input coordinates must be
reported with their derived effective values. Full arrays, grid hashes,
resource data, staged drivers, and the initial res2 near-miss are retained in
[`fast_cr1_lowT_fullgrid_prufer_20261105.json`](fast_cr1_lowT_fullgrid_prufer_20261105.json),
[`fast_cr1_lowT_res4_prufer_20261105.json`](fast_cr1_lowT_res4_prufer_20261105.json),
and [`fast_cr1_lowT_res8_prufer_20261105.json`](fast_cr1_lowT_res8_prufer_20261105.json).

For the same positive-cr low-T point, a fixed-point sensitivity screen used
`x = DN_eff_input + DN_gw(x)` with a `+1e-3` relative probe. The adaptive-grid
interval secant estimates `x=3.1059995117e-9`; native res4 estimates
`x=3.1062249717e-9`. Their relative difference is `7.2583e-5`, and Fast
differs from the native res4 estimate by `7.2285e-5`. At the native-res4
secant candidate, an additional adaptive-grid map gives a slope-corrected
relative residual of `7.2582e-5`. The res4 map changes sign between sampled
endpoints; adaptive residuals remain positive at both endpoints, so the grids
do not share a sign-changing bracket. The independent res8 map at Fast `x0`
also has positive residual and is consistent with the adaptive screen.

This meets the `2e-4` two-grid secant screening gate at this point. Secant
slopes are finite-interval estimates, and there is no common grid bracket, so
this is not a rigorous root enclosure or a full iterative Reference solution.
Raw maps and calculations are in
[`fast_cr1_lowT_root_secant_20261105.json`](fast_cr1_lowT_root_secant_20261105.json),
reproduced by
[`probe_fast_cr1_lowT_root_secant.py`](../scripts/probe_fast_cr1_lowT_root_secant.py).

## 2026-11-05 High-blue Fast integration screen at `N_inf=222.502`

A nonphysical stress point with `cr=0`, `r=10`, `n_t=2`, `kappa10=0.01`,
`T_re=1e-3`, and `DN_re=200` is directly evaluable. Fast returns
`DN_gw=0.1580247` on 51 adaptive nodes at target `1e-5`; its same-grid Prüfer
map passes the spectrum gate (`1.4542e-4 dex` max), but the Fast map differs
from native res4 by `7.575e-4`. Tightening the target to `1e-6` gives 57 nodes
and reduces Fast-to-res4 to `2.329e-4`, still just outside the `2e-4` gate.

At target `1e-7`, Fast uses 73 nodes and returns `DN_gw=0.1579128211` in
`0.353 s`. The exact-array Prüfer map has maximum spectrum error
`1.4639e-4 dex` and `DN_gw` map difference `1.525e-5`; native res2→res4
changes the map by `1.712e-4`, and Fast differs from res4 by `6.416e-5`. All
four fixed-background spectrum/map gates pass. The local panel indicator's
maximum is `9.32e-8`, while the diagnostic sum of those panels is `1.10e-6`;
the sum is not a rigorous error bound. This point shows that the looser target
can underresolve the bolometric integral even when the exact-grid spectrum
check passes.

The result is specific to this extreme blue point and does not define an
automatic target rule or certify its self-consistent Reference root. The raw
target `1e-5`, `1e-6`, and `1e-7` maps and drivers are retained in
[`fast_high_blue_fullgrid_prufer_20261105.json`](fast_high_blue_fullgrid_prufer_20261105.json),
[`fast_high_blue_target1e6_prufer_20261105.json`](fast_high_blue_target1e6_prufer_20261105.json),
[`fast_high_blue_target1e7_prufer_20261105.json`](fast_high_blue_target1e7_prufer_20261105.json),
[`validate_fast_high_blue_fullgrid_prufer.py`](../scripts/validate_fast_high_blue_fullgrid_prufer.py),
[`validate_fast_high_blue_target1e6.py`](../scripts/validate_fast_high_blue_target1e6.py),
and [`validate_fast_high_blue_target1e7.py`](../scripts/validate_fast_high_blue_target1e7.py).

## 2026-11-05 High-blue Fast root sensitivity screen

At the target-`1e-7` Fast fixed point, `x0=0.1579128211`. A `+1e-3`
relative `DN_eff` probe gives adaptive-grid and native-res4 interval secant
root estimates `0.1579152800` and `0.1579024733`, respectively; their
relative difference is `8.110e-5`. Fast differs from the res4 estimate by
`6.553e-5`. A further adaptive-grid map at the res4 candidate gives a
slope-corrected relative residual `8.114e-5`. The adaptive interval brackets
the residual sign change; the res4 interval does not, so there is no shared
sign-changing bracket.

This passes the `2e-4` two-grid secant screen for this point. The secant
slopes are interval estimates, not rigorous local derivatives, and this is not
a full iterative Reference root solve. Raw maps, residuals, hashes, and
resource metadata are in
[`fast_high_blue_root_secant_20261105.json`](fast_high_blue_root_secant_20261105.json),
reproduced by
[`probe_fast_high_blue_root_secant.py`](../scripts/probe_fast_high_blue_root_secant.py).
