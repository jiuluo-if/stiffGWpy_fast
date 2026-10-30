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

The separate blue-tilt res4 `h=.02, transition_refine=True` fixed-background
Reference map is still running; do not use this endpoint correction's Oracle C
checks to claim that candidate's self-consistent root is verified.
