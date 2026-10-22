# Accuracy and evidence boundaries

Status: current summary of scoped evidence; no full parameter-space release certification

Date: 2026-10-22

The formal fast profile is the single `goal-kink-hybrid` mode. The independent
continuous-sigma DOP853 pipeline is the precision reference; LSODA is a
regression/runtime path, not a truth oracle. Every numerical result below is
bounded by the artifact's commit, reference convention, parameter points, and
frequency grid.

## Latest scoped oracle evidence

The latest detailed accuracy audit is the dated
[`fast_v02_audit_report.md`](fast_v02_audit_report.md), generated from its
recorded 2026-09-12 evidence. It is not a full parameter-space certification or
a fresh precision rerun at the current repository HEAD.

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

## 2026-10-22 fresh `DN_re` edge probe

A current-HEAD full-grid comparison adds two reheating-parameter edges to the
three previously sampled Oracle C points. Oracle WKB and deep-reference
`DN_gw` agree within `2.4e-6` and `5.7e-6`, respectively. The formal fast
profile is within `1.2e-5` on the earlier three points, but has relative
differences of `2.74e-3` (`edge_dnre_lo`) and `2.50e-3` (`edge_dnre_hi`) on the
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

## Current full-path solver profile is separate

The 2026-10-21 full-path profile used solver source at commit
`cdbdb387479e11cd87a38d76a0e58348ca13e70f`. It covers default, high-T, stiff,
and high-kappa with 25 warm repeats per case; the first process call is
reported separately. The warm medians are `6.790`, `10.284`, `10.095`, and
`10.084 ms`, respectively. This stage profile did not run a matched-resource
LSODA comparison or a new oracle-accuracy sweep, so it adds no new precision
coverage and establishes no LSODA speed ratio. Full per-case values and
resource provenance are in [`benchmarks.md`](benchmarks.md).

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
