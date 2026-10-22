# Cobaya adapter

Status: current adapter contract; validation results remain artifact-bound
Date: 2026-10-22
Code version: see the run commit in each validation artifact

`stiffgwpy_fast.cobaya.stiffGW.stiffGW` is a Cobaya `Theory`. It exposes the
derived parameters `Delta_Neff_GW`, `Delta_Neff_total`, `log10hc_prim_fyr`,
`f_end`, and `Delta_Neff_GW_error`; it supplies `f`, `omGW_stiff`, `hubble`,
`kappa_s`, and `kappa_r` to likelihoods.

## Mode mapping

```yaml
accuracy_mode: fast      # the single formal production/MCMC path
```

The resolution order is strict:

```
accuracy_mode  ->  preset defaults  ->  explicit user overrides only
```

In `stiffGW.yaml`, `h`/`col_step`/`z_tail`/`freq_res` default to `0` (a
sentinel meaning "use the selected accuracy_mode"). This keeps adapter
defaults from overriding the selected preset. The current `fast` preset uses
`z_tail=5.0`; setting a non-zero value is an explicit override.

When `accuracy_mode` is omitted, the high-level API selects the single `fast`
path. Historical `production` and transition-refine names are deprecated
aliases that map to `fast`; validation scripts may still select internal
presets. Use `accuracy_mode: null` only when deliberately managing legacy
module settings. The adapter passes the resolved settings as an immutable,
per-call configuration, so mode selection does not change process-global
defaults.

For concurrent evaluations within one process, `workqueue` parallel calls are
serialized because that backend is not thread-safe. Before Numba initializes
the backend, calls are serialized conservatively. Explicit per-call configs on
initialized TBB/OpenMP backends can run concurrently; legacy calls without a
config remain serialized. See the [backend throughput comparison](benchmarks.md).


## `eval_freqs`

Set `eval_freqs: [log10(f1), ...]` or `eval_freqs: /path/to/file` to add native
solve nodes. This passes likelihood bins to the fast solver through
`SGWB_iter_fast(..., eval_freqs=...)` and avoids interpolation error at steep
spectral features. By default (`null`), the solver uses its own grid and the
likelihood interpolates over the returned spectrum; Layer C measured a maximum
per-bin dex interpolation error of `3.1e-4`.

## Telemetry

Each run exposes the following fields through `theory.engine_stats`:
`fast_evals`, `fast_failures`,
`fast_guard_rejections`, `fast_physical_rejections`, `lsoda_evals`,
`lsoda_fallbacks`, `reference_evals`, `escalations`, `fallback_fraction`,
`escalation_fraction`, `last_eval_status`, `eval_status_counts`
(`FAST`, `FAST_ESCALATED`, `FAST_ESCALATION_FAILED`, `FAST_FAILED`,
`FAST_GUARD_REJECTED`, `FAST_PHYSICAL_REJECTED`, `REFERENCE`, `LSODA`,
`LSODA_FALLBACK`, `LSODA_FALLBACK_FAILED`), and the estimated `|Delta logL|`.
The derived `Delta_Neff_GW_error` is `NaN` when the selected engine does not
provide an error estimate. `close()` logs the summary and warns above a 5%
fallback/escalation fraction.
The `engine_stats` fields `DN_gw_error` and `spectrum_error` use different
quantities: `DN_gw_error` is the combined relative estimate for integrated
`Delta N_eff`, while `spectrum_error` is the selected accuracy mode's calibrated
log-spectrum error anchor. `spectrum_error` is `None` when no accuracy mode is
selected; it is not a point-local spectrum measurement.

遥测状态会区分快速求解失败、物理拒绝、保护条件拒绝及回退/升级失败。引擎没有提供误差估计时，
派生参数 `Delta_Neff_GW_error` 为 `NaN`，不以零冒充已知误差。
`engine_stats` 中的 `DN_gw_error` 是积分 `Delta_Neff` 的组合相对误差估计；`spectrum_error` 是所选
精度档位的频谱对数误差校准值。未选择精度档位时 `spectrum_error` 为 `None`，它不是当前参数点的频谱实测误差。

## Engine options

Choose `engine: fast | lsoda | reference`. `fallback: True` retries with LSODA
after a numerical failure and tags the result `LSODA_FALLBACK`; it never retries
a deterministic `shared_Neff_guard` rejection. `auto_escalate` with
`likelihood_sigma`/`dlogl_tol`
escalates when the estimated `|Delta logL|` exceeds the budget.
It also escalates when the fast error budget marks its solver-path calibration
as unsupported. The retry resets the model to its original `DN_eff` before
running the validation preset or independent reference engine. A failed retry is
reported as `FAST_ESCALATION_FAILED`; stale fast error-budget fields are cleared.

启用 `auto_escalate` 时，除似然误差超限外，若当前快速路径没有对应的误差校准，也会升级求解。
升级前会将模型恢复到原始 `DN_eff`，再运行验证档或独立参考引擎。
重算失败时状态为 `FAST_ESCALATION_FAILED`，旧快速求解的误差字段会清除。

The serial adapter only requires the `cobaya` extra. Install the separate
`mpi` extra when the execution environment actually needs `mpi4py`.
