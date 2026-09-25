# Cobaya adapter

Status: current adapter contract; validation results remain artifact-bound
Date: 2026-09-24
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
(`FAST` / `FAST_ESCALATED` / `REFERENCE` / `LSODA` / `LSODA_FALLBACK`), and the
estimated `|Delta logL|`.  `close()` logs the summary and warns above a 5%
fallback/escalation fraction.

## Engine options

Choose `engine: fast | lsoda | reference`. `fallback: True` retries with LSODA
after a numerical failure and tags the result `LSODA_FALLBACK`; it never retries
a deterministic `shared_Neff_guard` rejection. `auto_escalate` with
`likelihood_sigma`/`dlogl_tol`
escalates when the estimated `|Delta logL|` exceeds the budget.

The serial adapter only requires the `cobaya` extra. Install the separate
`mpi` extra when the execution environment actually needs `mpi4py`.
