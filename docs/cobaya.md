# Cobaya adapter

Status: current
Date: 2026-09-03
Code version: see manifest `commit`

`stiffgwpy.cobaya.stiffGW.stiffGW` is a Cobaya `Theory`.  It exposes the derived
params `Delta_Neff_GW`, `Delta_Neff_total`, `log10hc_prim_fyr`, `f_end`,
`Delta_Neff_GW_error`, and provides `f`, `omGW_stiff`, `hubble`, `kappa_s`,
`kappa_r` to the likelihoods.

## Mode mapping

```yaml
accuracy_mode: fast      # plain-grid, speed-first
accuracy_mode: production # transition-refine, precision-first (default)
```

The resolution order is strict:

```
accuracy_mode  ->  preset defaults  ->  explicit user overrides only
```

In `stiffGW.yaml`, `h`/`col_step`/`z_tail`/`freq_res` default to `0` (a
sentinel meaning "use the selected accuracy_mode").  This prevents the Cobaya
default values from silently overriding the preset (the historical bug where a
`z_tail: 7.0` YAML default masked the preset's `z_tail: 8.0`).  Setting any of
them to a non-zero value is an explicit override.

The high-level fast solver also defaults to `production` when `accuracy_mode`
is omitted. `accuracy_mode: null` is reserved for compatibility with callers
that deliberately manage legacy module settings. The adapter passes resolved
settings as a per-call immutable configuration, so selecting a mode does not
mutate process-global solver defaults.

## `eval_freqs`

Set `eval_freqs: [log10(f1), ...]` or `eval_freqs: /path/to/file` to force-add
native solve nodes.  This is the path by which likelihood frequency bins reach
the fast solver as native nodes (`SGWB_iter_fast(..., eval_freqs=...)`),
removing interpolation error at steep spectral features.  By default (`null`)
the solver uses its own grid and the likelihood interpolates over the returned
spectrum (Layer C measured that interpolation per-bin dex error at ≤3.1e-4).

## Telemetry

`theory.engine_stats` exposes `fast_evals`, `fast_failures`,
`fast_guard_rejections`, `fast_physical_rejections`, `lsoda_evals`,
`lsoda_fallbacks`, `reference_evals`, `escalations`, `fallback_fraction`,
`escalation_fraction`, `last_eval_status`, `eval_status_counts`
(`FAST` / `FAST_ESCALATED` / `REFERENCE` / `LSODA` / `LSODA_FALLBACK`), and the
estimated `|Delta logL|`.  `close()` logs the summary and warns above a 5%
fallback/escalation fraction.

## Engine options

`engine: fast | lsoda | reference`.  `fallback: True` reruns with LSODA on a
numerical failure (tagged `LSODA_FALLBACK`); a deterministic `shared_Neff_guard`
rejection is never retried.  `auto_escalate` with `likelihood_sigma`/`dlogl_tol`
escalates when the estimated `|Delta logL|` exceeds the budget.

The serial adapter only requires the `cobaya` extra. Install the separate
`mpi` extra when the execution environment actually needs `mpi4py`.
