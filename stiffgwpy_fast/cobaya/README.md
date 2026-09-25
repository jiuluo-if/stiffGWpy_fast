# Cobaya adapter package

[中文说明](README_zh.md) | English

This package provides the Cobaya theory adapter and its packaged likelihoods.
The detailed option contract is documented in
[`docs/cobaya.md`](../../docs/cobaya.md).

## Theory entry point

Use the fully qualified class name in a Cobaya YAML file:

```yaml
theory:
  stiffgwpy_fast.cobaya.stiffGW.stiffGW:
    engine: fast
    fallback: True
    accuracy_mode: fast
```

The serial adapter requires the `cobaya` extra. Install the separate `mpi`
extra only when the run actually uses `mpi4py`.

## Numerical contract

- The adapter resolves `accuracy_mode`, preset defaults, then explicit user
  overrides. Sentinel zero values in the YAML do not mask the selected preset.
- `fast` is the single formal user profile. Numerical failures can use LSODA
  fallback only when enabled; deterministic physical guards are not retried.
- `eval_freqs` adds likelihood frequencies as native fast-solver nodes. Those
  nodes do not enter the bolometric integration support grid.
- `engine_stats` reports per-run engine, failure, guard, fallback, and
  escalation telemetry.

The package also contains the LVK and PTA likelihood adapters and their data
resources. See the [likelihood and resource index](likelihoods/README.md) for
the adapter-to-data mapping. Use the reproducibility and package-verification
instructions before publishing a wheel.
