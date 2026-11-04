# Repository scripts

[中文说明](README_zh.md) | English

This directory contains validation, benchmark, profiling, packaging, and
historical diagnostic scripts. Run commands from the repository root so that
package imports and artifact paths resolve as documented.

## Script groups

| Prefix or file | Purpose |
|---|---|
| `validate_*.py` | Replay fast/reference comparisons, parameter edges, sweeps, and mode checks |
| `build_*.py` | Assemble validation matrices and manifests from existing artifacts |
| `bench_*.py`, `benchmark_*.py` | Measure runtime, accuracy, or a named diagnostic hypothesis |
| `profile_*.py` | Attribute solver-stage time or measure property/call overhead on a named path |
| `smoke_*.py`, `verify_distribution.py` | Check installed-wheel resources and package boundaries |
| `plot_*.py` | Render plots from saved validation or benchmark data |
| `_resource_budget.py`, `_smoke_probe.py` | Shared support for thread limits and isolated smoke checks |

`validate_fast.py`, `validate_random.py`, and `convergence_study.py` select
signal-region channels from the reference spectrum (`log10(Omega_GW) >= -30`).
Their linear relative-error helper uses the candidate spectrum as denominator.

`benchmark_phase_a.py` compares the plain fixed grid with exact reheating-kink
splitting for one built-in parameter case. Its spectrum-relative metrics use
the shared helper with the candidate spectrum as denominator. The default
reference grid uses all candidate frequencies; `--subset-reference` is a
historical sparse-grid diagnostic that can change the reference `DN_gw` value.
This script is a scoped comparison, not a parameter-space accuracy certification.

`check_random_freq.py` compares both engines at one fixed parameter point and
marks ten seeded frequency channels while reporting whole-grid residuals. Those
channels are a visualization sample, not ten independent parameter cases. Its
full-grid relative residual includes near-zero channels and should be read with
the dex residual; its PNG defaults to the `scripts/` directory unless an output
path is supplied.

`bench_fast.py` compares one LSODA runtime with cold and warm fast timings; it
is a runtime/regression benchmark, not an accuracy comparison. Timed calls
use the fast preset grid and quadrature, with preset tolerance and frequency
resolution unless `--tol` or `--freq-res` overrides them. By default timing
omits the exact reheating-kink split; `--kink-split` includes it. The separate
safety-wrapper diagnostics use the formal fast kink-split path and record their
configuration in JSON. Its P95 uses linear percentile interpolation, so the
default three repeats are only a rough tail estimate.

`benchmark_transition_frontier_ab.py` alternates formal fast and a matched
transition-refinement configuration on the two `DN_re` edges. It records raw
solver-call timings and per-variant output digests; it does not establish
accuracy by itself. Read the matched Oracle C records separately.

`benchmark_res4_scalar_step.py` runs paired full-root Fast comparisons between
the default h=.005 path and the aligned h=.02 transition-refinement path. It
uses the saved tail10 Reference spectrum and records fixed-frequency errors,
thread affinity, first-call and warm P50/P95 timings, post-call RSS, batch
throughput, repeat digests, and all 260 per-mode errors for the three variants.
`benchmark_res4_reference_root.py` runs one
independent fixed-background Reference map at an explicit `DN_eff` on the same
260 native frequencies; it can take tens of minutes and does not perform the
outer self-consistency iteration itself.
`benchmark_sample11_reference_map.py` repeats the exact-input map check for
the retained 94-mode sample11 case using its saved Fast-root frequency grid.
It records the map residual, quadrature/interpolation estimates, and same-DN
per-mode spectrum errors, and saves the raw Reference spectrum beside its JSON.
`benchmark_sample11_sigma_exact.py` alternates complete h=.02 transition-refined
sample11 roots using the continuous-background primitive by default and the
explicit node primitive. It records matched-grid errors against the saved
tail10 Reference spectrum and warm P50/P95 timings. Pass `--reference` and
`--reference-map` to select the same-input map at the continuous-primitive root.
`benchmark_background_primitive_guard.py` alternates six normal Fast regimes
with the pre-kernel primitive finite check enabled/disabled to quantify its
runtime overhead and confirm identical spectra.

`benchmark_lvk_likelihood_cache.py` requires the Cobaya extra and compares the
previous LVK likelihood calculation with the cached-data/NumPy-reduction path
on the packaged 54,593-bin table. It records raw per-call samples, likelihood
differences, and a lower-support regression probe; it does not measure full
Cobaya or MCMC runtime.
By default it writes a timestamped JSON file under `docs/`; `--json` selects an
explicit output path.

Many `benchmark_*_spike.py` files are historical, standalone candidate tests;
they do not select or modify a production solver path. Read a script's CLI
arguments (or run `--help` when available) and check its output path first.
Preserve dated artifacts and
use a new output name for a new experiment unless replacement is intentional.
The SageNet comparison driver also requires a separate SageNet checkout and
runtime; SageNet is not installed by this package's extras.

## Reproducible runs

Use the documented environment, thread budget, seeds, warmup policy, and repeat
count for the question being tested. Separate cold/JIT startup from warm timing.
Record the source commit and keep numerical accuracy comparisons separate from
speed comparisons. See [`docs/reproducibility.md`](../docs/reproducibility.md),
[`docs/benchmarks.md`](../docs/benchmarks.md), and
[`docs/experiment_catalog.md`](../docs/experiment_catalog.md).
