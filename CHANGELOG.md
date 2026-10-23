# Changelog

## Unreleased — 2026-10-23

- Recorded the LVK likelihood's six-case, 200-pair A/B on the packaged data;
  overlapping cases show `1.90x`–`9.75x` lower per-call latency, with maximum
  absolute log-likelihood change `1.21e-10`. The no-overlap case is scoped
  separately; these are not full MCMC timings.

- Cache fixed LVK frequency/data columns during likelihood initialization, use
  NumPy's float64 pairwise chi-square reduction, and leave unsupported bins at
  zero without querying the cubic interpolator out of range.

- Reject non-boolean `kink_split` values so strings such as `"false"` cannot
  silently activate the exact-grid path.

- Made benchmark CPU-affinity setup best-effort and kept resource telemetry
  passive before Numba execution. Failed pinning is no longer reported as a
  successful selected CPU set, an empty CPU set is never sent to psutil's
  reset-all operation, and requested versus active thread settings are recorded
  separately.

- Refreshed both experiment catalogues with the 2026-10-23 transition-refinement
  frontier, its mixed paired timing, and direct links to the matched Oracle and
  timing records.

- Corrected `bench_fast.py` to forward the selected fast preset's convergence
  tolerance and frequency resolution, record its solver configuration, and
  compute P95 with NumPy's linear percentile method. The script now describes
  LSODA as a runtime baseline and distinguishes the unsplit grid from the
  exact-kink path. Existing result files remain historical and are not
  retroactively recalculated.

- Added a reproducible paired benchmark for the best passing sampled
  transition-refinement point. Matched Oracle C anchors show errors below
  `1e-4` on both tested `DN_re` edges; a fixed-affinity 30-pair A/B shows a
  latency tradeoff between edges. The candidate remains validation-only; see
  the accuracy, Oracle C, and benchmark records.

- Updated the English and Chinese root READMEs to use the corrected
  `tol=1e-6` performance profile and to distinguish the 2026-09-12 accuracy
  baseline from later sampled `DN_re` edge checks. The documentation index now
  points readers to the matching protocol and evidence scope.
- Reduced repeated work in the Phase A benchmark by reusing measured candidate
  spectra and cached dex differences; simplified duplicate exception types in
  the packaged Fermi-Dirac cache loader. The Phase A report now uses the shared
  linear-spectrum relative-error convention for signal, full-grid, and
  transition metrics; the Cobaya adapter also caches parsed `eval_freqs` across
  samples to avoid repeated file reads. Solver outputs are unchanged.
- Aligned the signal masks in the fixed, randomized, and resolution-convergence
  validation scripts with the shared metric contract: the reference spectrum
  selects signal channels, while the candidate spectrum remains the relative-
  error denominator. This can change which near-threshold channels enter the
  reported signal metric; numerical solver outputs are unchanged.
- Reworked the shared linear-spectrum relative-error helper to evaluate the
  equivalent log-difference form with `expm1`. This avoids cancellation for
  small spectral differences and false zero errors when both linear spectra
  underflow, while preserving the candidate-denominator definition. The random
  frequency diagnostic now reuses the shared relative and dex helpers and only
  materializes linear Omega values for its ten printed sample channels. Its
  reporting threshold is now clearly distinguished from solver convergence,
  and its default plot location and near-zero metric scope are documented.

- Corrected the fixed-resource profile harness to use the formal `fast` preset
  tolerance by default and record the full solver configuration, Numba layer,
  and resource telemetry. The refreshed `tol=1e-6` profile reports warm
  medians of `5.987/10.049/9.179/9.458 ms` for
  default/high-T/stiff/high-kappa; the older
  `tol=1e-4` measurements are now labeled historical. See
  [benchmarks](docs/benchmarks.md) and [accuracy](docs/accuracy.md).
- Routed Cobaya's reference engine through the model dispatcher so the reported
  engine and evaluation status match the calculation; reference runs are now
  counted once. The theory docstring also lists the absolute
  `Delta_Neff_GW_error` estimate and its `NaN` unavailable value.
- Exposed unsupported fast error-budget configurations and made
  `auto_escalate` move to the validation/reference path when calibration is
  missing. Failed retries clear stale error values and report a distinct
  failure status. See [accuracy](docs/accuracy.md) and [Cobaya](docs/cobaya.md).
- Cached fixed PTA likelihood preparation and replaced per-step tensor-output
  modulo checks with a next-index counter. The documented local A/B records
  show identical outputs; the 25-call direct likelihood medians changed from
  34.033 to 27.568 ms (IPTA) and 4.574 to 0.203 ms (NANOGrav), while the
  six-regime tensor sample-index full-path median improved 5.44% over 80 paired
  rounds. See [benchmarks](docs/benchmarks.md); these results have the scope
  limits stated there.
- Added backend-aware solver call synchronization. Workqueue calls remain
  serialized, while explicit-config TBB calls improved two-call batch median
  time 14.93% in the documented local comparison; this is not a single-call
  speedup claim.
- Extended the Oracle C full-grid audit to two `DN_re` edges. The formal fast
  profile misses the scoped `2e-4` accuracy gate there; a matched fine-step
  candidate passes but is 11.0–15.8x slower. A matched transition-refined
  candidate also passes at 8.2–16.3x latency. Both remain validation candidates;
  the production preset is unchanged. See [accuracy](docs/accuracy.md) and the
  [Oracle C report](docs/oracle_c_wkb_assessment.md).
- Made transition-refinement radius and factor explicit in immutable solver
  configuration and both Oracle reproduction tools. A matched radius-0.25
  candidate passes the two `DN_re` edge gates and is 1.29x/1.42x faster than
  radius 2.0, while still 6.3–11.4x slower than standard fast. A 500,000-node
  estimate cap bounds the configured refinement grid; the maximum local
  refinement factor is now 128 for safely coarser outer grids. The candidate
  remains validation-only; Oracle output directories are created automatically.
- Derived config override keys from the dataclass and reused `replace()` for
  immutable snapshots, removing repeated field-copy lists. Moved the
  transition-refinement explanation into the public solver API docstring and
  aligned its English and Chinese descriptions.
- Replaced workstation-specific attachment paths in two historical plans with
  portable provenance notes; their original plan content and outcome records
  remain in the repository.
- Added a matched Oracle sweep that keeps the local transition step fixed while
  coarsening the outer grid. Narrowing the passing `h=0.02/factor=64` window
  to radius `0.0625` reaches `1.14e-4/8.30e-5` error at `6.49/8.74 ms`;
  radii `0.03125`, `0.015625`, and `0.0078125` also pass. A fixed-affinity
  paired timing found no stable speed winner between the two narrowest radii.
  `h=0.04/factor=64` is faster but fails both gates. A wider
  `h=0.04/radius=0.25/factor=128` run
  passes only the high edge and has asymmetric repeated medians of
  `47.15/7.75 ms`; it is also rejected. No preset change.
- Batched `H2_vec` density-spline evaluations on the transition-refined path.
  The full-path A/B is bitwise identical and 3.86–4.09x faster on the two
  `DN_re` edge cases. The implementation retains scalar `math.exp` to avoid
  platform-level ULP differences; it does not change the formal exact-kink
  path. See [benchmarks](docs/benchmarks.md).
- Added Chinese explanations to standalone comment sections across the solver,
  reference, grid, background, and Cobaya modules. The grouped-comment scan
  found no remaining English-only blocks in the package, and the code ASTs are
  unchanged. Clarified that the archived `baseline_54d65e3.md` links are for
  historical discovery, and documented the rejected 10-point `paramsweep_ref`
  summary as scoped evidence.
- Removed an unused private phase-error coefficient and its stale reference to
  a nonexistent calibration function. Clarified the remaining changed inline
  comments in Chinese; numerical statements are unchanged.
- Clarified documentation navigation so the current accuracy summary is read
  separately from the dated 2026-09-12 baseline and older validation manifest;
  the 2026-10-22 two-edge probe is now bound to its actual source commits.
- Added Chinese explanations to active validation and plotting scripts for
  parameter branches, oracle node matching, worker reuse, and failure classes.
  Runtime statements are unchanged.
- A repository-wide grouped scan found no English-only substantive standalone
  comment blocks across 231 tracked Python files. Comment-only edits preserve
  runtime ASTs.
- Reviewed the remaining inline comments across solver, likelihood, benchmark,
  and test code. Explanatory comments now have Chinese context; conventional
  unit abbreviations and serialized field names remain unchanged.
- Expanded the bilingual experiment catalogue to cover additional dated
  oracle/tail, accuracy, profile, grid, kernel, and outer-solve artifact
  families, with explicit retention and interpretation rules.
- Linked the 2026-09-03 Numba parallel-kernel A/B and production output records
  to their report, and documented that the separate thread-scaling JSONs were
  removed by an earlier cleanup commit while their summary remains.

## Unreleased — 2026-09-25

- Clarified the current `fast` profile and aligned packaging instructions with
  the supported installation, verification, and release workflow. Expanded
  the Chinese README to cover the same core usage and evidence limits as the
  English guide. No solver behavior or numerical results changed.
- Added paired developer README guides for the package, Cobaya adapter, scripts,
  and tests. Updated core module docstrings/comments to describe current engine
  roles, legacy-setting boundaries, and calibrated-error scope in both English
  and Chinese. No runtime behavior or scientific data changed.
- Added paired guides for the LVK/PTA likelihood packages and their nested data
  resources, and linked the likelihood index from the Cobaya package guide.
  Likelihood code/YAML edits are documentation-only; resource files are unchanged.
- Removed a verified empty, unreferenced local MCMC smoke directory and recorded
  the cleanup in both experiment catalogues. Recent local chain data with
  incomplete provenance remains excluded from formal comparison claims.

## Unreleased — 2026-09-24

- Refreshed the bilingual project and documentation indexes to describe the
  single user-facing `fast` profile, its latest fixed-resource runtime
  profile, current MCMC comparison, and the commit-bound limits of older
  validation artifacts. No solver behavior or numerical evidence was changed.

## Unreleased — 2026-09-12

- Removed redundant work in the fast Python preparation layer with no change to
  any numerical output: the goal-grid builder no longer sorts and
  de-duplicates the native frequency grid twice, the `N`-independent
  `H2_vec`/`raw_last`/`raw_re`/`Delta_f` terms in `grid_independent_freqs` are
  hoisted out of the per-`N` closure, freshly allocated `Ogw`/`Oj`/`Opgw`
  buffers skip an immediate `fill(0.0)`, and repeated `derived_param` property
  evaluations are bound once.  All six probe regimes keep bit-identical
  `f`/`log10OmegaGW`/`DN_gw`/`g2`/`w2` SHA256 digests (648 checked comparisons)
  and unchanged `converged`/`n_freq`/`fast_failure_reason`.
- Paired A,B,B,A warm-runtime A/B at a fixed affinity
  (`docs/fast_pyoverhead_ab_symmetric.json`): the default point improves `5.4%`
  at 20 threads and `8.2%` at 2 threads, and no probe regime regresses by more
  than the `2%` noise margin in the symmetric-order median.

## Unreleased — 2026-09-11

- Switched the single formal fast mode's default `frequency_quadrature` to
  shape-preserving PCHIP.  Against the Oracle C WKB anchor the four named
  points move from `4.31e-4 / 1.24e-2 / 7.50e-4 / 1.29e-3` (Simpson) to
  `1.07e-5 / 1.66e-4 / 9.88e-6 / 1.15e-5` (PCHIP); `simpson` stays selectable
  and remains the embedded reference estimator.  The PCHIP hot path degrades a
  non-finite integrand to the established `nonfinite` abort/restore guard.
- Regenerated the HEAD validation evidence (runtime matrix, stability screen,
  same-grid reference, Oracle C anchors, error-budget probes, estimator
  coverage) for the new default and refreshed the validation manifest.
- Re-calibrated the `fast`/`ultra-fast` frequency-quadrature error budget from
  `1.0e-3` to `2.0e-4` (the measured PCHIP residual against the Oracle C WKB
  anchor is `1.07e-5`–`1.66e-4`), and refreshed `README.md`, `README_zh.md` and
  `docs/fast_v02_audit_report.md`.  The same-grid reference residual is
  `2.93e-4`; it is reported as a combined tail/transfer and quadrature term, so
  the branch stays `PARTIALLY VERIFIED`.

## Unreleased — 2026-09-10

- Separated native `eval_freqs` output nodes from bolometric integration support
  nodes, preserving the self-consistent `DN_gw` invariant across six probe
  regimes.
- Added Simpson/PCHIP embedded quadrature telemetry and current-HEAD evidence
  for frequency quadrature, node-count convergence, runtime, and fast-only
  stability. PCHIP and DN-driven refinement remain opt-in/rejected where the
  oracle evidence is not sufficient.
- Refreshed the release evidence for outer-reuse safety and recorded the
  rejected 89-node/PCHIP candidate; the manifest remains explicitly
  `PARTIALLY VERIFIED` because the DN and stable-runtime gates are not met.

## 0.2.1 — 2026-09-04

- Renamed the distribution metadata and Python import package to
  `stiffgwpy_fast`; PyPI normalizes the distribution display name to
  `stiffgwpy-fast`.
- Updated all Cobaya qualified names, resource paths, scripts, tests, CI
  commands, and bilingual installation instructions to use the new name.

## 0.2.0 — 2026-09-04

- Added immutable per-call fast-solver configuration with `fast` plain-grid as
  the default high-level profile; `production` remains an explicit precision mode.
- Added package-resource loading and an installed-wheel smoke test covering core,
  Cobaya, LIGO, and PTA data.
- Added Python-version CI coverage, a scheduled/on-demand slow numerical job,
  coverage reporting, and validation-manifest schema checks.
- Fixed optional-Cobaya test collection and capped adapter thread requests to
  the runner's available Numba budget.
- Fixed the installed-wheel smoke probe so optional Cobaya resources are
  checked only when Cobaya is installed.
- Kept the independent continuous-sigma reference as the precision anchor;
  existing `NOT VERIFIED` and `FAIL` validation results remain unchanged.

Verification for this release line is documented in
`docs/engineering_audit.md` and `docs/reproducibility.md`.
