# Experiment catalogue

[中文目录](experiment_catalog_zh.md) | English

Status: current map of research notes and run artifacts

Date: 2026-09-25

Code version: each artifact is commit-bound. The validation manifest is dated `2026-09-17`; newer profiles and reports are linked by their own run commit.

This catalogue groups the experiment notes in `docs/` by scientific question.
The Markdown reports explain protocols and decisions; JSON/JSONL files retain
the measurements. A candidate result is not a production result merely
because it is faster or numerically close on one case.

## How to read the evidence

1. Use [`benchmarks.md`](benchmarks.md) for the latest fixed-resource runtime
   profile, then read [`fast_v02_audit_report.md`](fast_v02_audit_report.md)
   and [`validation/validation_manifest.json`](validation/validation_manifest.json)
   as date-bound scoped accuracy evidence and explicit gaps.
2. Use the topic reports below for hypotheses, controls, negative findings,
   and the reason a method was accepted or rejected.
3. Use the linked JSON/JSONL artifacts for pointwise values, environment,
   commit provenance, repeat counts, and replay information.
4. Treat `SMOKE`, `SPIKE`, and `NOT VERIFIED` as exploratory or incomplete.
   A failure, explicit physical guard, or rejected candidate is retained as
   evidence; it is not a solver failure unless the report classifies it so.

## Research themes

| Theme | Main question and synthesis | Reports and data |
|---|---|---|
| Physics model and solver contract | Defines the cosmological model, output quantities, solver modes, accuracy limits, and Cobaya adapter. These documents describe the contract; they do not certify the full parameter space. | [`physics.md`](physics.md), [`numerical_method.md`](numerical_method.md), [`accuracy.md`](accuracy.md), [`cobaya.md`](cobaya.md), [`reproducibility.md`](reproducibility.md) |
| Parameter-space validation and oracle agreement | Separates successful comparisons, physical guards, and numerical failures. The artifact-backed 400-point LHS screen records 254 successful points, 146 explicit guards, and no numerical failures; the 240-point Sobol run records 212 successes and 28 guards. These are dated screens, not full certification of the formal fast profile. | [`parameter_validation.md`](parameter_validation.md), [`parameter_validation/parameter_validation_report.md`](parameter_validation/parameter_validation_report.md), [`paramsweep_plain/validation_summary.md`](paramsweep_plain/validation_summary.md), [`paramsweep_z8/validation_summary.md`](paramsweep_z8/validation_summary.md), [`paramsweep_z8b/validation_summary.md`](paramsweep_z8b/validation_summary.md), [`parameter_validation/`](parameter_validation/), [`paramsweep_plain/`](paramsweep_plain/), [`paramsweep_z8/`](paramsweep_z8/), [`paramsweep_z8b/`](paramsweep_z8b/), [`paramsweep_ref/`](paramsweep_ref/), [`validation/`](validation/) |
| DN error budget and frequency quadrature | Early A/B work found the SciPy PCHIP path more accurate but initially too slow for the pre-registered runtime budget. Sharing/vectorization reduced its overhead; later WKB-anchored error decomposition attributed the larger residual to frequency quadrature, and the default-switch assessment records the subsequent PCHIP decision. The latest detailed accuracy audit is dated 2026-09-12; it reports scoped measurements and remains a historical evidence snapshot. | [`fast_true_error_assessment.md`](fast_true_error_assessment.md), [`fast_residual_decomposition_assessment.md`](fast_residual_decomposition_assessment.md), [`fast_quadrature_ab_assessment.md`](fast_quadrature_ab_assessment.md), [`fast_quadrature_reuse_assessment.md`](fast_quadrature_reuse_assessment.md), [`fast_quadrature_default_switch_assessment.md`](fast_quadrature_default_switch_assessment.md), [`fast_v02_audit_report.md`](fast_v02_audit_report.md), `quadrature_*`, `fast_quadrature_*`, and `dn_*` JSON artifacts |
| Independent oracles and tail handoff | Oracle C's higher-order WKB correction is verified for the documented comparison and explains the frozen-tail defect. Oracle B's phase-averaged and finite-window prototypes were not promoted because they do not yet provide an independent truth anchor. Prüfer is accepted as a standalone reference prototype but remains reference-only pending broader certification. | [`oracle_c_wkb_assessment.md`](oracle_c_wkb_assessment.md), [`oracle_b_phase_averaged_assessment.md`](oracle_b_phase_averaged_assessment.md), [`oracle_b_phase_window_assessment.md`](oracle_b_phase_window_assessment.md), [`oracle_prufer_assessment.md`](oracle_prufer_assessment.md), `oracle_*`, `tail_*`, and `wkb_*` JSON artifacts |
| Numerical-method candidates | The linear-z Magnus transfer-map candidate was rejected for production: its measured propagation time exceeded the midpoint baseline and did not meet the experiment gate. Other transfer, phase, Riccati, Bessel, Airy, sparse-frequency, and small-angle files are candidate-specific evidence and must be read with their recorded status and commit. | [`transfer_map_phase_assessment.md`](transfer_map_phase_assessment.md), `*transfer*`, `*phase*`, `*riccati*`, `*bessel*`, `*airy*`, `*frequency*`, and `*small_angle*` JSON artifacts |
| Performance and hot-path experiments | Benchmark documents distinguish warm runtime, cold/JIT cost, thread scaling, and accuracy gates. `profile_fast_breakdown_*.json` is the largest raw family (220 records in this 2026-09-24 inventory); it is a set of per-round/per-regime profiles, not 220 independent production improvements. The `derived_param_*`, `phi_s2_*`, `fd_lookup_*`, `grouped_soa_*`, and `outer_snapshot_*` families are repeated cache, workspace, and kernel candidates; compare digest/status fields and all regimes before accepting a speed claim. | [`benchmarks.md`](benchmarks.md), [`performance_comparison_20260903.md`](performance_comparison_20260903.md), [`baseline_54d65e3.md`](archive/baselines/baseline_54d65e3.md), `benchmark_*`, `profile_*`, `derived_param_*`, `phi_s2_*`, `fast_phi_*`, `fd_lookup_*`, `grouped_soa_*`, `outer_snapshot_*`, and `count_assembly_*` artifacts |
| Posterior and external-model comparison | The 2026-09-24 report compares three MCMC methods on shared LVK data. Current-fast sampling steps are about 2.9–4.9 times faster than SageNet+ across the three reported scenarios. The low-reheating case covers 0% of the LVK band for every method and cannot support a data-fit comparison. This is a sampling-step benchmark, not isolated solver runtime. The report's canonical chain archive contains four chains × 10,000 retained samples; a separate ignored local file under `mcmc/chains/` contains three chains × 2,000 samples and is an earlier, unreferenced run retained pending provenance review. | [`mcmc_sagenet_compare/report.md`](mcmc_sagenet_compare/report.md), [`mcmc_posterior/posterior_validation.md`](mcmc_posterior/posterior_validation.md), [`mcmc_sagenet_compare/`](mcmc_sagenet_compare/), [`mcmc_posterior/`](mcmc_posterior/) |
| Engineering and test coverage | These are maintenance audits, not physics-result reports. The engineering audit remains cross-linked from the root README, changelog, and reproducibility guide, so it is retained in the active documentation tree. | [`engineering_audit.md`](engineering_audit.md), [`test_coverage_matrix.md`](test_coverage_matrix.md), [`test_duplication_audit.md`](test_duplication_audit.md) |

## Raw run families by location

| Location or filename family | What it contains | Handling |
|---|---|---|
| `paramsweep_*`, `parameter_validation/`, `validation/` | Matched-grid, parameter-space, convergence, and oracle-independence records | Retain: referenced by reports, manifests, or reproducibility instructions |
| `oracle_prufer_fullgrid*.json`, `oracle_prufer_*_sobol*.json` | Full-native-grid Prüfer comparisons and edge/Sobol samples | Retain: supports the reference-only certification boundary |
| `profile_fast_breakdown_*.json`, `benchmark_*.json` | Round-specific performance, thread, stability, and stage timing runs | Retain as dated measurements; use the matching report/commit before comparison |
| `derived_param_*.json`, `fd_lookup_*.json`, `grouped_soa_*.json`, `outer_snapshot_*.json` | Cache/workspace/lookup and outer-solve optimization variants | Exploratory family; compare status/digest equality and regime coverage |
| `phi_s2_*`, `fast_phi_*`, `phase_*` | Phase and `phi_s2` implementation candidates, including smoke runs | Exploratory; do not infer production promotion from a smoke result |
| `*transfer*`, `*wkb*`, `*riccati*`, `*bessel*`, `*airy*`, `*small_angle*`, `*frequency*` | Numerical-method prototypes, oracle studies, and follow-up rounds | Keep each recorded status; rejected and reference-only results are part of the research trail |
| `mcmc_posterior/`, `mcmc_sagenet_compare/`, `mcmc/` | Posterior validation, comparison report, figures, and local chain output | Retain current reports and their inputs; large/local chain files remain governed by `.gitignore`. The canonical comparison archive is `mcmc_sagenet_compare/posterior_chains_20260924.npz` (4 × 10,000 per method/context); `mcmc/chains/sagenet_compare_20260924.npz` is a distinct 3 × 2,000 run, not a duplicate. Keep it local and out of formal claims until its run settings/provenance are recorded. |

## Archived material and cleanup decisions

- [`archive/baselines/baseline_54d65e3.md`](archive/baselines/baseline_54d65e3.md)
  is an earlier 2026-09-11 two-thread snapshot with no active in-repository
  path references. It remains available as historical evidence and is no
  longer presented as a current baseline.
- Three Prüfer raw `.log` captures dated 2026-09-11 are kept locally under
  `archive/oracle_prufer/raw_logs/`. They are ignored by Git's repository-wide
  `*.log` rule and have no active path references; the later JSON artifacts
  and assessment are the shared, navigable records. The logs are retained
  locally because they were generated from an earlier commit and are not
  byte-for-byte substitutes for the later JSON. A fresh clone will not
  contain these local-only logs; see [`archive/README.md`](archive/README.md).
- No dated parameter-validation or posterior data were deleted: the oldest
  2026-09-02/03 records are still cited by the current validation/report
  chain. No current-round experiment output was discarded by age alone.
- A 2026-09-25 review found the ignored local file
  `mcmc/chains/sagenet_compare_20260924.npz` (3 chains × 2,000 samples) differs
  from the report archive (4 × 10,000); neither the report nor benchmark script
  references it. Keep this recent independent run locally, but exclude it from
  formal comparison claims until its settings and provenance are recorded. The
  old empty `mcmc_smoke2/chains/` directory had no files, history, or references
  and was removed during the final cleanup on 2026-09-25.

## Tracked Markdown inventory

This list covers the experiment notes and supporting documents currently
under `docs/`; use the thematic table above to follow each thread.

- Core: `physics.md`, `numerical_method.md`, `accuracy.md`, `parameter_validation.md`, `cobaya.md`, `benchmarks.md`, `reproducibility.md`.
- Accuracy and quadrature: `fast_v02_audit_report.md`, `fast_true_error_assessment.md`, `fast_residual_decomposition_assessment.md`, `fast_quadrature_ab_assessment.md`, `fast_quadrature_reuse_assessment.md`, `fast_quadrature_default_switch_assessment.md`, `performance_comparison_20260903.md`.
- Oracles and methods: `oracle_b_phase_averaged_assessment.md`, `oracle_b_phase_window_assessment.md`, `oracle_c_wkb_assessment.md`, `oracle_prufer_assessment.md`, `transfer_map_phase_assessment.md`.
- Parameter/posterior: `parameter_validation/parameter_validation_report.md`, `paramsweep_plain/validation_summary.md`, `paramsweep_z8/validation_summary.md`, `paramsweep_z8b/validation_summary.md`, `mcmc_posterior/posterior_validation.md`, `mcmc_sagenet_compare/report.md`.
- Engineering and plans: `engineering_audit.md`, `test_coverage_matrix.md`, `test_duplication_audit.md`, `superpowers/plans/2026-09-10-fast-v02-dn-gw.md`, `superpowers/plans/2026-09-23-fast-v02-tensor-propagation.md`.
- Navigation and history: `README.md`, `README_zh.md`, `experiment_catalog.md`, `experiment_catalog_zh.md`, `archive/README.md`, `archive/README_zh.md`, `archive/baselines/baseline_54d65e3.md`.
