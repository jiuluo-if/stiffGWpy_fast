# Experiment catalogue

[中文目录](experiment_catalog_zh.md) | English

Status: current map of research notes and run artifacts

Date: 2026-10-22

Code version: each artifact is commit-bound. The validation manifest is dated `2026-09-17`; newer profiles and reports are linked by their own run commit.

This catalogue groups the experiment notes in `docs/` by scientific question.
The Markdown reports explain protocols and decisions; JSON/JSONL files retain
the measurements. A candidate result is not a production result merely
because it is faster or numerically close on one case.

## How to read the evidence

1. Use [`benchmarks.md`](benchmarks.md) for the latest fixed-resource runtime
   profile and [`accuracy.md`](accuracy.md) for current scoped accuracy
   evidence. Read [`fast_v02_audit_report.md`](fast_v02_audit_report.md) and
   [`validation/validation_manifest.json`](validation/validation_manifest.json)
   as date-bound baselines with explicit coverage gaps.
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
| DN error budget and frequency quadrature | Early A/B work found the SciPy PCHIP path more accurate but initially too slow for the pre-registered runtime budget. Sharing/vectorization reduced its overhead; later WKB-anchored error decomposition attributed the larger residual to frequency quadrature, and the default-switch assessment records the subsequent PCHIP decision. The four-regime audit dated 2026-09-12 is a historical baseline; the newer two-edge `DN_re` accuracy scan and transition-refinement results are in `accuracy.md` and `oracle_c_wkb_assessment.md`. | [`fast_true_error_assessment.md`](fast_true_error_assessment.md), [`fast_residual_decomposition_assessment.md`](fast_residual_decomposition_assessment.md), [`fast_quadrature_ab_assessment.md`](fast_quadrature_ab_assessment.md), [`fast_quadrature_reuse_assessment.md`](fast_quadrature_reuse_assessment.md), [`fast_quadrature_default_switch_assessment.md`](fast_quadrature_default_switch_assessment.md), [`fast_v02_audit_report.md`](fast_v02_audit_report.md), `quadrature_*`, `fast_quadrature_*`, and `dn_*` JSON artifacts |
| Independent oracles and tail handoff | Oracle C's higher-order WKB correction agrees with the deep reference on five sampled full grids. The formal fast solver passes the `2e-4` comparison on three points but fails at two `DN_re` edges (`2.5e-3..2.7e-3`). Matched `h=0.0003125` and transition-refine candidates pass those two points. A radius-`0.125` sweep with fixed `0.0003125` local spacing passes at `h=0.0025`, `0.005`, and `0.02`; the fastest passing point is `h=0.02/factor=64` with errors `1.16e-4`/`8.25e-5` at `7.63/9.49 ms`. It remains `1.34x/2.09x` slower than standard fast. Coarsening to `h=0.04/factor=64` fails both gates. None is promoted. This is sampled evidence, not full-space certification. Oracle B's phase-averaged and finite-window prototypes were not promoted because they do not yet provide an independent truth anchor. Prüfer remains a standalone reference prototype. | [`oracle_c_wkb_assessment.md`](oracle_c_wkb_assessment.md), [`accuracy.md`](accuracy.md), [`oracle_c_wkb_corner_20261022/`](oracle_c_wkb_corner_20261022/), [`oracle_b_phase_averaged_assessment.md`](oracle_b_phase_averaged_assessment.md), [`oracle_b_phase_window_assessment.md`](oracle_b_phase_window_assessment.md), [`oracle_prufer_assessment.md`](oracle_prufer_assessment.md), `oracle_*`, `tail_*`, and `wkb_*` JSON artifacts |
| Numerical-method candidates | The linear-z Magnus transfer-map candidate was rejected for production: its measured propagation time exceeded the midpoint baseline and did not meet the experiment gate. Other transfer, phase, Riccati, Bessel, Airy, sparse-frequency, and small-angle files are candidate-specific evidence and must be read with their recorded status and commit. | [`transfer_map_phase_assessment.md`](transfer_map_phase_assessment.md), `*transfer*`, `*phase*`, `*riccati*`, `*bessel*`, `*airy*`, `*frequency*`, and `*small_angle*` JSON artifacts |
| Performance and hot-path experiments | Benchmark documents distinguish warm runtime, cold/JIT cost, thread scaling, and accuracy gates. The corrected 2026-10-22 formal `fast` profile uses solver source at commit `6707572e73d576a0ccf294172579cedc26e5e394`, preset `tol=1e-6`, and 25 warm repeats for default, high-T, stiff, and high-kappa; the 2026-10-21 profile at `cdbdb387479e11cd87a38d76a0e58348ca13e70f` omitted the low-level `tol` argument and used `1e-4`, so it remains lower-tolerance historical evidence. Both profiles separate the first process call; see `benchmarks.md`. The PTA likelihood-cache experiment measures direct-call latency for the IPTA and NANOGrav adapters and reports exact output equality on fixed inputs plus 10 varied nuisance/theory cases per adapter; these are not end-to-end MCMC speedups. The 2026-10-21 tensor sample-index optimization was compared against commit `68b17ac` over six regimes and 80 paired rounds; its raw record and rerun script are linked in `benchmarks.md`. A separate two-call TBB throughput comparison is recorded there; it does not claim a single-evaluation speedup. The 2026-10-22 reference-spectrum quadrature experiment cut `integrate_spectrum` function calls by about 25% across synthetic and formal-fast-shaped input and fixed the coarse-grid endpoint on even-size grids; its timing is for the integration helper only, not full solver speed. The tracked `profile_fast_breakdown_*.json` family contains 228 files as of 2026-10-22. These are per-round/per-regime profiles, not 228 independent production improvements. The `derived_param_*`, `phi_s2_*`, `fd_lookup_*`, `grouped_soa_*`, and `outer_snapshot_*` families are repeated cache, workspace, and kernel candidates; compare digest/status fields and all regimes before accepting a speed claim. | [`benchmarks.md`](benchmarks.md), [`benchmark_pta_likelihood_cache_20261021.json`](benchmark_pta_likelihood_cache_20261021.json), [`benchmark_pta_likelihood_cache.py`](../scripts/benchmark_pta_likelihood_cache.py), [`benchmark_reference_quad_breakpoints_20261022.json`](benchmark_reference_quad_breakpoints_20261022.json), [`benchmark_reference_quad_breakpoints.py`](../scripts/benchmark_reference_quad_breakpoints.py), [`performance_comparison_20260903.md`](performance_comparison_20260903.md), [`baseline_54d65e3.md`](archive/baselines/baseline_54d65e3.md), `benchmark_*`, `profile_*`, `derived_param_*`, `phi_s2_*`, `fast_phi_*`, `fd_lookup_*`, `grouped_soa_*`, `outer_snapshot_*`, and `count_assembly_*` artifacts |
| Posterior and external-model comparison | The 2026-09-24 report compares three MCMC methods on shared LVK data. Current-fast sampling steps are about 2.9–4.9 times faster than SageNet+ across the three reported scenarios. The low-reheating case covers 0% of the LVK band for every method and cannot support a data-fit comparison. This is a sampling-step benchmark, not isolated solver runtime. The report's canonical chain archive contains four chains × 10,000 retained samples; a separate ignored local file under `mcmc/chains/` contains three chains × 2,000 samples and is an earlier, unreferenced run retained pending provenance review. | [`mcmc_sagenet_compare/report.md`](mcmc_sagenet_compare/report.md), [`mcmc_posterior/posterior_validation.md`](mcmc_posterior/posterior_validation.md), [`mcmc_sagenet_compare/`](mcmc_sagenet_compare/), [`mcmc_posterior/`](mcmc_posterior/) |
| Engineering, tests, and historical plans | These are maintenance records, not physics-result reports. The engineering audit remains cross-linked from the root README, changelog, and reproducibility guide. The 2026-09-10 DN/runtime plan is a historical plan whose same-grid accuracy and `<4 ms` runtime gates remain open in the current evidence. The 2026-09-23 tensor-propagation plan is also historical: its exact-exp candidate was not adopted, while the later `next_assembly_k` implementation and measured outcome are recorded in the 2026-10-21 benchmark. | [`engineering_audit.md`](engineering_audit.md), [`test_coverage_matrix.md`](test_coverage_matrix.md), [`test_duplication_audit.md`](test_duplication_audit.md), [`superpowers/plans/2026-09-10-fast-v02-dn-gw.md`](superpowers/plans/2026-09-10-fast-v02-dn-gw.md), [`superpowers/plans/2026-09-23-fast-v02-tensor-propagation.md`](superpowers/plans/2026-09-23-fast-v02-tensor-propagation.md), [`accuracy.md`](accuracy.md), [`benchmark_kernel_assembly_stride_20261021.json`](benchmark_kernel_assembly_stride_20261021.json) |

## 2026-10-22 transition-refined `H2_vec` optimization

The transition-refined path called the SciPy density spline once per grid
node. Batching those queries preserved all solver output arrays exactly and
reduced warm median runtime by `3.86x`–`4.09x` on the two `DN_re` edges.
This is an implementation-speed result for the validation path; it does not
change the formal exact-kink profile or add accuracy coverage. See the
[benchmark summary](benchmarks.md#transition-refined-h2_vec-batching-2026-10-22),
[raw timing record](benchmark_h2_vectorization_20261022.json), and
[reproduction script](../scripts/benchmark_h2_vectorization.py).

## 2026-10-22 matched local-refinement sweep

Three radius-`0.125` configurations preserve a `0.0003125` local step while
coarsening the outer background grid. At `h=0.02/factor=64`, both sampled
`DN_re` edges pass the `2e-4` gate with medians `7.63/9.49 ms`; the candidate
is still `1.34x/2.09x` slower than standard fast. Increasing to `h=0.04` cuts
the median by about a third but fails both accuracy gates. None is promoted. Protocol,
matched Oracle residuals, and raw records are in
[`oracle_c_wkb_assessment.md`](oracle_c_wkb_assessment.md) and the linked
`oracle_c_wkb_corner_20261022/transition_h*` directories.
The wider `radius=0.25` run at factor 64 fails both edges; increasing to factor
128 restores the same local spacing as the `h=0.02` pass but still fails at the
low edge (`4.51e-4`). Its repeated `47.15/7.75 ms` medians are asymmetric. Both
configurations are recorded as rejected candidates in the assessment report
and raw-data family.
Narrowing the passing `h=0.02/factor=64` window to radius `0.0625` improves the
errors to `1.14e-4/8.30e-5` and medians to `6.49/8.74 ms`; this remains slower
than standard fast and is not promoted.
Halving the radius again to `0.03125` passes both gates at `1.12e-4/8.29e-5`,
with `6.41/8.09 ms` medians. The unpinned separate-run timing difference is
small; no stable gain over radius `0.0625` is established.
The fixed-affinity paired comparison of radii `0.015625` and `0.0078125` also
found no stable speed winner; see the raw paired record linked from the Oracle
C assessment.

## Raw run families by location

| Location or filename family | What it contains | Handling |
|---|---|---|
| `paramsweep_*`, `parameter_validation/`, `validation/` | Matched-grid, parameter-space, convergence, and oracle-independence records | Retain: referenced by reports, manifests, or reproducibility instructions |
| `oracle_a_*`, `oracle_b_*`, `oracle_prufer_*`, `oracle_tail_*`, `tail_*`, `wkb_*` | Same-grid Oracle A/B/C, Prüfer reference comparisons, phase/tail sensitivity, WKB corrections, and diagnostic subsets | Retain as scoped oracle and error-attribution evidence; Prüfer stays reference-only and no artifact implies full-space certification |
| `phase_recurrence_oracle_20260917/` | Five Oracle C WKB records for default, high-kappa, high-T, low-T, and stiff points; each contains 12 selected frequencies and was generated from commit `58af509` | Retain as a sparse phase-recurrence sensitivity subset; it is not a native full-grid comparison or parameter-space certification |
| `validation_h005_round_20260916/`, `validation_h0075_round_20260916/`, `validation_h0075_round_20260916b/` | Three dated 10-point step-size scans from source commit `1bf840d`; `h=0.0075` first recorded 10 `fast_exception` rows, while the completed `h=0.005` and rerun `h=0.0075` summaries fail the signal, transition, and `DN_gw` acceptance gates | Retain as rejected historical evidence; completed summaries report maximum signal relative errors around `6.78e-3` and maximum `DN_gw` relative errors around `3.60e-3`, not current accuracy evidence |
| `profile_fast_breakdown_*.json`, `benchmark_*.json` | Round-specific performance, thread, stability, and stage timing runs | Retain as dated measurements; use the matching report/commit before comparison |
| Other `profile_*.json`, `*_round*.json`, `*_head*.json`, `*_candidate*.json`, `*_spike*.json` | Per-round profiles and standalone implementation experiments; many record rejected or reference-only candidates | Retain the raw record with its originating report; determine status from the report, not the filename or date |
| `ab_baseline_njit_*.json`, `ab_candidate_njit_*.json`, `production_ab_*.txt` | Raw baseline/candidate measurements and production-output fingerprints for the 2026-09-03 Numba parallel-kernel change; thread scaling is summarized in the report, but its `scaling_njit_t*.json` samples were removed by cleanup commit `2619065` | Retain the remaining raw records with [`performance_comparison_20260903.md`](performance_comparison_20260903.md); neither they nor the report summary are a current runtime baseline |
| `fast_true_error_*`, `fast_residual_*`, `error_budget_probe_*`, `fast_quadrature_*`, `quadrature_*` | Dated Oracle error attribution, estimator calibration, and quadrature comparisons | Retain as scoped accuracy evidence; use `accuracy.md` and linked assessments for current interpretation |
| `derived_param_*.json`, `fd_lookup_*.json`, `grouped_soa_*.json`, `outer_snapshot_*.json` | Cache/workspace/lookup and outer-solve optimization variants | Exploratory family; compare status/digest equality and regime coverage |
| `gen_fast_*`, `h2_*`, `outer_*`, `kernel_*`, `prep_*`, `probe_*`, `canonical_fast_*`, `fast_pyoverhead_*` | Grid, background, outer-solve, kernel, and Python-overhead candidates across dated rounds | Retain accepted and rejected comparisons; check each report and resource record before using a timing result |
| `phi_s2_*`, `fast_phi_*`, `phase_*` | Phase and `phi_s2` implementation candidates, including smoke runs | Exploratory; do not infer production promotion from a smoke result |
| `*transfer*`, `*wkb*`, `*riccati*`, `*bessel*`, `*airy*`, `*small_angle*`, `*frequency*` | Numerical-method prototypes, oracle studies, and follow-up rounds | Keep each recorded status; rejected and reference-only results are part of the research trail |
| `mcmc_posterior/`, `mcmc_sagenet_compare/`, `mcmc/` | Posterior validation, comparison report, figures, and local chain output | Retain current reports and their inputs; large/local chain files remain governed by `.gitignore`. The canonical comparison archive is `mcmc_sagenet_compare/posterior_chains_20260924.npz` (4 × 10,000 per method/context); `mcmc/chains/sagenet_compare_20260924.npz` is a distinct 3 × 2,000 run, not a duplicate. Keep it local and out of formal claims until its run settings/provenance are recorded. |

## Archived material and cleanup decisions

- [`paramsweep_ref/validation_summary.json`](paramsweep_ref/validation_summary.json) is a 10-point summary generated from commit `f87e969` with `z_tail=7`, `rtol=1e-9`, and `freq_res=2`. Its signal and integrated `DN_gw` acceptance gates are false. No code or benchmark script consumes it; sibling JSONL files retain the pointwise inputs. Keep it as a compact record of this rejected/partial sweep, not as current certification.
- [`archive/baselines/baseline_54d65e3.md`](archive/baselines/baseline_54d65e3.md)
  is an earlier 2026-09-11 two-thread snapshot. No solver or benchmark script
  consumes it; the README and catalog link to it only for historical discovery.
  It remains evidence of that dated run and is not a current baseline.
- [`archive/reference/head_20260910/benchmark_default.jsonl`](archive/reference/head_20260910/benchmark_default.jsonl)
  is a separate single-point fast/reference record from commit `767056d`; no
  code or benchmark script consumes it. Keep it as dated point-level evidence,
  not as a replacement for the later two-thread snapshot or a current baseline.
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
