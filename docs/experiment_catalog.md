# Experiment catalogue

[中文目录](experiment_catalog_zh.md) | English

Status: current map of research notes and run artifacts

Last full inventory: 2026-10-22

Focused refresh: 2026-11-04 (formal-fast profile, Prüfer provenance, transition-refinement frontier, LVK likelihood cache, and res4 root closure)

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
| Independent oracles and tail handoff | Oracle C's higher-order WKB correction agrees with the deep reference on five sampled full grids. The formal fast solver passes the `2e-4` comparison on three points but fails at two `DN_re` edges (`2.5e-3..2.7e-3`). Matched `h=0.0003125` and transition-refine candidates pass those two points. A radius-`0.125` sweep with fixed `0.0003125` local spacing passes at `h=0.0025`, `0.005`, and `0.02`; the fastest passing point is `h=0.02/factor=64` with errors `1.16e-4`/`8.25e-5` at `7.63/9.49 ms`. It remains `1.34x/2.09x` slower than standard fast. Coarsening to `h=0.04/factor=64` fails both gates. The 2026-10-23 frontier candidate passes both sampled accuracy gates but has mixed paired timing (candidate/formal-fast ratios `0.812` at the low edge and `1.177` at the high edge), so it remains validation-only. Oracle B's phase-averaged and finite-window prototypes were not promoted because they do not yet provide an independent truth anchor. Prüfer is now available as an opt-in Reference backend: on the exact res4 Reference root/full 260-mode grid it differs from Cartesian by `1.34e-8` in `DN_gw` and at most `4.28e-5 dex` in spectrum, with four-mode paired aggregate speedup `4.27x`. The backend remains opt-in and is not global-domain certified. The author ARDC full-grid attempt exceeded its 45-minute CPU budget without a result. See [`accuracy.md`](accuracy.md#2026-11-04-opt-in-reference-prufer-backend) and [`reference_optimization_20261104.md`](reference_optimization_20261104.md) for the method assumptions, raw data, and limits. | [`oracle_c_wkb_assessment.md`](oracle_c_wkb_assessment.md), [`accuracy.md`](accuracy.md), [`reference_optimization_20261104.md`](reference_optimization_20261104.md), [`oracle_c_wkb_corner_20261022/`](oracle_c_wkb_corner_20261022/), [`oracle_c_wkb_corner_20261023/`](oracle_c_wkb_corner_20261023/), [`benchmark_transition_frontier_ab_20261023.json`](benchmark_transition_frontier_ab_20261023.json), [`oracle_b_phase_averaged_assessment.md`](oracle_b_phase_averaged_assessment.md), [`oracle_b_phase_window_assessment.md`](oracle_b_phase_window_assessment.md), [`oracle_prufer_assessment.md`](oracle_prufer_assessment.md), [`reference_prufer_production_res4_refroot_20261104.json`](reference_prufer_production_res4_refroot_20261104.json), [`reference_prufer_production_ab_smoke_20261104.json`](reference_prufer_production_ab_smoke_20261104.json), [`reference_ardc_res4_timeout_20261104.json`](reference_ardc_res4_timeout_20261104.json), `oracle_*`, `tail_*`, and `wkb_*` JSON artifacts |
| Numerical-method candidates | The linear-z Magnus transfer-map candidate was rejected for production: its measured propagation time exceeded the midpoint baseline and did not meet the experiment gate. Other transfer, phase, Riccati, Bessel, Airy, sparse-frequency, and small-angle files are candidate-specific evidence and must be read with their recorded status and commit. | [`transfer_map_phase_assessment.md`](transfer_map_phase_assessment.md), `*transfer*`, `*phase*`, `*riccati*`, `*bessel*`, `*airy*`, `*frequency*`, and `*small_angle*` JSON artifacts |
| Performance and hot-path experiments | Benchmark documents distinguish warm runtime, cold/JIT cost, thread scaling, and accuracy gates. The corrected 2026-10-22 formal `fast` profile uses solver source at commit `6707572e73d576a0ccf294172579cedc26e5e394`, preset `tol=1e-6`, and 25 warm repeats for default, high-T, stiff, and high-kappa; the 2026-10-21 profile at `cdbdb387479e11cd87a38d76a0e58348ca13e70f` omitted the low-level `tol` argument and used `1e-4`, so it remains lower-tolerance historical evidence. Both profiles separate the first process call; see `benchmarks.md`. The PTA likelihood-cache experiment measures direct-call latency for the IPTA and NANOGrav adapters and reports exact output equality on fixed inputs plus 10 varied nuisance/theory cases per adapter; these are not end-to-end MCMC speedups. The 2026-10-23 LVK likelihood benchmark on the packaged 54,593-bin table shows 1.90x–9.75x lower per-call latency on five overlapping synthetic support cases; the sixth case is below LVK data support and measures only the unsupported-bin reduction path (`25.2x`). Maximum absolute log-likelihood change is `1.21e-10`; this is not a full MCMC speedup. See `benchmarks.md`. The 2026-10-21 tensor sample-index optimization was compared against commit `68b17ac` over six regimes and 80 paired rounds; its raw record and rerun script are linked in `benchmarks.md`. A separate two-call TBB throughput comparison is recorded there; it does not claim a single-evaluation speedup. The 2026-10-22 reference-spectrum quadrature experiment cut `integrate_spectrum` function calls by about 25% across synthetic and formal-fast-shaped input and fixed the coarse-grid endpoint on even-size grids; its timing is for the integration helper only, not full solver speed. The tracked `profile_fast_breakdown_*.json` family contains 228 files as of 2026-10-22. These are per-round/per-regime profiles, not 228 independent production improvements. The `derived_param_*`, `phi_s2_*`, `fd_lookup_*`, `grouped_soa_*`, and `outer_snapshot_*` families are repeated cache, workspace, and kernel candidates; compare digest/status fields and all regimes before accepting a speed claim. | [`benchmarks.md`](benchmarks.md), [`benchmark_pta_likelihood_cache_20261021.json`](benchmark_pta_likelihood_cache_20261021.json), [`benchmark_pta_likelihood_cache.py`](../scripts/benchmark_pta_likelihood_cache.py), [`benchmark_lvk_likelihood_cache_20261023.json`](benchmark_lvk_likelihood_cache_20261023.json), [`benchmark_lvk_likelihood_cache.py`](../scripts/benchmark_lvk_likelihood_cache.py), [`benchmark_reference_quad_breakpoints_20261022.json`](benchmark_reference_quad_breakpoints_20261022.json), [`benchmark_reference_quad_breakpoints.py`](../scripts/benchmark_reference_quad_breakpoints.py), [`performance_comparison_20260903.md`](performance_comparison_20260903.md), [`baseline_54d65e3.md`](archive/baselines/baseline_54d65e3.md), `benchmark_*`, `profile_*`, `derived_param_*`, `phi_s2_*`, `fast_phi_*`, `fd_lookup_*`, `grouped_soa_*`, `outer_snapshot_*`, and `count_assembly_*` artifacts |
| Posterior and external-model comparison | The 2026-09-24 report compares three MCMC methods on shared LVK data. Current-fast sampling steps are about 2.9–4.9 times faster than SageNet+ across the three reported scenarios. The low-reheating case covers 0% of the LVK band for every method and cannot support a data-fit comparison. This is a sampling-step benchmark, not isolated solver runtime. The report's canonical chain archive contains four chains × 10,000 retained samples; a separate ignored local file under `mcmc/chains/` contains three chains × 2,000 samples and is an earlier, unreferenced run retained pending provenance review. | [`mcmc_sagenet_compare/report.md`](mcmc_sagenet_compare/report.md), [`mcmc_posterior/posterior_validation.md`](mcmc_posterior/posterior_validation.md), [`mcmc_sagenet_compare/`](mcmc_sagenet_compare/), [`mcmc_posterior/`](mcmc_posterior/) |
| Engineering, tests, and historical plans | These are maintenance records, not physics-result reports. The engineering audit remains cross-linked from the root README, changelog, and reproducibility guide. The 2026-09-10 DN/runtime plan is a historical plan whose same-grid accuracy and `<4 ms` runtime gates remain open in the current evidence. The 2026-09-23 tensor-propagation plan is also historical: its exact-exp candidate was not adopted, while the later `next_assembly_k` implementation and measured outcome are recorded in the 2026-10-21 benchmark. | [`engineering_audit.md`](engineering_audit.md), [`test_coverage_matrix.md`](test_coverage_matrix.md), [`test_duplication_audit.md`](test_duplication_audit.md), [`superpowers/plans/2026-09-10-fast-v02-dn-gw.md`](superpowers/plans/2026-09-10-fast-v02-dn-gw.md), [`superpowers/plans/2026-09-23-fast-v02-tensor-propagation.md`](superpowers/plans/2026-09-23-fast-v02-tensor-propagation.md), [`accuracy.md`](accuracy.md), [`benchmark_kernel_assembly_stride_20261021.json`](benchmark_kernel_assembly_stride_20261021.json) |

## 2026-10-22 TBB chunk-size scheduling screen

The TBB screen covered five chunk sizes across four cases. Output digests match
the formal workqueue reference, but the fastest chunk differs by case and the
small timing differences are unpaired; no global scheduling change was
selected. The 20 preliminary JSON files were consolidated into one traceable
summary. See [`benchmarks.md`](benchmarks.md#tbb-dynamic-chunk-size-screen-2026-10-22)
and [`benchmark_tbb_chunksize_screen_20261022.json`](benchmark_tbb_chunksize_screen_20261022.json).

## 2026-10-23 base-HEAD formal-fast profile

A five-regime profile at base HEAD `ccea67b`, before the exact horizon-boundary
correction, reports warm medians of `5.984/6.106/9.346/9.970/9.529 ms` for
default, low-T, high-T, stiff, and high-kappa. Every case retained identical
output digests across 26 calls. The tensor propagation stage accounts for
about 49% of default runtime and 54%–58% in the three slow cases. This is a
dated base-HEAD profile, not an A/B speed claim or new accuracy evidence; the
full protocol and raw records are in the
[benchmark report](benchmarks.md#2026-10-23-base-head-formal-fast-profile).

## 2026-10-23 PTA theory-frequency support boundary

IPTA and NANOGrav avoid spline construction on an empty theory-knot set and
avoid extrapolating below the first supported knot. With configured packaged
frequencies and theory interval `[-7, -6]`, all 9 IPTA bins and 14 NANOGrav
bins fall below support and receive zero primordial power; the SMBHB component
remains active. Full-support and exact single-knot helper probes match the
prior formula exactly. See the
[benchmark note](benchmarks.md#pta-theory-frequency-support-boundary-2026-10-23)
and [source-bound probe record](pta_frequency_support_boundary_20261023.json).

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
C assessment and the
[paired timing data](benchmark_transition_radius_paired_20261022.json).

## 2026-10-23 transition-refinement frontier

The `h=0.0296875`, radius-`0.015625`, factor-`95` candidate passes both
sampled `DN_re` edges below `1e-4`; the tested `h=0.03` neighbor fails the high
edge. A 30-pair fixed-affinity comparison finds the candidate faster on the
low edge and slower on the high edge, so it remains validation-only. Matched
Oracle records, accuracy interpretation, and raw paired timings are in
[`oracle_c_wkb_assessment.md`](oracle_c_wkb_assessment.md),
[`accuracy.md`](accuracy.md), and
[`benchmark_transition_frontier_ab_20261023.json`](benchmark_transition_frontier_ab_20261023.json);
the replay driver is [`benchmark_transition_frontier_ab.py`](../scripts/benchmark_transition_frontier_ab.py).

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
  `archive/oracle_prufer/raw_logs/`. Their leading JSON records identify source
  commit `91598e7`. They are ignored by Git's repository-wide `*.log` rule and
  have no active code or benchmark-script consumers; the later JSON artifacts
  and assessment are the shared, navigable records. The logs are retained
  locally because they are from an earlier commit and are not byte-for-byte
  substitutes for the later JSON. A fresh clone will not contain these
  local-only logs; see [`archive/README.md`](archive/README.md).
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
- Fast domain pressure screen: [`fast_domain_pressure_20261030.md`](fast_domain_pressure_20261030.md) records the historical 514-point screen and current-head rechecks. The 2026-11-04 402-point boundary run reports 278 success statuses but only two analytic exact-zero accuracy certifications; the other 276 finite outputs remain uncertified. One `A_s=100` point was explicitly stopped after 240.8 CPU seconds, while a prior slow coupled point now returns finite output with an uncertified local error estimate. The 128-point Sobol/LHS pilots each find four finite spectra, all uncertified. Raw records: [`fast_domain_pressure_boundary_20261104.jsonl`](fast_domain_pressure_boundary_20261104.jsonl), [`fast_domain_pressure_sobol128_20261104.jsonl`](fast_domain_pressure_sobol128_20261104.jsonl), and [`fast_domain_pressure_latin128_20261104.jsonl`](fast_domain_pressure_latin128_20261104.jsonl).
- Navigation and history: `README.md`, `README_zh.md`, `experiment_catalog.md`, `experiment_catalog_zh.md`, `archive/README.md`, `archive/README_zh.md`, `archive/baselines/baseline_54d65e3.md`.

## 2026-11-04 res4 root closure and aligned-grid scalar step

The old pre-terminal-step artifact with 53 modes above `1e-3 dex` does not
describe current HEAD. A fresh current-source decomposition finds 18 default
uniform-path modes above `5e-4 dex` and none above `1e-3`; the h=.02
transition-refined candidate is below `5e-4 dex` on all 260 native frequencies.
Its Reference fixed-point root was independently closed at relative map
residual `-7.92e-10`; Fast and Reference roots differ by `5.77e-5` relative.
An independent continuous-coefficient DOP853 replay on the 18 remaining
baseline edge modes agrees with refined Fast to max `2.56e-6 dex` and with
tail9 Reference to max `2.47e-5 dex`.

A fresh fixed-background 260-mode factorial run at the candidate DN finds all
three tested propagation preparations (uniform fast primitive, uniform
`sigma_exact`, and transition-refined) below `4.20e-4 dex` maximum against
tail9. The current full-root default discrepancy is therefore dominated by the
outer DN-map/root difference; the old 53-mode tensor-propagation interpretation
is not current evidence. A fresh 15-pair 8-thread benchmark measured candidate
P50/P95 at `62.0/68.7 ms` versus `67.7/70.7 ms` default. The three-variant
batch rate is `13.82 cases/s`; scalar-step vs array-step spectra differ by at
most `2.85e-14 dex`.

The second hard point, retained random sample11, also passes on the current
h=.02 path: its exact-input tail10 Reference DN-map residual is `-9.48e-6`
relative, and same-DN spectrum maximum error is `9.149e-5 dex` on 94 modes.
Raw Fast/Reference arrays and a rerunnable driver are linked in
[`accuracy.md`](accuracy.md).

For mode-level attribution, the repository now includes all 260 current-root
errors, all 260 same-background primitive/path comparisons, the historical 53
and current 18 continuous-coefficient DOP853 records, and the matching tail9
Reference array. Their provenance summary is
[`res4_mode_decomposition_20261104.json`](res4_mode_decomposition_20261104.json);
the detailed data files are linked from [`accuracy.md`](accuracy.md).

A targeted `freq_res=4/8/16` check at fixed `h=.005,z_tail=8` leaves the same
18 modes above `5e-4 dex`; increasing frequency density is not the main source
of their error. The raw comparison is
[`res4_default_frequency_density_20261104.json`](res4_default_frequency_density_20261104.json).

Commit `6f5f472` makes an omitted `sigma_exact` select the continuous
background primitive only on the opt-in transition-refined path; an explicit
false value and the regular Fast path retain the node primitive. The updated
res4 candidate has zero modes above `5e-4 dex` and a maximum error of
`3.456e-4 dex`; the sample11 exact-input Reference map residual is `-6.21e-6`
with maximum same-input spectrum error `7.057e-5 dex`. The paired runs remain
below 100 ms P50. Raw comparisons and reproduction commands are in
[`accuracy.md`](accuracy.md#2026-11-04-continuous-background-default-for-transition-refinement).
The exact-input res4 tail10 Reference map now also directly closes at the
continuous candidate root: map residual `+3.466e-5`, same-input spectrum max
`3.458e-4 dex` on 260 modes.

The current `A_s=1e20` trace now stops at its first bad `S2=exp(Psi)` value
before tensor propagation, reporting the primitive and index; no extended
finite spectrum is claimed. A 501-pair six-regime A/B confirms the diagnostic
guard preserves output exactly and keeps the worst P50/P95 ratio below `1.025`.
See the staged failure details and raw timing samples in
[`fast_domain_pressure_20261030.md`](fast_domain_pressure_20261030.md) and
[`background_primitive_guard_benchmark_20261104.json`](background_primitive_guard_benchmark_20261104.json).

Commit `671aa11` adds a guarded scalar terminal-step path for transition grids
whose interior intervals are uniform within floating-point roundoff. The
previous exact interval-array path is retained for real refinements and tiny
physical anchors. At the target point, it preserves the array-path root to
`1.67e-15` and spectra to `6.04e-14 dex`; the candidate meets the 100 ms
8-thread target. The focused code tests are `61 passed, 1 deselected`; GitHub
Actions run `41759168617` passed all nine jobs. Full point-level results and
resource measurements are in [`accuracy.md`](accuracy.md) and
[`res4_scalar_step_profile_20261104.json`](res4_scalar_step_profile_20261104.json).
The Reference spectrum used for reproduction is
[`res4_reference_secant_root_tail10_20261104.npz`](res4_reference_secant_root_tail10_20261104.npz),
[`res4_reference_secant_root_tail10_20261104.json`](res4_reference_secant_root_tail10_20261104.json),
and the benchmark drivers are
[`benchmark_res4_scalar_step.py`](../scripts/benchmark_res4_scalar_step.py)
and [`benchmark_res4_reference_root.py`](../scripts/benchmark_res4_reference_root.py).
These are sampled hard-point results, not global domain certification.

## 2026-11-04 Reference DOP853 dense-output removal

The Cartesian Reference now uses event states from `y_events` and disables
global DOP853 dense output. Two fixed res4 hard modes matched bitwise and ran
`1.275x/1.292x` faster with reduced sampled peak RSS. Raw paired samples and
reproduction instructions are in
[`reference_dense_output_benchmark_20261104.json`](reference_dense_output_benchmark_20261104.json),
[`reference_optimization_20261104.md`](reference_optimization_20261104.md),
and [`benchmark_reference_dense_output.py`](../scripts/benchmark_reference_dense_output.py).
After changing the Cartesian baseline, the Prüfer A/B was repeated on the same
four modes: aggregate speedup `3.349x`, per-mode P50/P95 `3.085x/3.990x`, and
maximum spectrum error `1.401e-5 dex`; see
[`reference_prufer_production_ab_postdense_20261104.json`](reference_prufer_production_ab_postdense_20261104.json).

The Reference background now evaluates its unchanged SciPy cubic coefficients
through a Numba Horner helper for expensive modes and retains a no-JIT scalar
path for cheap cold single-mode requests. Same-input res4 A/B runs were
`1.61x–1.76x` faster; four archived spectrum modes matched within
`1.93e-8 dex`. See
[`reference_spline_eval_benchmark_20261104.json`](reference_spline_eval_benchmark_20261104.json)
and [`benchmark_reference_spline_eval.py`](../scripts/benchmark_reference_spline_eval.py).
Cold low-frequency single-mode calls do not compile; batched high-frequency
spectra warm the kernel before worker dispatch. Four-worker outputs matched
serial bitwise on four selected modes, but the small parallel run was slower,
so no thread speedup is claimed. The post-spline Prüfer comparison remains
`3.137x` on four modes with max spectrum error `5.22e-7 dex`; raw data are in
[`reference_prufer_production_ab_post_spline_20261104.json`](reference_prufer_production_ab_post_spline_20261104.json)
and [`reference_spline_workers_20261104.json`](reference_spline_workers_20261104.json).
Cold dispatch records are in
[`reference_spline_cold_low_20261104.json`](reference_spline_cold_low_20261104.json)
and [`reference_spline_cold_high_20261104.json`](reference_spline_cold_high_20261104.json);
reproduce them with [`benchmark_reference_spline_cold.py`](../scripts/benchmark_reference_spline_cold.py).
On the same updated source, a candidate-only 260-mode tail10 Prüfer run at the
independent Reference root had map difference `4.39e-7` relative and maximum
spectrum error `5.25e-5 dex`; runtime was `569.1 s` on one worker. It is not a
paired full-grid Cartesian speed comparison. See
[`reference_prufer_spline_res4_refroot_20261104.json`](reference_prufer_spline_res4_refroot_20261104.json).
The expanded eight-frequency A/B reports `3.626x` aggregate (`3.586x`
including JIT warmup), per-mode P50/P95 `2.778x/4.174x`, and maximum same-input
spectrum error `3.94e-6 dex`; see
[`reference_prufer_production_ab_8mode_20261104.json`](reference_prufer_production_ab_8mode_20261104.json).
The complete 260-mode alternating A/B at the independent res4 root is now
available in
[`reference_prufer_cartesian_paired_fullgrid_20261104.json`](reference_prufer_cartesian_paired_fullgrid_20261104.json):
all modes completed, paired spectrum error median/P95/max is
`1.575e-7/6.654e-6/5.246e-5 dex`, aggregate speedup is `3.470x` (`3.468x`
including JIT), and per-mode speedup P50/P95 is `2.632x/4.332x`. Reconstruct
integral maps and the secant-based root-shift estimate with
[`analyze_reference_paired_fullgrid.py`](../scripts/analyze_reference_paired_fullgrid.py);
the secant estimate is not a rigorous root bracket.

The standalone Newton–Kantorovich phase-transfer screen is also recorded at
[`nonosc_phase_res4_mode_20261104.json`](nonosc_phase_res4_mode_20261104.json).
Three res4 modes cover two tail handoffs and one no-tail propagation across the
`nu=.1` and `nu=100` background boundaries. At order 24 they match current
Cartesian within `1.7461e-8`, `4.1274e-9`, and `1.1137e-7 dex`; single-run
candidate times are `0.138/0.137/0.499 s` versus `10.935/6.457/9.520 s`.
These are exploratory single-mode timings, not distributional benchmarks. The
method has not passed full-grid, integral-map, self-consistency, or
parameter-domain tests and remains a research prototype. Reproduce with
[`prototype_nonosc_phase_res4_mode.py`](../scripts/prototype_nonosc_phase_res4_mode.py)
and `--frequency -10.772906114757099` or
`--frequency -14.6013 --phase-start 6`.
The subsequent complete-grid screen completed 213/260 modes (maximum/P95
order-32 error `1.548e-4/5.562e-8 dex`); 45 modes failed the high-frequency
screen and two had non-positive `Q`. The existing 20-step Newton limit was
retained; 14 modes passed the collocation and dense residual gates despite not
reaching the tighter update stopping test. The `50.8x` summed time ratio
applies only to the 213 supported modes and compares a single candidate run
with archived Cartesian timings; it is not full-grid coverage or a paired
performance claim.
See [`nonosc_phase_res4_fullgrid_20261104.json`](nonosc_phase_res4_fullgrid_20261104.json)
and reproduce with
[`benchmark_nonosc_phase_res4_fullgrid.py`](../scripts/benchmark_nonosc_phase_res4_fullgrid.py).
An isolated `--minimum-gamma 0` full-grid screen completed 219/260 modes, with
39 residual failures and two non-positive-`Q` failures; it recovered six
additional modes, so the default screen remains unchanged. Its artifact is
[`nonosc_phase_res4_fullgrid_gamma0_20261104.json`](nonosc_phase_res4_fullgrid_gamma0_20261104.json).
An independent order probe of its 39 residual failures found 5 no-tail modes
that converge between orders 16/24 and match Cartesian within `3.28e-9 dex`;
order 32 failed for all 39. These are recorded as two-order cross-checks, not
as part of the stricter 219-mode 16/24/32 count. See
[`nonosc_phase_res4_residual_order_probe_20261104.json`](nonosc_phase_res4_residual_order_probe_20261104.json)
and rerun with
[`analyze_phase_residual_order_failures.py`](../scripts/analyze_phase_residual_order_failures.py).
Three additional single-mode checks at `DN_eff=0.133732083004284` match
Cartesian within `8.42e-9`, `4.12e-9`, and `3.39e-9 dex` for two tail modes and
one no-tail mode. They vary the background normalization only; see
[`nonosc_phase_res4_dn133732_tail_cross_20261104.json`](nonosc_phase_res4_dn133732_tail_cross_20261104.json),
[`nonosc_phase_res4_dn133732_high_20261104.json`](nonosc_phase_res4_dn133732_high_20261104.json),
and
[`nonosc_phase_res4_dn133732_notail_20261104.json`](nonosc_phase_res4_dn133732_notail_20261104.json).
At a distinct minimized Fast pressure point (`r=.19649`, `n_t=1.06225`,
`T_re=1761.31 GeV`, `kappa10=.0115985`, fixed `DN_eff=0`), six tail/no-tail
modes match Cartesian within `3.86e-8 dex` maximum; this is not full-spectrum
or root validation.
See [`nonosc_phase_sample2_tail_cross_20261104.json`](nonosc_phase_sample2_tail_cross_20261104.json),
[`nonosc_phase_sample2_high_20261104.json`](nonosc_phase_sample2_high_20261104.json),
[`nonosc_phase_sample2_notail_20261104.json`](nonosc_phase_sample2_notail_20261104.json),
[`nonosc_phase_sample2_mid1_20261104.json`](nonosc_phase_sample2_mid1_20261104.json),
[`nonosc_phase_sample2_mid2_20261104.json`](nonosc_phase_sample2_mid2_20261104.json),
and
[`nonosc_phase_sample2_mid3_20261104.json`](nonosc_phase_sample2_mid3_20261104.json).
At the sample2 `freq_res=8`, `z_tail=5` Reference root (`DN_eff=2.2914118670`),
a fresh fixed Cartesian map over 800 frequencies closes to `-1.85e-10` relative
root residual. The initial phase screen matched 700/800 at `gamma_min=10` and
709/800 at0, with maximum error `8.668e-7 dex` on the strict three-order set,
but took about `145 s` versus `45 s` for Cartesian. It was not a speed win and
the 91 unsupported modes prevented a complete phase map. The initial
residual-order replay added 14 adjacent-order plus Cartesian checks (723 modes
at that stage); the follow-up results are recorded below. See
[`sample2_phase_reference_root_res8_tail5_20261104.json`](sample2_phase_reference_root_res8_tail5_20261104.json),
[`sample2_phase_reference_root_res8_tail5_gamma0_20261105.json`](sample2_phase_reference_root_res8_tail5_gamma0_20261105.json),
and
[`sample2_phase_root_residual_order_probe_20261105.json`](sample2_phase_root_residual_order_probe_20261105.json).
Two modes also exercise the positive `cr` branch; at `cr=.1`, `r=.01` the
model derives `n_t=-.00125`, and the tail/no-tail errors are `3.54e-8` and
`3.21e-9 dex`. See
[`nonosc_phase_cr_positive_tail_20261104.json`](nonosc_phase_cr_positive_tail_20261104.json)
and
[`nonosc_phase_cr_positive_notail_20261104.json`](nonosc_phase_cr_positive_notail_20261104.json).
On the distinct minimized sample2 cosmology's full 800-node Reference-root
grid, phase propagation passes 700/800 at `gamma_min=10` and737/800 at0, with
maximum error `8.668e-7 dex` on either strict set; gamma0 takes132.2s versus
43.0s for Cartesian and has no complete phase integral map. A tail-event
classification regression corrected four `z_today<z_tail` cases. Stopping
Newton at the unchanged collocation gate and conditional backtracking add nine
strict modes; least-squares fallback for steps that still fail to decrease the
residual adds19 more, `[656,671,672,687,700,701,702,703,714,716,718,726,730,
733,734,755,756,767,769]`. No prior strict mode is lost. Among57 remaining
residual failures,20 pass adjacent-order checks and70 completed order
candidates across44 modes pass tighter Cartesian/Prüfer references. All20
adjacent-order modes overlap that set, so781/800 are validated at the standard
handoff. A separate phase-start scan with 1/128 spacing dual-validates six more
modes `[757,759,780,782,783,786]` against tighter Cartesian/Prüfer references,
for787/800 unique frequencies; seven residual failures and six non-positive-`Q`
modes remain; refining from 1/64 to 1/128 adds no further modes. Lower
orders 8/12/16 add no unique coverage.
These checks still do not provide a complete phase integral map or root.
The Windows phase archive records mode687 as a pass, while Linux CI reproduces
an explicit order-16 residual rejection (`7.76e-10` vs the `1e-10` gate). This
candidate is platform-sensitive and the `786/800` count is Windows-specific;
no cross-platform count is claimed. The regression accepts either a result
that meets the original gates or the explicit rejection, with no threshold
change. See
[`sample2_phase_reference_root_res8_tail5_20261104.json`](sample2_phase_reference_root_res8_tail5_20261104.json),
[`sample2_phase_reference_root_res8_tail5_gamma0_20261105.json`](sample2_phase_reference_root_res8_tail5_gamma0_20261105.json),
[`sample2_phase_root_residual_order_probe_20261105.json`](sample2_phase_root_residual_order_probe_20261105.json),
[`sample2_phase_residual_orders_dual_reference_20261105.json`](sample2_phase_residual_orders_dual_reference_20261105.json),
[`sample2_phase_root_lower_order_probe_20261105.json`](sample2_phase_root_lower_order_probe_20261105.json),
[`sample2_phase_start_offset_probe_20261105.json`](sample2_phase_start_offset_probe_20261105.json),
[`sample2_phase_reference_root_res8_tail5_gamma10_20261105.json`](sample2_phase_reference_root_res8_tail5_gamma10_20261105.json),
and [`sample2_reference_root_res8_tail5_gamma0_20261105.json`](sample2_reference_root_res8_tail5_gamma0_20261105.json).
The lower-order probe completed0/0/15 of the latest57 residual modes at
orders8/12/16 and found no adjacent-order pair; keep it as an exploratory
result only. See
[`sample2_phase_root_lower_order_probe_20261105.json`](sample2_phase_root_lower_order_probe_20261105.json)
and `scripts/probe_sample2_lower_collocation_orders.py`.
Splitting each remaining residual interval into two or four subsegments produced
no screened candidate across the seven remaining modes and is not retained;
raw evidence is in
[`sample2_phase_residual_segment_split_probe_20261105.json`](sample2_phase_residual_segment_split_probe_20261105.json)
and `scripts/probe_sample2_residual_segment_splits.py`.
The targeted start-offset replay is `scripts/probe_sample2_phase_start_offsets.py`.
Two archived Cartesian tail10 maps also give a same-source, same-grid secant
slope estimate of `0.02113` and a first-order root shift estimate of
`-1.08e-10`. This is not a local derivative or root enclosure; see
[`res4_reference_map_sensitivity_20261104.json`](res4_reference_map_sensitivity_20261104.json)
and [`analyze_reference_map_sensitivity.py`](../scripts/analyze_reference_map_sensitivity.py).

## 2026-11-05 Reference cost experiments

The Reference full-grid profile identified continuous `H2`/`sigma` evaluation
as the main Python-level cost: 4.32M DOP853 RHS calls and 19.55s inclusive
background time in a 45.16s sample2 run. A Numba implementation of the same
piecewise equations and existing cubic coefficients matched 2,017 branch and
random points bitwise. On the exact 800-node sample2 fixed background, three
paired full-grid A/B runs matched the archived spectrum, tail mask and
`DN_gw` bitwise while reducing warm P50 runtime from 41.392s to 27.513s (33.53%).
The 0.579s JIT warm cost is separate; the candidate is batch-warmed only for
32 or more frequencies that include costly modes, and single-frequency requests keep
the no-new-JIT path. This is one cosmology/configuration only. Reproducible
records: [`reference_h2_numba_fullgrid_ab_20261105.json`](reference_h2_numba_fullgrid_ab_20261105.json),
[`reference_fullgrid_profile_20261105.json`](reference_fullgrid_profile_20261105.json),
[`reference_fullgrid_profile_compiled_20261105.json`](reference_fullgrid_profile_compiled_20261105.json),
[`reference_h2_numba_probe_20261105.json`](reference_h2_numba_probe_20261105.json),
and driver [`benchmark_reference_h2_numba_fullgrid.py`](../scripts/benchmark_reference_h2_numba_fullgrid.py).

The batch warmup crossover was also checked at 16, 32 and 64 stratified nodes.
At 32 nodes, the cold compiled call including 0.280s H2 compile was 0.272s faster
than the warm Python P50; all outputs matched bitwise. The 16-node path did not
compile H2. See [`reference_h2_numba_batch_crossover_20261105.json`](reference_h2_numba_batch_crossover_20261105.json)
and [`benchmark_reference_h2_batch_crossover.py`](../scripts/benchmark_reference_h2_batch_crossover.py).

The `_find_start_N` anchor-reuse change produced bitwise-identical start
coordinates across the full 800-node sample2 grid and bitwise-identical complete
outputs on five modes; paired with the compiled H2 evaluator, full-map warm P50
improved 0.88% (27.853s→27.608s). The result is recorded separately in
[`reference_start_search_fullgrid_ab_20261105.json`](reference_start_search_fullgrid_ab_20261105.json)
and [`benchmark_reference_start_search_fullgrid.py`](../scripts/benchmark_reference_start_search_fullgrid.py).
The five-frequency Prüfer ODE screen remains mixed at0.86–1.57x and was not
promoted to a full-grid speed claim.

The fused Numba tensor RHS is enabled only after a costly Reference batch of
at least32 frequencies. Three exact 800-frequency A/B pairs match the archived
Cartesian spectrum and `DN_gw` bitwise; warm P50/P95 is `27.258/27.361 s` for
the current compiled-H2 RHS and `23.851/23.959 s` for the fused RHS (12.5%
lower P50). The measured kernel/background warm cost is `0.679 s`. A 32-node
stratified screen estimates the first compiled batch at `1.036 s`, below the
current-path P50 `1.098 s`; outputs match bitwise. This is one cosmology and
one worker, so it is not a multi-background throughput or memory claim. Drivers
and records:
[`benchmark_reference_numba_rhs_fullgrid.py`](../scripts/benchmark_reference_numba_rhs_fullgrid.py),
[`reference_numba_rhs_fullgrid_ab_20261105.json`](reference_numba_rhs_fullgrid_ab_20261105.json),
[`benchmark_reference_numba_rhs_batch_crossover.py`](../scripts/benchmark_reference_numba_rhs_batch_crossover.py),
and [`reference_numba_rhs_batch_crossover_20261105.json`](reference_numba_rhs_batch_crossover_20261105.json).

A positive-`cr` cross-check extends the compiled-vs-Python exact-equality
comparison to two modes and compares Cartesian with Prüfer propagation. Both
formulations agree within `3.80e-8` relative in `Omega_GW` on the selected
frequencies; the sample is not a full-grid or root result. Driver and raw
record: [`validate_reference_h2_numba_cr_positive.py`](../scripts/validate_reference_h2_numba_cr_positive.py)
and [`reference_h2_numba_cr_positive_20261105.json`](reference_h2_numba_cr_positive_20261105.json).

## 2026-11-05 broad Fast parameter pilot

Two512-case combined Sobol/LHS pilots returned13 finite-positive success
outputs and4 nonfinite outputs; all13 successes remain uncertified. Most
samples were rejected by the resource preflight or ended in a failure or
exception. This is a sampling screen, not a measured success rate or a global
boundary. See the raw files and statuses listed in
[`fast_domain_pressure_20261030.md`](fast_domain_pressure_20261030.md#2026-11-05-256-point-sobol-and-latin-hypercube-pilots).

## 2026-11-05 subnormal integral scaling and extreme frequency-grid audit

At QMC sample146, the fixed-background Prüfer Reference spectrum reached
subnormal values. SciPy PCHIP warned while forming reciprocal secant slopes.
`integrate_spectrum` now scales only this small-signal case and restores the
same scale to the integral and error estimates. A warning-as-error regression
and the focused Reference suite pass (`57 passed, 1 deselected`).

With `rtol=1e-9`, `z_tail=5`, and one worker, full Prüfer maps at frequency
resolutions 1/2/4/8 use 295/590/1181/2362 nodes. Res=4 to res=8 `DN_gw` changes
by `6.48e-5` relative, but res=1 to res=8 changes by `1.066e-2`. The native
77-node Fast integral at this extreme point is therefore not certified. This
is a fixed-background single-point convergence check, not a parameter-domain
claim. See [accuracy evidence](accuracy.md#2026-11-05-subnormal-reference-integral-and-extreme-grid-convergence),
the [driver](../scripts/validate_qmc_extreme_reference_grid_convergence.py),
and the [full result](qmc_sobol146_reference_grid_convergence_20261105.json).
On the same native 77 frequencies, Cartesian versus Prüfer tail5 differs by at
most `2.89e-7 dex`, while default Fast versus tail10 differs by `0.00645 dex`
and `1.47%` in `DN_gw`. `sigma_exact` plus tail10 reduces these same-grid
differences to `0.000346 dex` and `5.09e-5` relative, but does not establish a
resolved bolometric integral. Per-mode and tail-depth evidence is stored in
[`qmc_sobol146_reference_crosscheck_20261105.json`](qmc_sobol146_reference_crosscheck_20261105.json),
[`qmc_sobol146_full_prufer_tail_convergence_20261105.json`](qmc_sobol146_full_prufer_tail_convergence_20261105.json),
[`qmc_sobol146_mode_tail_probe_20261105.json`](qmc_sobol146_mode_tail_probe_20261105.json),
and [`qmc_sobol146_prufer_tail_depths_20261105.json`](qmc_sobol146_prufer_tail_depths_20261105.json).
The separate 800-node normal-scale integrator A/B returned bitwise-identical
integrals and error estimates, with P50 `0.079655 s` before and `0.079686 s`
after (0.04% measured overhead); see
[`reference_integrate_spectrum_scale_ab_20261105.json`](reference_integrate_spectrum_scale_ab_20261105.json)
and [`benchmark_reference_integrate_spectrum_scale.py`](../scripts/benchmark_reference_integrate_spectrum_scale.py).
Reintegrating the archived native 77-mode tail5/tail10 Prüfer outputs with the
current code changed `DN_gw` by only `1.74e-16` relative and produced no
warnings; see
[`qmc_sobol146_post_scale_reintegration_20261105.json`](qmc_sobol146_post_scale_reintegration_20261105.json)
and [`reintegrate_qmc_sobol146_post_scale.py`](../scripts/reintegrate_qmc_sobol146_post_scale.py).

The Fast custom PCHIP slope now uses a scaled weighted harmonic-mean formula
that avoids reciprocal overflow on subnormal tails. The new warning-as-error
regression passes. A 30-pair normal Fast A/B is bitwise identical with no warm
P50/P95 regression. The high-level Fast API gained optional adaptive-grid
arguments without changing its default. On sample146, the default 77-node goal
result differs from res=8 by `84.6%` relative to Reference (`45.8%` with the
symmetric max-value denominator). Adaptive target `3e-4` still leaves a
`2.58%` relative map difference; the explicitly validated
`adaptive`, `target=1e-6`, `sigma_exact`, `z_tail=10` route used 321 nodes and
passed an exact-array Prüfer comparison (`5.06e-5` relative `DN_gw`, maximum
`3.72e-4 dex`). It took 660.6s on one worker and applies to this single point.
See [accuracy evidence](accuracy.md#2026-11-05-stable-fast-pchip-slopes-and-adaptive-grid-point-check),
the [sweep](qmc_sobol146_fast_adaptive_grid_20261105.json),
the [same-grid Reference result](qmc_sobol146_fast_adaptive_reference_20261105.json),
and [normal-path benchmark](fast_pchip_stable_slopes_ab_20261105.json).
The public wrapper reproduces the same frequency hash and `DN_gw`; its
continuity check and spectrum-copy limitation are recorded in
[`qmc_sobol146_adaptive_api_continuity_20261105.json`](qmc_sobol146_adaptive_api_continuity_20261105.json).

## 2026-11-05 positive-cr adaptive Fast fixed-background validation

The opt-in adaptive Fast grid was evaluated on a standard-parameter
`cr=1` background (effective `n_t=-0.00125`) at target `1e-4`. All 256
positive spectrum bins agree with same-array Cartesian and Prüfer Reference
within `5.3041e-4 dex`; the fixed-map `DN_gw` differences are below `1.0e-5`
relative. A stratified Cartesian/Prüfer subset agrees within `2.14e-5 dex`.
Full Fast and Prüfer raw arrays and environment/source hashes are in
[`adaptive_cr1_prufer_fullgrid_20261105.json`](adaptive_cr1_prufer_fullgrid_20261105.json),
driven by
[`validate_adaptive_cr1_prufer_fullgrid.py`](../scripts/validate_adaptive_cr1_prufer_fullgrid.py).

The 484-node native Prüfer `freq_res=2` map differs by `9.9973e-5` relative in
`DN_gw` from the exact 256-node adaptive-grid Prüfer map, and the Fast output
differs from the dense map by `1.0989e-4`; both are within the `2e-4` gate.
Full native-grid arrays are in
[`adaptive_cr1_prufer_res2_20261105.json`](adaptive_cr1_prufer_res2_20261105.json),
driven by
[`validate_adaptive_cr1_prufer_res2.py`](../scripts/validate_adaptive_cr1_prufer_res2.py).
This is a two-grid convergence screen, not a formal error bound or self-consistent
Reference root certification. The higher Prüfer `freq_res=4` and alternate
Cartesian `freq_res=2` maps were cancelled at the 1200 s CPU budget without
results; these are resource timeouts, not numerical failures. The Fast target
sweep is also non-monotone. See the scoped partial record
[`adaptive_standard_background_validation_20261105.json`](adaptive_standard_background_validation_20261105.json)
and the detailed [accuracy note](accuracy.md#2026-11-05-adaptive-fast-positive-cr-fixed-background-cross-check).
Reference workers=4 were slower than workers=1 (132.56 s versus 84.36 s on
the 32-bin subset) with bitwise-identical output; no parallel acceleration is
claimed.

The follow-up Fast fixed-point screen at `N_inf=122.502` evaluates native-grid
Reference maps around the Fast self-consistent value. Native and adaptive-grid
interval secant roots differ by `8.95e-5` relative, and the native secant
estimate differs from Fast by `1.15e-4`. This passes the `2e-4` pointwise root
screen, but is not a formal bound or a full iterative Reference root solve; the
two grids do not share a common narrow sign-changing bracket. See
[`fast_ninf122_reference_root_secant_20261105.json`](fast_ninf122_reference_root_secant_20261105.json)
and [`probe_fast_ninf122_reference_root_secant.py`](../scripts/probe_fast_ninf122_reference_root_secant.py).

## 2026-11-05 Fast tail-depth repair above `N_inf=120`

A standard-positive-density Fast point at `N_inf=122.502` was skipped by the
old pressure harness's `N_inf>120` resource cutoff but runs directly. Default
75-node and adaptive target-`1e-4` 242-node Fast maps differ by about `1.11%`;
both are finite but initially lacked accuracy certification. Same-frequency
Cartesian and Prüfer agree, locating the `z_tail=5` spectral error in Fast's
tail handoff: the selected-bin maximum is `2.89e-3 dex`; raising `z_tail` to8
reduces it to `1.46e-4 dex`. The full tail8/target-`1e-4` map passed its
same-grid spectrum screen but failed the dense res2 DN screen at `3.84e-4`.

The tail8/target-`1e-5` Fast configuration uses299 nodes and passes the
same-grid Prufer spectrum check (`1.656e-4 dex` max), the `DN_gw` map screen
(`2.50e-5` relative), and the 451-node native res2 screen (`8.95e-5` grid
change; `1.145e-4` Fast-to-dense map difference). This is one fixed-background
Fast configuration, not a Reference root or a global dispatch rule. Full raw
data and reproduction files:
[`fast_ninf122_tail8_target1e5_validation_20261105.json`](fast_ninf122_tail8_target1e5_validation_20261105.json),
[`validate_fast_ninf122_tail8_target1e5.py`](../scripts/validate_fast_ninf122_tail8_target1e5.py),
[`fast_ninf122_tail_sensitivity_20261105.json`](fast_ninf122_tail_sensitivity_20261105.json),
and [`probe_fast_ninf122_tail_sensitivity.py`](../scripts/probe_fast_ninf122_tail_sensitivity.py).

With `DN_re=200`, the same Fast settings reach `N_inf=222.502`. The full
299-node same-grid Prüfer spectrum passes at `1.4962e-4 dex` maximum; the
451-node native res2 map changes by `8.95e-5` relative, and Fast differs from
the native map by `1.1453e-4`. This is a fixed-background validation only; the
independent self-consistent root at this second point remains unchecked. See
[`fast_ninf222_fullgrid_prufer_20261105.json`](fast_ninf222_fullgrid_prufer_20261105.json)
and [`validate_fast_ninf222_fullgrid_prufer.py`](../scripts/validate_fast_ninf222_fullgrid_prufer.py).

The low-temperature positive-cr point (`T_re=1e-3`) has derived `N_inf=80.4482`,
`n_t=-0.00125`, and effective `DN_re=57.9462`. Fast passes the full spectrum
screen (`1.5079e-4 dex` max); native res2 narrowly misses the DN gate, while
res4/res8 converge within `8.971e-5` and Fast-to-res8 is `1.743e-5`. This is a
fixed-background result; the Reference self-consistent root remains open. See
[`fast_cr1_lowT_fullgrid_prufer_20261105.json`](fast_cr1_lowT_fullgrid_prufer_20261105.json),
[`fast_cr1_lowT_res4_prufer_20261105.json`](fast_cr1_lowT_res4_prufer_20261105.json),
[`fast_cr1_lowT_res8_prufer_20261105.json`](fast_cr1_lowT_res8_prufer_20261105.json),
and [`validate_fast_cr1_lowT_fullgrid.py`](../scripts/validate_fast_cr1_lowT_fullgrid.py).

The follow-up fixed-point screen estimates adaptive and native res4 secant
roots that differ by `7.258e-5`; the Fast value differs from the res4 estimate
by `7.229e-5`. An adaptive map at the res4 candidate gives a corrected
relative residual `7.258e-5`. The two grids do not share a sign-changing
bracket, so this is a two-grid secant screen, not a formal root enclosure. See
[`fast_cr1_lowT_root_secant_20261105.json`](fast_cr1_lowT_root_secant_20261105.json)
and [`probe_fast_cr1_lowT_root_secant.py`](../scripts/probe_fast_cr1_lowT_root_secant.py).

The high-blue fixed-background point `cr=0, r=10, n_t=2, kappa10=.01`,
`T_re=.001`, `DN_re=200` passes the full spectrum and res4 integral screens
only at the tighter adaptive target `1e-7`: maximum spectrum error is
`1.464e-4 dex`, Fast-to-res4 DN difference `6.416e-5`. Targets `1e-5` and
`1e-6` miss the Fast-to-res4 integral gate despite passing the same-grid
spectrum check. The cumulative panel indicator is recorded as a diagnostic,
not a certified error bound. See
[`fast_high_blue_fullgrid_prufer_20261105.json`](fast_high_blue_fullgrid_prufer_20261105.json),
[`fast_high_blue_target1e6_prufer_20261105.json`](fast_high_blue_target1e6_prufer_20261105.json),
[`fast_high_blue_target1e7_prufer_20261105.json`](fast_high_blue_target1e7_prufer_20261105.json),
and the three corresponding validation scripts in `scripts/`.

The target-`1e-7` high-blue point also passes a two-grid secant root screen:
adaptive and native-res4 estimates differ by `8.110e-5`, Fast differs from
the res4 estimate by `6.553e-5`, and the corrected adaptive residual at the
res4 candidate is `8.114e-5`. This is not a formal root enclosure; the sampled
grids do not share a sign-changing bracket. See
[`fast_high_blue_root_secant_20261105.json`](fast_high_blue_root_secant_20261105.json)
and [`probe_fast_high_blue_root_secant.py`](../scripts/probe_fast_high_blue_root_secant.py).

The joint `r=10,n_t=2,kappa10=100` point at `N_inf=222.502` also passes the
target-`1e-7` fixed-background spectrum and res2/res4 map screens: maximum
spectrum error `1.472e-4 dex`, res2→res4 difference `1.717e-4`, and
Fast→res4 difference `6.248e-5`. Its self-consistent root remains open. See
[`fast_high_blue_kappa100_prufer_20261105.json`](fast_high_blue_kappa100_prufer_20261105.json)
and [`validate_fast_high_blue_kappa100.py`](../scripts/validate_fast_high_blue_kappa100.py).

Its two-grid secant root estimates differ by `8.098e-5`; the corrected
adaptive residual at the native res4 candidate is `8.099e-5`. The adaptive
interval brackets the residual sign change, but native res4 does not, so the
result is a screening check rather than a formal root enclosure. See
[`fast_high_blue_kappa100_root_secant_20261105.json`](fast_high_blue_kappa100_root_secant_20261105.json)
and [`probe_fast_high_blue_kappa100_root.py`](../scripts/probe_fast_high_blue_kappa100_root.py).

The high-r strongly red positive-cr point (`r=10`, derived `n_t=-1.25`,
`kappa10=100`, `N_inf=82.6947`) passes fixed-background spectrum and res4/res8
map checks at target `1e-7`: maximum spectral error `1.476e-4 dex`, res4→res8
`6.437e-5`, and Fast→res8 `4.526e-6`. Native res2 is under-resolved by
`9.663e-4`; the self-consistent Reference root was not yet checked in this
fixed-background report. See
[`fast_cr1_red_extreme_fullgrid_prufer_20261105.json`](fast_cr1_red_extreme_fullgrid_prufer_20261105.json),
[`fast_cr1_red_extreme_res8_prufer_20261105.json`](fast_cr1_red_extreme_res8_prufer_20261105.json),
and [`validate_fast_cr1_red_extreme_fullgrid.py`](../scripts/validate_fast_cr1_red_extreme_fullgrid.py).

Its two-grid secant root estimates differ by `6.712e-5`; Fast-to-res4-root
shift is `6.890e-5` and the corrected adaptive residual at the res4 candidate
is `6.713e-5`. Neither grid shows a sampled endpoint sign change, and the
native res4 frequency hash shifts with `DN_eff`; this is a scalar-map secant
screen, not a root enclosure. See
[`fast_cr1_red_extreme_root_secant_20261105.json`](fast_cr1_red_extreme_root_secant_20261105.json)
and [`probe_fast_cr1_red_extreme_root.py`](../scripts/probe_fast_cr1_red_extreme_root.py).

The subsequent self-consistent Prüfer Reference solve on the same validated
743-frequency Fast grid converged in 3 iterations. Its root is
`DN_eff=0.0017772220202209284` (fixed-point residual `-7.35e-14`); Fast is
`1.779e-6` relative away. The solve took `331.16 s`, one worker, with no
warnings. This closes the fixed-grid root check only; res4/res8 convergence at
the iterative root was subsequently checked with a direct native res4 map:
`DN_gw=0.0017771027289385473`, relative root residual `6.712e-5`, and only
`1.47e-9` relative change from the x0 res4 map. The 1,572-point res8 map is
reused at the nearly identical Fast root; no res8 propagation was repeated at
the iterative root. See
[`fast_cr1_red_extreme_reference_root_20261105.json`](fast_cr1_red_extreme_reference_root_20261105.json)
and [`solve_fast_cr1_red_extreme_reference_root.py`](../scripts/solve_fast_cr1_red_extreme_reference_root.py), plus
[`fast_cr1_red_extreme_root_res4_20261105.json`](fast_cr1_red_extreme_root_res4_20261105.json)
and [`validate_fast_cr1_red_extreme_root_res4.py`](../scripts/validate_fast_cr1_red_extreme_root_res4.py).
