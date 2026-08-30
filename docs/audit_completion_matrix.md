# Audit completion matrix (2026-08-30)

This matrix separates implementation evidence from scientific certification. A
row marked **open** is intentionally not promoted to a success claim.

| Requirement | Evidence | Status |
|---|---|---|
| Distinguish engine/grid/physical/likelihood/posterior errors | `audit_phase1.md`, `audit_phase2.md`, `audit_phase3.md`, `audit_mcmc.md` | **done** |
| Canonical Cobaya derived names and contract | `stiffgwpy/cobaya/stiffGW.py`, `stiffGW.yaml`, `tests/test_cobaya_adapter.py` | **done** |
| Qualified import after install | `tests/test_cobaya_adapter.py`, `scripts/_smoke_probe.py` (Cobaya 3.6.2) | **done** |
| Fast/fallback observability | `LCDM_SG` telemetry fields (including guard/input-vs-numerical reason), adapter `engine_stats`, Cobaya `close()` log/warning | **done** |
| Production defaults | `stiffGW.yaml`: production, h=0.01, col_step=4, z_tail=7, freq_res=1, 8 threads; regression-tested in `test_cobaya_adapter.py` | **done** |
| 1030-point parameter sweep | `docs/paramsweep/sweep_phase3.jsonl`, `summary.json`, `worst_cases.json`, `failure_map.json` | **done** |
| Cache/state isolation | `tests/test_engine.py` A→B→A and thread 1→8→1 regression | **done** |
| Cold/warm benchmark and fallback fraction | `scripts/bench_fast.py`, `docs/audit_summary.md` | **done** |
| Cobaya posterior equivalence with >=2000 ESS | `scripts/mcmc_compare.py` has ESS/covariance/distance gates, records YAML SHA-256 and an explicit same-prior/seed/sampler contract, and supports a deterministic shared initial point; no qualifying chains yet | **open** |
| Final posterior shift decision (<0.1 sigma) | `mcmc_compare.py` now enforces a per-parameter `max_abs_posterior_shift` gate (default 0.1); requires the open long-chain run | **open** |

## Current quantitative evidence

- Same-grid engine difference: approximately 1e-5 in final Delta N_eff at the
  production-scale grid; this is not a continuum accuracy bound.
- Shared-grid/continuum-reference bias: approximately 0.7% at h=0.01 in the
  documented convergence study.
- Full sweep (latest record per point): 782 `ok` and 248 shared N_eff guard
  rejections. Three transient LSODA `MemoryError` records occurred during an
  overloaded initial attempt and were reproduced successfully in a serial
  retry; they remain an infrastructure-risk note rather than final failures.
- Existing short-chain comparison: max absolute pointwise Delta logL about
  0.090, but N=20 chains have minimum ESS about 1 and are not posterior
  certification evidence.
- A 10,000-accepted fast-only pilot (production, seed 20260830; compact record in
  `docs/mcmc/fast_pilot_10000_summary.json`) reached minimum
  ESS 4.44; it intentionally skipped LSODA and pointwise checks, so it remains
  diagnostic only. Telemetry: 30014 fast evaluations, 16 numerical failures,
  2460 guard rejections, and zero fallbacks with `--no-fallback`.
- The comparison runner now accepts `--initial-point`; the baseline point is
  stored in `docs/mcmc/initial_point_baseline.json` and is applied identically
  to both engine configurations without changing their priors.
