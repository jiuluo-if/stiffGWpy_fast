# Parameter-space validation, matched z8 (fast vs continuous-sigma reference)

Date: 2026-09-03. Engine truth = `stiffgwpy/reference.py` (continuous-sigma + DOP853).
LSODA is not a truth anchor anywhere in this document. Every point is solved by BOTH
engines on the SAME `grid_independent` frequency grid and the SAME `z_tail=8.0`
(242-249 nodes), so the residual isolates engine error (no grid/tail mismatch).

Artifacts:
- `docs/paramsweep_z8/reference_points.jsonl` (9 matched single points)
- `docs/paramsweep_z8/validation_summary.json` (aggregate gates)
- `docs/paramsweep_ref/fast_sweep.jsonl` (240-point Sobol production sweep, z8)

## 1. Matched single points (Layer A physics validation, 9/9 completed)

fast: production z8 (phase_max sub-stepping, transition_refine, grid_independent f_res=1).
ref: DOP853 rtol=1e-9, continuous sigma, z_tail=8. rel = (DN_fast-DN_ref)/DN_ref.
dex/rel are per-mode statistics in the signal band log10 f in [-6, +1].

| point | DN_fast | DN_ref | DN rel | signal dex max | signal rel max | fast s | ref s |
|---|---|---|---:|---:|---:|---:|---:|
| default | 0.002263514 | 0.002262833 | +3.01e-4 | 2.60e-4 | 5.98e-4 | 4.1 | 359.4 |
| stiff (kappa10=1) | 0.225529102 | 0.225626633 | -4.32e-4 | 2.93e-4 | 6.75e-4 | 2.6 | 358.4 |
| lowT (T_re=10) | 5.2e-8 | 5.2e-8 | +1.46e-3* | 2.61e-4 | 6.02e-4 | 1.6 | 368.3 |
| highT (T_re=1e4) | 0.056413060 | 0.056405888 | +1.27e-4 | 3.00e-4 | 6.92e-4 | 3.5 | 362.0 |
| rad_dominant (kappa10=1e-6) | 2.308e-7 | 2.306e-7 | +7.18e-4 | 2.79e-4 | 6.43e-4 | 1.9 | 375.0 |
| tiny_r (r=1e-6) | 2.3692e-7 | 2.3691e-7 | +3.82e-5 | 3.08e-4 | 7.09e-4 | 1.2 | 349.7 |
| transition-sensitive | 0.000723010 | 0.000723557 | -7.56e-4 | 2.88e-4 | 6.64e-4 | 3.4 | 357.3 |
| cr0_blue (cr=0, n_t=0.2, DN_re=5) | 0.057022070 | 0.056987328 | +6.10e-4 | 2.85e-4 | 6.58e-4 | 1.5 | 253.2 |
| extreme (r=3e-2, kappa10=1e-1) | 1.524201531 | 1.524863254 | -4.34e-4 | 3.00e-4 | 6.92e-4 | 3.4 | 361.7 |

* lowT DN ~ 5.2e-8 is a DN-of-DN amplification artifact: both engines agree to
  5.5e-10 absolute; the relative number has a tiny denominator. The physically
  meaningful per-mode signal error (rel max 6.02e-4, dex 2.6e-4) passes the gate.

Aggregate (`validation_summary.json`, n=9):
- signal band: rel max 7.09e-4, rel p95 7.02e-4, dex max 3.08e-4 -> PASS (<1e-3)
- transition band (-2..0): rel max 7.09e-4, rel p95 7.02e-4 -> PASS (<1e-3)
- integrated DN rel abs: median 4.34e-4, p95 1.18e-3, max 1.46e-3 (lowT artifact)
  -> NOT <1e-4. The engine architectural residual is ~3e-4..7.6e-4 relative on
  Delta_Neff at z8 matched settings (median 4.3e-4); 1e-4 is below the current
  per-mode Magnus/grid architecture residual. See honest limits below.

## 2. 240-point Sobol sweep (Layer B), production fast engine, z8

`docs/paramsweep_ref/fast_sweep.jsonl` -- Sobol (scrambled, seed 20260831) over
r in [1e-6,1e-1] log, n_t in [-0.5,0.5], cr in {0,1}, T_re in [10,1e6] log,
DN_re in [0,30], kappa10 in [1e-6,1] log (LCDM anchors fixed).

- 240 draws: 212 ok (88.3%), 28 rejected by the shared_Neff_guard (extreme
  (r, DN_re, kappa10) corners where the shared background Delta_Neff grows past
  the documented guard; rejection is explicit and traceable, never silent).
- fast runtime: median 5.34 s, p95 8.82 s, max 14.2 s per production point.
- adaptive frequency grid: median 236 nodes, p95 266, max 286 (sparse in smooth
  regions, refined at knee/stiff/cutoff).
- WKB/adiabatic handoff defect (per-mode max eps): median 6.84e-4, max 6.84e-4
  (flat: the dominant handoff defect is set by the z~8 envelope, not by parameter
  extremes).
- local a-posteriori combined relative Delta-Neff error estimate (11-category
  budget, `estimate_local_error`): median 4.0e-4 for physically significant
  points; the estimate saturates at ~1.0 (100% relative) only in the DN -> 0
  corner (DN < ~1e-5, absolute error <= 1e-5, physically unobservable).

## 3. Why matched z8 (not z7, not LSODA)

- z7 matched reference is NOT converged per mode: reference itself moves
  ~3.6e-4 dex between z_tail 7 and 8, so z7 overstates the engine error
  (dex ~8e-4, rel ~1.8e-3). z8 matched is the physics verdict.
- LSODA (default anchor) sits at -1.68e-3 relative DN vs the continuous-sigma
  reference on the default point; fast production is -2.9e-4..+3.0e-4 (6-30x
  closer) at matched settings. LSODA numbers are reported for context only.

## 4. Honest limits (not threshold-tuned)

- The integrated-Delta-Neff 1e-4 acceptance is NOT met: measured residual is
  median 4.3e-4 (z8, f_res=1) and ~2.9e-4 (default point, deep oracle f_res=2,
  rtol=1e-10). The residual is engine-architecture level: per-mode dex residuals
  ~2.6-3.1e-4 are internally converged against h/z_tail/freq_res/phase variations
  (<1e-4 dex), i.e. the remaining term is the frozen-z Magnus + transition-refine
  envelope, not a tuning artifact.
- Per-mode signal-region Omega_GW relative error < 1e-3 IS met at all 9 matched
  points (max 7.09e-4). Per-mode error < 1e-4 is NOT met (would require a
  higher-order per-mode integrator or continuous-sigma fast sigma).
- 100x-vs-LSODA production runtime: the honest comparison at matched accuracy
  settings (z8 production) is 3-16x vs the documented 18.6 s LSODA anchor; the
  100x+ figure is only true for the plain-grid coarser mode (12 ms, -0.16% DN),
  not for the accuracy mode that meets the 1e-3 physics gate.


## 5. Layer C posterior validation (fast production vs continuous-sigma reference)

Run 2026-09-03 on this host, fixed seed 20260903. Driver: `scripts/importance_posterior.py`
(IS phase) and `scripts/cobaya_posterior_fast_vs_reference.py` (bounded Cobaya MCMC scaffold).
Artifacts in `docs/mcmc_posterior/`: `mock_truth.json`, `is_draws.npz`, `is_posterior.json`,
`is_pointwise.json`, `is_report.json`.

Mock likelihood: 11 frequency bins (4 PTA + 4 knee + 3 LVK, log10 f/Hz in [-8.6, 2.6]),
truth = reference spectrum at fiducial (log10r, n_t) = (-2.0, 0.0); Gaussian sigma_dex = 0.05;
priors log10r in [-5, 0], n_t in [-0.5, 0.5]. 9000 Gaussian-proposal IS draws (mu_r=-2,
sig_r=0.03, mu_nt=0, sig_nt=0.25).

Pointwise engine certification (the likelihood-aware gate): 240 posterior-bulk points
(drawn from the IS posterior, seed+7), each solved by BOTH fast production z8
(grid_independent freq_res=1, transition_refine, phase_max sub-stepping) and the
continuous-sigma reference (DOP853 rtol=1e-9, z_tail=8) at the SAME 11 bins. 240/240 ok.

Key fix in this run: fast now receives `eval_freqs=BINS` so the 11 likelihood bins are
NATIVE solve-grid nodes (`SGWB_iter_fast(..., eval_freqs=...)`, added 2026-09-03). The
driver previously CubicSpline-interpolated the fast spectrum onto the bins across
0.105-dex nodes of the steep high-frequency wall, producing a 1e-2-dex-level interpolation
artifact at r ~ 1e-2 (dex max 6.86e-3, |dll| max 0.1035 -> gate FAIL). With native-node
bins the residual drops to the solve-node level:

| statistic (240 points, fast vs reference) | value |
|---|---:|
| per-bin dex max (11 bins) | 3.10e-4 (at logf=-2.0) |
| dex max over all points/bins | 3.10e-4 |
| dex p95 | 3.02e-4 |
| dex median | 8.0e-5 |
| Delta logL (ref - fast) max_abs | 7.30e-3 |
| Delta logL p95_abs | 4.70e-3 |
| Delta logL mean | -4.62e-4 |

Posterior statistics (IS, ess = 4167.4 >= 2000):
- log10r: mean -1.99974, std 0.01586, median -1.99966, p16/p84 -2.01550/-1.98388
- n_t: mean 0.00798, std 0.28878 (n_t is NOT constrained by the cr=1 mock: with
  `cr=1` the spectrum is n_t-independent, so n_t posteriors are prior-dominated and
  reported as documentation only)
- posterior shift after e^{Delta logL} reweighting (reference-consistent posterior):
  log10r -0.0011 sigma, n_t +0.0002 sigma (both < 0.1 sigma)

Verdicts (`is_report.json`): ess_is>=2000 PASS, dll_max_abs<0.1 PASS, shift<0.1sigma PASS.

Honest notes: (1) the ESS>=2000 and shift claims rest on the IS posterior from 9000 fast
production draws (ess 4167.4) with the pointwise fast-vs-reference dex gate at the
posterior bulk; full two-chain KS/Wasserstein/KL would need a reference-engine MCMC
chain at ~350-935 s/pt on this host and is not attempted here - the IS e^{Delta logL}
reweighting is the exact estimator used instead (valid because |Delta logL| is tiny,
max 7.3e-3). (2) The bounded real-Cobaya chains in `docs/mcmc_posterior/chains/` are
~30-row scaffold runs (documented as not converged) - the production-grade posterior
claim is the IS chain above. (3) n_t shift is meaningless for cr=1 (flat direction);
the physically meaningful parameter shift is log10r = -0.0011 sigma.
