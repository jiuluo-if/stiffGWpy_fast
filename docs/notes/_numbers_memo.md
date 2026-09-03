# session numbers memo (authoritative; folded into audit docs)

## engines / settings
- fast production (h=0.01, col_step=4, z_tail=8, freq_res=1, freq_grid=adaptive or
  grid_independent, transition_refine=True, phase_max=0.5, threads=8):
  measured 1.2-14.2 s/self-consistent pt this host (median 5.34 s over 212-point sweep;
  default point 4.1 s).
- reference: continuous-sigma DOP853 per mode rtol=1e-9; z8 242-freq ~350-375 s/point single;
  deep oracle default (z8 484 pt, rtol 1e-10) ~935 s.
- LSODA anchor (documented, z5): 18.56 s/point. NOT used as truth anywhere.

## deep oracle default (matched z8, 484 grid, rtol 1e-10) docs/reference/deep_oracle_default.json
- fast DN_eff=0.002262265; ref DN_gw=0.0022629309 -> DN rel = -2.94e-4 (abs -6.66e-7)
- signal rel p95 6.35e-4 max 7.12e-4 ; dex max 3.09e-4
- transition rel p95 6.74e-4 max 7.12e-4
- all-band rel max 1.2e-3 (low-f tail dominated)

## matched z8 sweep COMPLETE (rtol 1e-9, freq_res 1, 242-grid; 9 singles) -> docs/paramsweep_z8/
- default +3.01e-4 (sig rel max 5.98e-4), stiff -4.32e-4 (6.75e-4), lowT +1.46e-3* (6.02e-4),
  highT +1.27e-4 (6.92e-4), rad_dom +7.18e-4 (6.43e-4), tiny_r +3.8e-5 (7.09e-4),
  transition -7.56e-4 (6.64e-4), cr0_blue +6.10e-4 (6.58e-4), extreme -4.34e-4 (6.92e-4)
  * lowT: DN-of-DN artifact (DN~5.2e-8); signal gate still passes.
- DN rel abs: median 4.34e-4, p95 1.18e-3, max 1.46e-3 (lowT) -> NOT <1e-4 (honest limit)
- signal: rel max 7.09e-4, p95 7.02e-4 -> PASS <1e-3; transition same -> PASS
- validation_summary.json + validation_summary.md written.

## matched z7 sweep (superseded; ref not converged at z7)
- dex ~7.9-8.1e-4 -> rel_sig_max ~1.8e-3 overstates engine error; ref z7-vs-z8 moves ~3.6e-4 dex.
  z8 matched is the physics verdict.

## 240-point Sobol sweep (production z8) docs/paramsweep_ref/fast_sweep.jsonl
- 240 pts: 212 ok (88.3%), 28 shared_Neff_guard (extreme corners, explicit not silent)
- runtime median 5.34 s, p95 8.82 s, max 14.2 s
- adaptive freq grid median 236 nodes (p95 266, max 286)
- handoff eps per-mode max median 6.84e-4 (flat across parameters)
- local a-posteriori combined rel DN error estimate: median 4.0e-4 (11-category budget);
  saturates ~1.0 only in DN->0 corner (abs <=1e-5, unobservable)

## runtime pareto (documented, default point)
- reference z5 171.7 s / z8 deep 935 s / z8 rtol1e-9 359 s
- LSODA z5 18.56 s
- fast_transition z5 0.28 s (66x vs LSODA z5; 613x vs ref z5)
- fast plain grid 0.012-0.24 s (77x-1541x vs LSODA; DN ~-0.16%)
- fast production z8 4.1 s -> ~87x vs z8 ref 359 s; 4.5x vs LSODA at 6-30x smaller DN error
- HONEST: 100x-vs-LSODA only for plain-grid mode; production z8 is 4-16x vs LSODA anchor.

## local error budget (estimate_local_error, 11 categories)
background_model / sigma_transition / ode_integration / horizon_crossing / wkb_handoff /
interpolation / frequency_grid / quadrature / tail_approximation / floating_point /
self_consistency; combined = max(model,transition) + RSS(rest). Anchors: phase_max
Richardson 8e-6 (pm^2); transition_refine 1e-4 (plain 2.4e-3, sigma_exact 4.2e-3);
grid_uniform f_res1 5.6e-4 / f_res2 2e-5 / construct 1.6e-3; tail z8 2e-5 (z10 2.4e-7);
interp 1e-9; self-cons floor 1e-5.
