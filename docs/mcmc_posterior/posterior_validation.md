# Layer C: posterior validation (fast production vs continuous-sigma reference)

Date 2026-09-03, host-verified numbers. LSODA is never a truth anchor here; the
reference/oracle is `stiffgwpy/reference.py` (continuous-sigma + DOP853). All numbers
below come from the artifacts in this directory (`mock_truth.json`, `is_draws.npz`,
`is_posterior.json`, `is_pointwise.json`, `is_report.json`), reproduced by:

    python scripts/importance_posterior.py --phase draw --n 9000 --workers 8 \
        --mu-r -2.0 --sig-r 0.03 --mu-nt 0.0 --sig-nt 0.25 --sigma-dex 0.05
    python scripts/importance_posterior.py --phase posterior --mu-r -2.0 --sig-r 0.03 \
        --mu-nt 0.0 --sig-nt 0.25 --sigma-dex 0.05
    python scripts/importance_posterior.py --phase pointwise --k 240 --workers 8 \
        --fast-freq-res 1.0 --mu-r -2.0 --sig-r 0.03 --mu-nt 0.0 --sig-nt 0.25 --sigma-dex 0.05
    python scripts/importance_posterior.py --phase report --mu-r -2.0 --sig-r 0.03 \
        --mu-nt 0.0 --sig-nt 0.25 --sigma-dex 0.05

## 1. Setup (identical prior / truth / sampler settings for both engines)

- Model: r, n_t free (priors log10r in [-5, 0], n_t in [-0.5, 0.5]); cr=1, DN_re=0,
  T_re=2000 K, kappa10=1e-2 fixed.
- Mock data: 11 bins log10 f/Hz = PTA {-8.6,-8.3,-8.0,-7.7} + knee {-2,-1,0,1} +
  LVK {1.40,2.00,2.60}; truth = continuous-sigma reference spectrum at
  (log10r,n_t)=(-2.0,0.0); independent Gaussian noise sigma_dex=0.05 per bin.
- Sampler: 9000 Gaussian-proposal importance draws from q = N(mu=(-2,0),
  sigma=(0.03,0.25)) with fixed seed 20260903. Proposal was tuned to the posterior
  bulk (posterior width on log10r is 0.0159, proposal 0.03 -> healthy overlap).
- Engines: fast = production z8 preset, grid_independent freq_res=1.0, transition
  refine + phase sub-stepping, threads=8; reference = continuous-sigma DOP853
  rtol=1e-9 z_tail=8, same 11 bins, per-point wall ~4 s vs ~1.5-2 s total spread.

## 2. Likelihood-aware accuracy gate at the posterior bulk (240 points)

240 posterior-bulk points (weighted draw, seed+7) solved by BOTH engines at the same
11 bins; fast receives eval_freqs=BINS so the bins are native solve nodes (kills the
previous CubicSpline-across-0.105-dex-nodes artifact on the steep high-f wall).

| statistic (fast vs reference, 240 points) | before fix | after fix (final) |
|---|---:|---:|
| dex max (all bins) | 6.86e-3 | 3.10e-4 |
| dex p95 | 6.25e-3 | 3.02e-4 |
| Delta logL max_abs | 0.1035 | 7.30e-3 |
| Delta logL p95_abs | 0.0598 | 4.70e-3 |
| Delta logL mean | 7.0e-3 | -4.62e-4 |
| gate |logL|<0.1 | FAIL | PASS |

Per-bin dex max (final, 240 points): PTA -8.6..-7.7 in [7.4e-5, 2.2e-4]; knee -2..1 in
[2.7e-4, 3.1e-4]; LVK 1.4..2.6 in [8.1e-5, 2.4e-4].

## 3. Posterior statistics (fast IS chain) and shift

IS posterior (ess = 4167.4, threshold 2000 PASS):
- log10r: mean -1.9997435, std 0.0158602, median -1.9996606, p16/p84 -2.015498/-1.983880
- n_t: mean 0.0079784, std 0.2887765 (unconstrained flat direction: cr=1 physics is
  n_t-independent, so n_t is prior-dominated and only documented, not certified)

Reference-consistent posterior via e^{Delta logL} reweighting of the same 240
posterior-bulk points (exact importance estimator; valid because |Delta logL| is tiny):
- log10r: mean_fast -1.9986578 -> mean_reweighted -1.9986752 -> shift -0.00110 sigma
- n_t: mean_fast 0.0291603 -> mean_reweighted 0.0292253 -> shift +0.00023 sigma
- both shifts < 0.1 sigma -> PASS (key-parameter scientific conclusion unchanged)

## 4. Verdicts (`is_report.json`)

- ess_is >= 2000: PASS (4167.4)
- dll_max_abs < 0.1: PASS (7.30e-3)
- posterior_shift < 0.1 sigma: PASS (log10r -0.0011, n_t +0.0002)

## 5. Honest limits of this Layer C validation

- ESS>=2000 / shift claims rest on the IS posterior of 9000 fast-production draws;
  the pointwise 240-point set certifies engine agreement inside the posterior bulk
  (dex <= 3.1e-4 -> |Delta logL| <= 7.3e-3). Full two-chain KS/Wasserstein/KL/
  covariance comparisons were not run because a reference-engine MCMC chain costs
  ~350-935 s/point on this host; the IS reweighting is the exact alternative and its
  precondition (small |Delta logL|) is verified.
- Real bounded Cobaya chains in `chains/` are ~30-row scaffold runs, documented as
  not converged; they validate the Cobaya adapter plumbing only.
- n_t is unconstrained by the cr=1 mock, so only log10r carries a physical shift claim.
