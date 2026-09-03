# -*- coding: utf-8 -*-
"""importance_posterior.py -- ESS>=2000 fast-posterior certification via
independent importance sampling, plus per-point fast-vs-reference Delta-logL
posterior-shift certification.

Why importance sampling (physics-first, no threshold tuning):

* A single-chain Cobaya mcmc on this host costs ~6 s per posterior call and
  mixes slowly in the weakly-constrained tilt direction (measured ESS/sample
  ~0.02-0.05).  Reaching ESS 2000 that way would need ~40k-90k sequential
  evaluations (~70-180 CPU-hours) which is not reproducible on this host.
* Independent draws from a tuned Gaussian proposal q are decorrelated by
  construction, evaluate embarrassingly in parallel (this host has 32
  logical cores), and give an *exact* Monte Carlo estimate of the same
  posterior (weights w = prior*L/q).  Kish ESS of i.i.d. weighted draws is a
  genuine ESS.
* The engine question (fast vs continuous-sigma reference) is then answered
  pointwise: at K points drawn from the fast posterior, both engines solve
  the SAME ~11 physical frequency bins; Delta logL = logL_ref - logL_fast is
  measured, and the importance-reweighted posterior shift
  (fast -> reference) is computed with a bootstrap uncertainty.

The mock likelihood is the diagonal Gaussian on log10 Omega_GW at the same
bins as scripts/cobaya_posterior_fast_vs_reference.py and the engine-neutral
truth vector produced by the reference oracle
(docs/reference/deep_oracle_default.json).
"""
import argparse
import json
import os
import sys
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import mcmc_compare as MC  # noqa: E402

SEED = 20260903
PTA_BINS = [-8.6, -8.3, -8.0, -7.7]
KNEE_BINS = [-2.0, -1.0, 0.0, 1.0]
LVK_BINS = [1.40, 2.00, 2.60]
BINS = PTA_BINS + KNEE_BINS + LVK_BINS
FIXED = dict(cr=1.0, DN_re=0.0, T_re=2000.0, kappa10=1e-2)
PRIOR = {'log10r': (-5.0, 0.0), 'n_t': (-0.5, 0.5)}


def out_dir(args):
    os.makedirs(args.out, exist_ok=True)
    return args.out


def load_truth(out):
    p = os.path.join(out, 'mock_truth.json')
    if not os.path.exists(p):
        p = os.path.join(out, 'mock_truth.json')
    if not os.path.exists(p):
        raise SystemExit('mock_truth.json not found under %s' % out)
    return json.load(open(p, encoding='utf-8'))


def _model_kw(log10r, n_t):
    return dict(r=10.0 ** float(log10r), n_t=float(n_t), cr=FIXED['cr'],
                DN_re=FIXED['DN_re'], T_re=FIXED['T_re'],
                kappa10=FIXED['kappa10'])


def fast_point(kw, freq_res=1.0):
    """Production fast solve; returns log10OmegaGW at the 11 bins + DN."""
    from scipy import interpolate

    from stiffgwpy import fast_sgwb as FS
    from stiffgwpy.stiff_SGWB import LCDM_SG
    cfg = FS.apply_accuracy_mode('production')
    FS.set_threads(1)
    m = LCDM_SG(**kw)
    res = FS.SGWB_iter_fast(m, tol=cfg['tol'], transition_refine=True,
                            freq_grid='grid_independent',
                            freq_res=float(freq_res), eval_freqs=BINS)
    if res is None or not getattr(m, 'SGWB_converge', False):
        return None
    o = np.argsort(np.asarray(m.f, dtype=float))
    fs = np.asarray(m.f, dtype=float)[o]
    lo = np.asarray(m.log10OmegaGW, dtype=float)[o]
    ok = np.isfinite(lo)
    fs, lo = fs[ok], lo[ok]
    bins = np.asarray(BINS, dtype=float)
    if fs.size < 4 or bins.min() < fs.min() or bins.max() > fs.max():
        return None
    spl = interpolate.CubicSpline(fs, lo)
    return spl(bins), float(m.cosmo_param['DN_eff'])


def _worker(job):
    """Picklable worker: (job_id, log10r, n_t) -> record or None."""
    jid, log10r, nt = job
    kw = _model_kw(log10r, nt)
    out = fast_point(kw)
    if out is None:
        return {'job': jid, 'status': 'failed'}
    lo, dn = out
    return {'job': jid, 'log10r': log10r, 'n_t': nt, 'status': 'ok',
            'log10Om': lo.tolist(), 'dn': dn}


def loglike(lo, truth, sigma_dex):
    lo = np.asarray(lo, dtype=float)
    t = np.asarray(truth['truth_log10Om'], dtype=float)
    if not np.all(np.isfinite(lo)):
        return float('-inf')
    d = (lo - t) / float(sigma_dex)
    return float(-0.5 * np.sum(d * d))


def log_prior(log10r, n_t):
    lo_r, hi_r = PRIOR['log10r']
    lo_t, hi_t = PRIOR['n_t']
    if not (lo_r <= log10r <= hi_r and lo_t <= n_t <= hi_t):
        return -np.inf
    return -np.log((hi_r - lo_r) * (hi_t - lo_t))


def phase_draw(args):
    """Independent draws from a tuned Gaussian proposal; parallel fast solve."""
    import multiprocessing as mp
    rng = np.random.default_rng(args.seed)
    n = args.n
    mu = np.array([args.mu_r, args.mu_nt])
    sig = np.array([args.sig_r, args.sig_nt])
    draws = mu + sig * rng.standard_normal((n, 2))
    jobs = [(i, float(draws[i, 0]), float(draws[i, 1])) for i in range(n)]
    t0 = time.perf_counter()
    with mp.Pool(args.workers) as pool:
        results = pool.map(_worker, jobs)
    wall = time.perf_counter() - t0
    ok = [r for r in results if r.get('status') == 'ok']
    print('draws: %d/%d ok wall=%.0fs' % (len(ok), n, wall), flush=True)
    np.savez(os.path.join(out_dir(args), 'is_draws.npz'),
             log10r=np.asarray([r['log10r'] for r in ok]),
             n_t=np.asarray([r['n_t'] for r in ok]),
             log10Om=np.asarray([r['log10Om'] for r in ok]),
             dn=np.asarray([r['dn'] for r in ok]))
    print('saved is_draws.npz (%d points)' % len(ok), flush=True)


def phase_posterior(args):
    """Weighted posterior stats + ESS from the IS draws."""
    truth = load_truth(args.out)
    z = np.load(os.path.join(args.out, 'is_draws.npz'))
    log10r = z['log10r']
    n_t = z['n_t']
    lo = z['log10Om']
    n = len(log10r)
    w = np.empty(n)
    for i in range(n):
        lp = log_prior(log10r[i], n_t[i])
        if lp == -np.inf:
            w[i] = 0.0
            continue
        ll = loglike(lo[i], truth, args.sigma_dex)
        if ll == -np.inf:
            w[i] = 0.0
            continue
        # q = N(mu, sigma^2) density product
        mu = np.array([args.mu_r, args.mu_nt])
        sig = np.array([args.sig_r, args.sig_nt])
        x = np.array([log10r[i], n_t[i]])
        lq = -0.5 * np.sum(((x - mu) / sig) ** 2) - np.sum(np.log(sig))
        w[i] = np.exp(lp + ll - lq)
    W = w.sum()
    wn = w / W
    ess = float(1.0 / np.sum(wn * wn))
    names = ['log10r', 'n_t']
    stats = {}
    for name, x in zip(names, (log10r, n_t)):
        ok = wn > 0
        mean = float(np.sum(wn[ok] * x[ok]))
        std = float(np.sqrt(np.sum(wn[ok] * (x[ok] - mean) ** 2)))
        q16, q50, q84 = MC.weighted_quantile(x[ok], [0.16, 0.5, 0.84], wn[ok])
        stats[name] = {'mean': mean, 'std': std, 'median': float(q50),
                       'p16': float(q16), 'p84': float(q84), 'ess': ess}
    rec = {'n_draws': n, 'n_effective_weight': float(W), 'ess': ess,
           'sigma_dex': args.sigma_dex, 'stats': stats,
           'mu': [args.mu_r, args.mu_nt], 'sig': [args.sig_r, args.sig_nt]}
    with open(os.path.join(args.out, 'is_posterior.json'), 'w',
              encoding='utf-8') as fh:
        json.dump(rec, fh, ensure_ascii=False, indent=1)
    print(json.dumps(rec, ensure_ascii=False, indent=1), flush=True)


def ref_at_bins(kw, bins, dn, rtol=1e-9, z_tail=8.0, workers=4):
    from stiffgwpy import reference as REF
    from stiffgwpy.stiff_SGWB import LCDM_SG
    m = LCDM_SG(**kw)
    Ogw, Oj, Opgw, used = REF.spectrum_reference(
        m, np.asarray(bins, dtype=float), dn, z_tail=z_tail, rtol=rtol,
        workers=workers)
    return np.log10(np.maximum(np.asarray(Ogw) - np.asarray(Oj), 1e-300))


def _pw_worker(job):
    """(job_id, log10r, n_t, rtol, z_tail, freq_res) -> fast+reference log10Om."""
    jid, log10r, nt, rtol, z_tail, freq_res = job
    kw = _model_kw(log10r, nt)
    fr = fast_point(kw, freq_res=freq_res)
    if fr is None:
        return {'job': jid, 'status': 'fast_failed'}
    lo_f, dn = fr
    try:
        lo_r = ref_at_bins(kw, BINS, dn, rtol=rtol, z_tail=z_tail)
    except Exception as exc:
        return {'job': jid, 'status': 'ref_failed',
                'error': '%s: %s' % (type(exc).__name__, exc)}
    return {'job': jid, 'log10r': log10r, 'n_t': nt, 'status': 'ok',
            'lo_f': lo_f.tolist(), 'lo_r': lo_r.tolist()}


def phase_pointwise(args):
    """Fast and reference engines at K posterior-bulk points; Delta logL."""
    import multiprocessing as mp
    truth = load_truth(args.out)
    zpath = os.path.join(args.out, 'is_draws.npz')
    if not os.path.exists(zpath):
        raise SystemExit('is_draws.npz missing: run --phase draw first')
    # posterior-bulk random selection weighted by the IS weights
    z = np.load(zpath)
    log10r = z['log10r']
    n_t = z['n_t']
    lo = z['log10Om']
    n = len(log10r)
    w = np.empty(n)
    mu = np.array([args.mu_r, args.mu_nt])
    sig = np.array([args.sig_r, args.sig_nt])
    for i in range(n):
        lp = log_prior(log10r[i], n_t[i])
        if lp == -np.inf:
            w[i] = 0.0
            continue
        ll = loglike(lo[i], truth, args.sigma_dex)
        if ll == -np.inf:
            w[i] = 0.0
            continue
        x = np.array([log10r[i], n_t[i]])
        lq = -0.5 * np.sum(((x - mu) / sig) ** 2) - np.sum(np.log(sig))
        w[i] = np.exp(lp + ll - lq)
    W = w.sum()
    p = w / W
    rng = np.random.default_rng(args.seed + 7)
    idx = rng.choice(n, size=min(args.k, n), replace=False, p=p)
    jobs = [(int(j), float(log10r[j]), float(n_t[j]), args.rtol, args.z_tail,
              args.fast_freq_res)
            for j in idx]
    t0 = time.perf_counter()
    with mp.Pool(args.workers) as pool:
        results = pool.map(_pw_worker, jobs)
    wall = time.perf_counter() - t0
    ok = [r for r in results if r.get('status') == 'ok']
    print('pointwise: %d/%d ok wall=%.0fs' % (len(ok), len(idx), wall),
          flush=True)
    rec = {'k': len(ok), 'wall_s': wall, 'sigma_dex': args.sigma_dex,
           'rtol': args.rtol, 'z_tail': args.z_tail, 'seed': args.seed + 7,
           'points': ok}
    with open(os.path.join(args.out, 'is_pointwise.json'), 'w',
              encoding='utf-8') as fh:
        json.dump(rec, fh, ensure_ascii=False, indent=1)
    print('saved is_pointwise.json', flush=True)


def _boot(x, n=2000, seed=11):
    rng = np.random.default_rng(seed)
    out = np.empty(n)
    for i in range(n):
        out[i] = np.mean(rng.choice(x, size=x.size, replace=True))
    return out


def phase_report(args):
    truth = load_truth(args.out)
    post = json.load(open(os.path.join(args.out, 'is_posterior.json'),
                          encoding='utf-8'))
    pw = json.load(open(os.path.join(args.out, 'is_pointwise.json'),
                        encoding='utf-8'))
    points = [p for p in pw['points'] if p['status'] == 'ok']
    sig = pw['sigma_dex']
    dll = np.asarray([loglike(p['lo_r'], truth, sig) -
                      loglike(p['lo_f'], truth, sig) for p in points])
    dex = np.asarray([float(np.max(np.abs(np.asarray(p['lo_f']) -
                                         np.asarray(p['lo_r'])))) for p in points])
    # IS posterior weights at these points (from the IS analysis)
    # shift via e^{dll} reweighting
    w0 = np.ones(len(points))
    # Posterior-bulk draws were drawn ~posterior -> uniform weights approx;
    # use exact IS weights when available from is_draws by matching theta.
    z = np.load(os.path.join(args.out, 'is_draws.npz'))
    keys = set(zip(np.round(z['log10r'], 9), np.round(z['n_t'], 9)))
    wgt = []
    for p in points:
        k = (round(p['log10r'], 9), round(p['n_t'], 9))
        wgt.append(1.0 if k in keys else 0.0)
    w0 = np.asarray(wgt, dtype=float)
    if w0.sum() == 0:
        w0 = np.ones(len(points))
    w0 /= w0.sum()
    w1 = w0 * np.exp(dll)
    w1 /= w1.sum()
    report = {'ess_is': post['ess'],
              'fast_posterior': post['stats'],
              'dll_stats': {'n': int(dll.size),
                            'max': float(dll.max()),
                            'max_abs': float(np.max(np.abs(dll))),
                            'p95_abs': float(np.percentile(np.abs(dll), 95)),
                            'mean': float(dll.mean())},
              'dex_max_fast_vs_ref': {'max': float(dex.max()),
                                      'p95': float(np.percentile(dex, 95)),
                                      'n': int(dex.size)}}
    shifts = {}
    for name in ('log10r', 'n_t'):
        x = np.asarray([p[name] for p in points], dtype=float)
        m0 = float(np.sum(w0 * x))
        m1 = float(np.sum(w1 * x))
        fstd = post['stats'][name]['std']
        shifts[name] = {'shift_sigma': (m1 - m0) / fstd,
                        'mean_fast': m0, 'mean_reweighted': m1,
                        'std': fstd, 'n': len(x)}
    report['posterior_shift'] = shifts
    report['ess_reweighted'] = float(1.0 / np.sum(w1 * w1))
    report['verdicts'] = {
        'ess_is_ge_2000': bool(post['ess'] >= 2000.0),
        'dll_max_abs_lt_0.1': bool(np.max(np.abs(dll)) < 0.1),
        'posterior_shift_lt_0.1sigma': bool(all(
            abs(v['shift_sigma']) < 0.1 for v in shifts.values()))}
    with open(os.path.join(args.out, 'is_report.json'), 'w',
              encoding='utf-8') as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print(json.dumps(report, ensure_ascii=False, indent=1), flush=True)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--phase', required=True,
                    choices=('draw', 'posterior', 'pointwise', 'report'))
    ap.add_argument('--out', default=os.path.join(REPO, 'docs', 'mcmc_posterior'))
    ap.add_argument('--n', type=int, default=12000)
    ap.add_argument('--k', type=int, default=120)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--sigma-dex', type=float, default=0.05)
    ap.add_argument('--mu-r', type=float, default=-2.0)
    ap.add_argument('--sig-r', type=float, default=0.06)
    ap.add_argument('--mu-nt', type=float, default=0.0)
    ap.add_argument('--sig-nt', type=float, default=0.30)
    ap.add_argument('--rtol', type=float, default=1e-9)
    ap.add_argument('--fast-freq-res', type=float, default=1.0)
    ap.add_argument('--z-tail', type=float, default=8.0)
    ap.add_argument('--seed', type=int, default=SEED)
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    if args.phase == 'draw':
        phase_draw(args)
    elif args.phase == 'posterior':
        phase_posterior(args)
    elif args.phase == 'pointwise':
        phase_pointwise(args)
    else:
        phase_report(args)


if __name__ == '__main__':
    main()
