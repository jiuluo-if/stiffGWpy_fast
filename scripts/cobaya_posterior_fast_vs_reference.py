# -*- coding: utf-8 -*-
"""cobaya_posterior_fast_vs_reference.py -- likelihood-aware posterior
certification of the fast engine against the independent continuous-sigma
DOP853 reference (never against LSODA).

Method (physics-first, and built around the measured per-point engine error,
which is ~3e-4 dex in the signal band / ~3e-4 relative in Delta N_eff):

1. ``truth``   : store the mock data vector.  Central values are the
   *reference* (continuous-sigma DOP853 oracle) spectrum of the fiducial
   point at a set of physical frequency bins (PTA nHz band, reheating-knee
   band and LVK band).  The mock is therefore engine-neutral: it is built
   from the reference engine, not from fast.
2. ``chain``   : Cobaya mcmc, engine=fast (production accuracy mode), fixed
   seed, fixed initial point, sampling (log10r, n_t).  The mock likelihood
   is a diagonal Gaussian on the log10-Omega_GW bins.
3. ``pointwise`` : K parameter points are drawn uniformly at random from the
   fast posterior (thinned chain).  At every point BOTH engines are solved
   at the SAME frequency bins; the per-point log-likelihood difference
   Delta logL = logL(reference) - logL(fast) is measured directly.
4. ``report`` : posterior stats (mean/std/MAP/16-84%/ESS/covariance) for the
   fast chain, the measured Delta-logL distribution over the posterior bulk,
   and the importance-reweighted posterior shift (fast -> reference) in
   sigma units with a bootstrap uncertainty.  Certification verdicts.

The reference leg is cheap here because only ~11 frequency bins per point
are solved (a full 242-frequency DOP853 solve costs ~3.5 min; the 11-bin
solve costs seconds).  The integrated Delta N_eff contribution is reported
separately using the full-grid measurements of scripts/validate_fast_vs_
reference.py (bounded |Delta N_eff| abs error), so the two error channels
stay explicit instead of being mixed into one number.

Usage (from the repo root):

    python scripts/cobaya_posterior_fast_vs_reference.py --phase truth
    python scripts/cobaya_posterior_fast_vs_reference.py --phase chain --samples 4000
    python scripts/cobaya_posterior_fast_vs_reference.py --phase pointwise --k 60
    python scripts/cobaya_posterior_fast_vs_reference.py --phase report
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cobaya.likelihood import Likelihood  # noqa: E402
from scipy import interpolate  # noqa: E402

# re-use the (LSODA-era but engine-agnostic) chain statistics helpers
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mcmc_compare as MC  # noqa: E402

SEED = 20260902
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Frequency bins (log10 f/Hz) that pin the physics of the default fiducial:
#   PTA nHz band  -> tensor tilt n_t (f^{n_t} plateau),
#   knee band     -> stiff/reheating transition region (fast-vs-ref worst),
#   LVK band      -> amplitude in the f^{-2}-like tail (amplitude/r).
PTA_BINS = [-8.6, -8.3, -8.0, -7.7]
KNEE_BINS = [-2.0, -1.0, 0.0, 1.0]
LVK_BINS = [1.40, 2.00, 2.60]
BINS = PTA_BINS + KNEE_BINS + LVK_BINS

FIDUCIAL = dict(log10r=-2.0, n_t=0.0)

# Fixed non-sampled physics of the mock experiment (default benchmark point).
FIXED = dict(cr=1.0, DN_re=0.0, T_re=2000.0, kappa10=1e-2)


def out_dir(args):
    os.makedirs(args.out, exist_ok=True)
    return args.out


def oracle_truth(bins):
    """Reference (DOP853 oracle) log10 Omega_GW at the bins, engine-neutral."""
    p = os.path.join(REPO, 'docs', 'archive', 'reference', 'deep_oracle_default.json')
    d = json.load(open(p, encoding='utf-8'))
    order = np.argsort(np.asarray(d['G'], dtype=float))
    G = np.asarray(d['G'], dtype=float)[order]
    lo = np.asarray(d['ref']['log10Om'], dtype=float)[order]
    ok = np.isfinite(G) & np.isfinite(lo)
    spl = interpolate.CubicSpline(G[ok], lo[ok])
    return spl(np.asarray(bins, dtype=float))


class SGWB_mock(Likelihood):
    """Diagonal-Gaussian mock likelihood on log10 Omega_GW at fixed bins.

    Data vector (truth_log10Om) is engine-neutral: it is produced by the
    continuous-sigma DOP853 reference at the fiducial point, so a fast chain
    is scored exactly like a reference chain and any per-point Delta logL is
    a genuine engine-accuracy effect.
    """
    bins: list
    truth_log10Om: list
    sigma_dex: float = 0.05

    def initialize(self):
        # Direct instantiation (unit tests) has no options injected yet.
        if getattr(self, 'bins', None) is None:
            return
        self.log.debug('mock bins=%d sigma_dex=%.3f',
                       len(self.bins), self.sigma_dex)

    def get_requirements(self):
        return {'f': None, 'omGW_stiff': None}

    def logp(self, _derived=None, **params_values):
        f_theory = self.provider.get_result('f')
        lo_theory = self.provider.get_result('omGW_stiff')
        order = np.argsort(np.asarray(f_theory, dtype=float))
        fs = np.asarray(f_theory, dtype=float)[order]
        lo = np.asarray(lo_theory, dtype=float)[order]
        ok = np.isfinite(lo)
        fs = fs[ok]
        lo = lo[ok]
        bins = np.asarray(self.bins, dtype=float)
        if fs.size < 4 or bins.min() < fs.min() or bins.max() > fs.max():
            return float('-inf')
        spl = interpolate.CubicSpline(fs, lo)
        lo_at_bins = spl(bins)
        return mock_loglike(lo_at_bins,
                            {'truth_log10Om': self.truth_log10Om},
                            self.sigma_dex)


def phase_truth(args):
    truth = oracle_truth(BINS)
    rec = {'fiducial': FIDUCIAL, 'fixed': FIXED, 'bins': BINS,
           'truth_log10Om': truth.tolist(),
           'source': 'deep_oracle_default (continuous-sigma DOP853, z8, rtol=1e-10)',
           'sigma_dex_default': 0.05}
    with open(os.path.join(out_dir(args), 'mock_truth.json'), 'w',
              encoding='utf-8') as fh:
        json.dump(rec, fh, ensure_ascii=False, indent=1)
    print(json.dumps(rec, ensure_ascii=False, indent=1), flush=True)


def _load_truth(out):
    p = os.path.join(out, 'mock_truth.json')
    if not os.path.exists(p):
        raise SystemExit('mock_truth.json missing: run --phase truth first')
    return json.load(open(p, encoding='utf-8'))


def model_kwargs(log10r, n_t, truth=None):
    r = 10.0 ** float(log10r)
    return dict(r=r, n_t=float(n_t), cr=FIXED['cr'], DN_re=FIXED['DN_re'],
                T_re=FIXED['T_re'], kappa10=FIXED['kappa10'])


def build_info(truth, samples, sigma_dex=0.05, seed=SEED):
    """Cobaya run-info: stiffGW(fast) + diagonal-Gaussian mock likelihood."""
    bins = np.asarray(truth['bins'], dtype=float)
    y = np.asarray(truth['truth_log10Om'], dtype=float)
    info = {
        'theory': {'stiffgwpy.cobaya.stiffGW.stiffGW': {
            'engine': 'fast', 'fallback': True, 'accuracy_mode': 'production',
            'fast_threads': 8, 'likelihood_sigma': 5e-4, 'dlogl_tol': 1e-3}},
        'likelihood': {'cobaya_posterior_fast_vs_reference.SGWB_mock': {
            'python_path': os.path.dirname(os.path.abspath(__file__)),
            'bins': bins.tolist(), 'truth_log10Om': y.tolist(),
            'sigma_dex': float(sigma_dex)}},
        'params': {
            'log10r': {'prior': {'min': -5.0, 'max': 0.0},
                       'ref': FIDUCIAL['log10r'], 'proposal': 0.2,
                       'latex': r'\log_{10}r'},
            'r': {'value': 'lambda log10r: np.power(10., log10r)',
                  'latex': r'r'},
            'n_t': {'prior': {'min': -0.5, 'max': 0.5},
                    'ref': FIDUCIAL['n_t'], 'proposal': 0.03,
                    'latex': r'n_\mathrm{t}'},
            'Omega_bh2': {'value': 0.0223828},
            'Omega_ch2': {'value': 0.1201075},
            'H0': {'value': 67.32117},
            'DN_eff': {'value': 0.0},
            'A_s': {'value': 2.100549e-9},
            'cr': {'value': 1.0},
            'T_re': {'value': 2000.0, 'latex': r'T_\mathrm{re}'},
            'DN_re': {'value': 0.0, 'latex': r'\Delta N_\mathrm{re}'},
            'kappa10': {'value': 1e-2, 'latex': r'\kappa_{10}'}},
        'sampler': {'mcmc': {'max_samples': int(samples), 'seed': int(seed),
                             'proposal_scale': 2.0,
                             'Rminus1_stop': None,
                             'max_tries': 1e6}}}
    return info


def phase_chain(args):
    truth = _load_truth(args.out)
    info = build_info(truth, args.samples, sigma_dex=args.sigma_dex,
                      seed=args.seed)
    from cobaya.run import run
    seed_dir = 'fast_production_s%d' % int(args.seed)
    chain_out = os.path.join(out_dir(args), 'chains', seed_dir)
    t0 = time.time()
    updated, sampler = run(info, output=chain_out, force=True, resume=False)
    products = sampler.products()
    wall_s = time.time() - t0
    sample = products['sample'].to_numpy()
    sample_params = list(products['sample'].columns)
    stats = MC.chain_stats(sample, sample_params)
    engine_stats = None
    model = getattr(sampler, 'model', None)
    theories = getattr(model, 'theory', {}) if model is not None else {}
    for theory in (theories.values() if hasattr(theories, 'values') else theories):
        if hasattr(theory, 'engine_stats'):
            engine_stats = theory.engine_stats
            break
    rec = {'seed': int(args.seed), 'samples': int(len(sample)),
           'wall_s': wall_s,
           'per_sample_s': wall_s / max(len(sample), 1),
           'sigma_dex': args.sigma_dex,
           'stats': stats, 'engine_stats': engine_stats}
    np.savez(os.path.join(out_dir(args), 'fast_chain_s%d.npz' % int(args.seed)),
             sample=sample, sample_params=np.asarray(sample_params))
    with open(os.path.join(out_dir(args),
                           'fast_chain_s%d.json' % int(args.seed)), 'w',
              encoding='utf-8') as fh:
        json.dump(rec, fh, ensure_ascii=False, indent=1, default=float)
    print('chain done: seed=%d n=%d wall=%.0fs (%.2fs/sample)' %
          (int(args.seed), len(sample), wall_s,
           wall_s / max(len(sample), 1)), flush=True)
    print('min_ess=%.0f  stats=%s' % (
        stats['min_ess'], {k: v['mean'] for k, v in stats['stats'].items()
                           if v}), flush=True)
    if engine_stats:
        print('engine_stats=%s' % json.dumps(engine_stats), flush=True)


def fast_log10_at_bins(kw, bins):
    """Production fast solve; returns log10OmegaGW at bins + full-grid dn."""
    from stiffgwpy import fast_sgwb as FS
    from stiffgwpy.stiff_SGWB import LCDM_SG
    cfg = FS.apply_accuracy_mode('production')
    m = LCDM_SG(**kw)
    res = FS.SGWB_iter_fast(m, tol=cfg['tol'], transition_refine=True,
                            freq_grid='grid_independent', freq_res=2.0)
    if res is None or not getattr(m, 'SGWB_converge', False):
        return None, None
    o = np.argsort(np.asarray(m.f, dtype=float))
    fs = np.asarray(m.f, dtype=float)[o]
    lo = np.asarray(m.log10OmegaGW, dtype=float)[o]
    lo = lo[np.isfinite(lo)]
    fs = fs[np.isfinite(fs)]
    bins = np.asarray(bins, dtype=float)
    if fs.size < 4 or bins.min() < fs.min() or bins.max() > fs.max():
        return None, None
    spl = interpolate.CubicSpline(fs, lo)
    return spl(bins), float(m.cosmo_param['DN_eff'])


def ref_log10_at_bins(kw, bins, dn, rtol=1e-9, z_tail=8.0):
    """Continuous-sigma DOP853 reference at the SAME bins and matched DN."""
    from stiffgwpy import reference as REF
    from stiffgwpy.stiff_SGWB import LCDM_SG
    m = LCDM_SG(**kw)
    Ogw, Oj, Opgw, used = REF.spectrum_reference(m, np.asarray(bins, dtype=float),
                                                 dn, z_tail=z_tail, rtol=rtol,
                                                 workers=8)
    return np.log10(np.maximum(np.asarray(Ogw) - np.asarray(Oj), 1e-300))


def mock_loglike(lo_test, truth, sigma_dex, n_bad=-1e10):
    """Diagonal Gaussian in dex on the 11 bins (already in log10 Omega)."""
    lo_test = np.asarray(lo_test, dtype=float)
    truth = np.asarray(truth['truth_log10Om'], dtype=float)
    if not np.all(np.isfinite(lo_test)):
        return float('-inf')
    d = (lo_test - truth) / float(sigma_dex)
    return float(-0.5 * np.sum(d * d))


def phase_pointwise(args):
    truth = _load_truth(args.out)
    chain_path = os.path.join(args.out, 'fast_chain.npz')
    if not os.path.exists(chain_path):
        raise SystemExit('fast chain missing: run --phase chain first')
    z = np.load(chain_path)
    sample = z['sample']
    sample_params = list(z['sample_params'])
    params, specials = MC.split_sample_columns(sample_params)
    ind = MC.param_indices(sample_params)
    mlpost = sample[:, specials['minuslogpost']].astype(float)
    ok = np.isfinite(mlpost)
    # burn-in: keep the last 60% (Cobaya already drops the burn-in chunk for
    # single chains only if the sampler converged; be explicit here).
    start = int(0.4 * ok.sum())
    pool = np.nonzero(ok)[0][start:]
    rng = np.random.default_rng(args.seed + 1)
    idx = rng.choice(pool, size=min(args.k, pool.size), replace=False)
    # Only evaluate bulk points (the posterior support that matters).
    logpost = mlpost[idx]
    logpost = logpost - np.nanmax(logpost)
    keep = np.exp(logpost) > 0.001  # within ~6.9 log-likelihood of the mode
    idx = idx[keep]
    if idx.size < 8:
        raise SystemExit('too few bulk chain points (%d)' % idx.size)
    rows = sample[idx]
    bins = np.asarray(truth['bins'], dtype=float)
    sig = args.sigma_dex
    out = []
    t0 = time.perf_counter()
    for i, row in enumerate(rows):
        log10r = float(row[ind['log10r']])
        nt = float(row[ind['n_t']])
        kw = model_kwargs(log10r, nt)
        lo_f, dn = fast_log10_at_bins(kw, bins)
        rec = {'log10r': log10r, 'n_t': nt, 'status': 'ok'}
        if lo_f is None:
            rec['status'] = 'fast_failed'
            out.append(rec)
            continue
        ll_fast = mock_loglike(lo_f, truth, sig)
        try:
            lo_r = ref_log10_at_bins(kw, bins, dn)
        except Exception as exc:
            rec['status'] = 'ref_exception:%s' % type(exc).__name__
            out.append(rec)
            continue
        ll_ref = mock_loglike(lo_r, truth, sig)
        rec['ll_fast'] = ll_fast
        rec['ll_ref'] = ll_ref
        rec['dll'] = ll_ref - ll_fast          # logL(ref) - logL(fast)
        rec['dex_max_fast_vs_ref'] = float(
            np.max(np.abs(lo_f - lo_r)))
        out.append(rec)
        print('pt %3d log10r=%+.3f n_t=%+.3f dll=%+.3e dexmax=%.2e' %
              (i + 1, log10r, nt, rec['dll'], rec['dex_max_fast_vs_ref']),
              flush=True)
    rec_all = {'k': len(out), 'wall_s': time.perf_counter() - t0,
               'seed': args.seed + 1, 'sigma_dex': sig,
               'burn_start_frac': 0.4, 'points': out}
    with open(os.path.join(args.out, 'pointwise_reference.json'), 'w',
              encoding='utf-8') as fh:
        json.dump(rec_all, fh, ensure_ascii=False, indent=1)
    print('pointwise done: %d points in %.0fs' % (len(out), rec_all['wall_s']),
          flush=True)


def _merge_chain_stats(out, chain_jsons):
    """Concatenate independent chains (rows already carry block weights)."""
    try:
        all_w = []
        all_x = {}
        first = None
        for jn in chain_jsons:
            cj = json.load(open(os.path.join(out, jn), encoding='utf-8'))
            if first is None:
                first = cj
            z = np.load(os.path.join(out, jn.replace('.json', '.npz')))
            sp = list(z['sample_params'])
            params, specials = MC.split_sample_columns(sp)
            ind = MC.param_indices(sp)
            samp = z['sample']
            w = samp[:, specials['weight']].astype(float)
            ml = samp[:, specials['minuslogpost']].astype(float)
            ok = np.isfinite(ml) & (w > 0)
            w = w[ok]
            all_w.append(w)
            for name in params:
                col = samp[ok, ind[name]].astype(float)
                all_x.setdefault(name, []).append(col)
        W = np.concatenate(all_w)
        n = int(W.sum())
        out_stats = {'n_total': n, 'n_finite': n,
                     'failure_rate': 0.0, 'min_ess': float('inf'),
                     'params': first['stats']['params'],
                     'covariance_params': first['stats']['covariance_params'],
                     'covariance': None, 'map': first['map']}
        merged = {}
        for name, cols in all_x.items():
            x = np.concatenate(cols)
            wsum = float(W.sum())
            mean = float(np.sum(W * x) / wsum)
            std = float(np.sqrt(np.sum(W * (x - mean) ** 2) / wsum))
            q16, q50, q84 = MC.weighted_quantile(x, [0.16, 0.5, 0.84], W)
            ess = (float(x.size) if std <= max(1e-15, abs(mean) * 1e-12)
                   else MC.effective_sample_size(x, W))
            merged[name] = {'mean': mean, 'std': std, 'median': float(q50),
                            'p16': float(q16), 'p84': float(q84), 'ess': ess}
        out_stats['stats'] = merged
        out_stats['min_ess'] = min((v['ess'] for v in merged.values()
                                    if v is not None), default=0.0)
        return out_stats
    except Exception as exc:
        print('merge failed: %s' % exc, flush=True)
        return None


def _bootstrap_shift(shifts, n=2000, seed=42):
    rng = np.random.default_rng(seed)
    means = []
    x = np.asarray(shifts, dtype=float)
    for _ in range(n):
        b = rng.choice(x, size=x.size, replace=True)
        means.append(np.mean(b))
    return float(np.mean(means)), float(np.std(means))


def phase_report(args):
    _load_truth(args.out)  # verify the artifact exists (chain files checked below)
    chain_jsons = sorted(f for f in os.listdir(args.out)
                         if f.startswith('fast_chain_s') and f.endswith('.json'))
    if not chain_jsons:
        raise SystemExit('no fast chain found: run --phase chain first')
    chain_path = os.path.join(args.out, chain_jsons[0])
    if not os.path.exists(chain_path):
        raise SystemExit('fast chain missing: run --phase chain first')
    chain = json.load(open(chain_path, encoding='utf-8'))
    stats = chain['stats']
    if len(chain_jsons) > 1:
        merged = _merge_chain_stats(args.out, chain_jsons)
        if merged is not None:
            stats = merged
    # ---- measured per-point Delta logL from the reference legs ----
    pw = None
    pw_path = os.path.join(args.out, 'pointwise_reference.json')
    if os.path.exists(pw_path):
        pw = json.load(open(pw_path, encoding='utf-8'))
    n_total = sum(int(json.load(open(os.path.join(args.out, jn),
                                    encoding='utf-8'))['samples'])
                  for jn in chain_jsons)
    report = {'sigma_dex': args.sigma_dex, 'n_chains': len(chain_jsons),
              'fast_chain': {'samples': n_total,
                             'wall_s': chain['wall_s'],
                             'min_ess': stats['min_ess']}}
    if chain.get('engine_stats'):
        es = chain['engine_stats']
        report['telemetry'] = {'eval_status_counts': es['eval_status_counts'],
                               'last_eval_status': es['last_eval_status'],
                               'fallback_fraction': es['fallback_fraction'],
                               'escalation_fraction': es['escalation_fraction']}
    verdicts = {}
    if pw and pw['points']:
        dll = np.asarray([p['dll'] for p in pw['points']
                          if p['status'] == 'ok' and 'dll' in p])
        dex = np.asarray([p['dex_max_fast_vs_ref'] for p in pw['points']
                          if p['status'] == 'ok' and 'dex_max_fast_vs_ref' in p])
        dll_stats = {'n': int(dll.size), 'max': float(dll.max()),
                     'min': float(dll.min()),
                     'max_abs': float(np.max(np.abs(dll))),
                     'p95_abs': float(np.percentile(np.abs(dll), 95)),
                     'mean': float(dll.mean())}
        report['dll_fast_vs_reference'] = dll_stats
        report['dex_max_fast_vs_ref'] = {
            'max': float(dex.max()), 'p95': float(np.percentile(dex, 95)),
            'n': int(dex.size)}
        # importance-reweighted posterior shift, per sampled parameter
        points = [p for p in pw['points'] if p['status'] == 'ok' and 'dll' in p]
        w = np.exp(np.asarray([p['dll'] for p in points], dtype=float))
        w = w / w.sum()
        shifts = {}
        ess_reweight = float(1.0 / np.sum(w * w))
        for name in ('log10r', 'n_t'):
            x = np.asarray([p[name] for p in points], dtype=float)
            fstd = stats['stats'].get(name, {}).get('std')
            fmean = stats['stats'].get(name, {}).get('mean')
            if not fstd:
                continue
            wmean = float(np.sum(w * x))
            # plain (unweighted) subset mean, for comparison
            pmean = float(np.mean(x))
            shifts[name] = {'mean_fast_full': fmean,
                            'mean_subset': pmean,
                            'mean_reweighted': wmean,
                            'shift_sigma': (wmean - fmean) / fstd,
                            'subset_sigma': (pmean - fmean) / fstd}
        report['posterior_shift'] = shifts
        report['ess_reweighted'] = ess_reweight
        # per-sigma strength error surface (reuse the same measured dll and
        # rescale to what the stronger mock would have produced)
        report['error_surface'] = {}
        for sig in (0.10, 0.05, 0.02, 0.01):
            factor = (args.sigma_dex / sig) ** 2
            dll_s = dll * factor
            ws = np.exp(dll_s - dll_s.max())
            ws = ws / ws.sum()
            report['error_surface']['sigma_dex_%.2f' % sig] = {
                'max_abs_dll': float(np.max(np.abs(dll_s))),
                'ess_kish': float(1.0 / np.sum(ws * ws))}
        ok = np.isfinite(dll)
        verdicts['dll_max_abs_lt_0.1'] = bool(np.max(np.abs(dll[ok])) < 0.1)
        verdicts['posterior_shift_lt_0.1sigma'] = bool(
            all(abs(float(v['shift_sigma'])) < 0.1
                for v in report['posterior_shift'].values()))
    report['verdicts'] = verdicts
    with open(os.path.join(args.out, 'posterior_report.json'), 'w',
              encoding='utf-8') as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1, default=float)
    print(json.dumps(report, ensure_ascii=False, indent=1), flush=True)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--phase', choices=('truth', 'chain', 'pointwise',
                                        'report'), required=True)
    ap.add_argument('--out', default=os.path.join(REPO, 'docs', 'mcmc_posterior'))
    ap.add_argument('--samples', type=int, default=4000)
    ap.add_argument('--sigma-dex', type=float, default=0.05)
    ap.add_argument('--k', type=int, default=60)
    ap.add_argument('--seed', type=int, default=SEED, help='RNG seed for a chain run')
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    if args.phase == 'truth':
        phase_truth(args)
    elif args.phase == 'chain':
        phase_chain(args)
    elif args.phase == 'pointwise':
        phase_pointwise(args)
    else:
        phase_report(args)


if __name__ == '__main__':
    main()
