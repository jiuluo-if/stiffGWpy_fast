# -*- coding: utf-8 -*-
"""
mcmc_compare.py -- Cobaya MCMC comparison: LSODA reference chain vs fast chain.

Runs two Cobaya MCMC chains with the same yaml run-info and the same seed --
one with ``engine: lsoda``, one with ``engine: fast`` (optionally a named
accuracy mode) -- then:

  * compares posterior mean / std / 16-84% CI / MAP per parameter,
  * evaluates *both* engines at a thinned subset of the chain points and
    reports the per-point ``Delta logL`` distribution (fast - lsoda),
  * reports posterior shift ``(mean_fast - mean_lsoda)/std_lsoda``,
  * reports failure rate (non-finite logpost) and wall-clock speedup.

Usage (run from the repo root so ``python_path`` entries resolve):

    python scripts/mcmc_compare.py --samples 2000 --fast-mode production --n-eval 50
    python scripts/mcmc_compare.py --skip-lsoda-chain --n-eval 200

Outputs ``mcmc_compare.jsonl`` under ``--out`` (default ``docs/mcmc``) plus the
raw Cobaya chain output under ``<out>/chains/<engine>/``.

Cost note: the per-point LSODA evaluation is the expensive part (tens of
seconds per point on this host).  Start with ``--n-eval 50``; the chain runs
themselves are cheap only when ``engine: fast`` is used, the LSODA chain is
the long pole.  Set ``SGWB_POOL_SIZE`` to limit LSODA worker processes on
memory-constrained hosts.
"""

import argparse
import copy
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np

# ---------------------------------------------------------------------------
# Pure analysis helpers (importable and unit-testable without cobaya)
# ---------------------------------------------------------------------------

def split_sample_columns(sample_params):
    """Split a Cobaya chain column list into parameter columns and specials."""
    specials = {}
    params = []
    for i, name in enumerate(sample_params):
        if (name == 'weight' or name.startswith('minuslog') or
                name == 'chi2' or name.startswith('chi2__')):
            specials[name] = i
        else:
            params.append(name)
    return params, specials


def param_indices(sample_params):
    """Map each parameter name to its column index in the full chain array."""
    params, _ = split_sample_columns(sample_params)
    pset = set(params)
    return {name: i for i, name in enumerate(sample_params) if name in pset}


def _finite_rows(mlpost):
    return np.isfinite(mlpost)


def effective_sample_size(values, weights=None):
    """Estimate ESS from Kish weights and the initial-positive ACF sequence."""
    x = np.asarray(values, dtype=float)
    if weights is None:
        w = np.ones(x.size, dtype=float)
    else:
        w = np.asarray(weights, dtype=float)
    valid = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x = x[valid]
    w = w[valid]
    n = x.size
    if n < 3:
        return float(n)
    kish = float((w.sum() ** 2) / np.sum(w ** 2))
    x = x - x.mean()
    var = float(np.dot(x, x) / n)
    if var <= 0:
        return min(float(n), kish)
    ac = np.correlate(x, x, mode='full')[n - 1:] / (var * n)
    tau = 1.0
    for lag in range(1, n):
        rho = float(ac[lag])
        if rho <= 0:
            break
        tau += 2.0 * rho
    acf_ess = n / tau
    return float(max(1.0, min(float(n), kish, acf_ess)))


def weighted_quantile(values, quantiles, weights=None):
    """Weighted quantiles using the left-continuous empirical CDF."""
    x = np.asarray(values, dtype=float)
    if weights is None:
        return np.percentile(x, np.asarray(quantiles) * 100.0)
    w = np.asarray(weights, dtype=float)
    mask = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x = x[mask]
    w = w[mask]
    if not x.size:
        return np.full(np.asarray(quantiles).shape, np.nan)
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cdf = np.cumsum(w) / np.sum(w)
    q = np.asarray(quantiles, dtype=float)
    indices = np.searchsorted(cdf, q, side='left')
    return x[np.clip(indices, 0, x.size - 1)]


def chain_stats(sample, sample_params):
    """Posterior summary from a Cobaya chain array.

    Returns per-parameter mean/std/median/16-84% CI and the MAP row
    (parameter with the minimum finite minuslogpost), plus the fraction of
    samples with non-finite minuslogpost (cobaya rejects those points).
    """
    params, specials = split_sample_columns(sample_params)
    indices = param_indices(sample_params)
    mlpost = sample[:, specials['minuslogpost']].astype(float)
    ok = _finite_rows(mlpost)
    if 'weight' in specials:
        weights = sample[:, specials['weight']].astype(float)
    else:
        weights = np.ones(sample.shape[0], dtype=float)
    total = int(sample.shape[0])
    failure_rate = float(1.0 - ok.sum() / total) if total else 0.0
    stats = {}
    for name in params:
        col = sample[ok, indices[name]].astype(float)
        wcol = weights[ok]
        valid = np.isfinite(col) & np.isfinite(wcol) & (wcol > 0)
        col = col[valid]
        wcol = wcol[valid]
        if col.size:
            wsum = float(wcol.sum())
            mean = float(np.sum(wcol * col) / wsum)
            std = float(np.sqrt(np.sum(wcol * (col - mean) ** 2) / wsum))
            q16, q50, q84 = weighted_quantile(col, [0.16, 0.5, 0.84], wcol)
            # Constant derived parameters carry no mixing information.  Do
            # not let a degenerate weighted ACF/Kish calculation report ESS=1
            # and block every certification run.
            ess = (float(col.size) if std <= max(1e-15, abs(mean) * 1e-12)
                   else effective_sample_size(col, wcol))
            stats[name] = {
                'mean': mean,
                'std': std,
                'median': float(q50),
                'p16': float(q16),
                'p84': float(q84),
                'ess': ess,
            }
        else:
            stats[name] = None
    map_row = {}
    if ok.any():
        i_map = int(np.nanargmin(mlpost[ok]))
        row = sample[ok, :][i_map]
        for name in params:
            map_row[name] = float(row[indices[name]])
        map_row['minuslogpost'] = float(mlpost[ok][i_map])
    finite_matrix = (sample[ok][:, [indices[name] for name in params]].astype(float)
                     if ok.any() else np.empty((0, len(params))))
    finite_weights = weights[ok] if ok.any() else np.empty(0)
    if finite_matrix.size:
        valid = (np.all(np.isfinite(finite_matrix), axis=1) &
                 np.isfinite(finite_weights) & (finite_weights > 0))
        finite_matrix = finite_matrix[valid]
        finite_weights = finite_weights[valid]
    if finite_matrix.shape[0] > 1 and finite_weights.sum() > 0:
        wn = finite_weights / finite_weights.sum()
        centered = finite_matrix - np.sum(wn[:, None] * finite_matrix, axis=0)
        covariance = centered.T @ (centered * wn[:, None])
    else:
        covariance = np.empty((len(params), len(params)))
    min_ess = min((v['ess'] for v in stats.values() if v is not None),
                  default=0.0)
    return {'params': params, 'n_finite': int(ok.sum()), 'n_total': total,
            'failure_rate': failure_rate, 'stats': stats, 'map': map_row,
            'covariance_params': params,
            'covariance': np.asarray(covariance).tolist(),
            'min_ess': float(min_ess)}


def posterior_shift(lsoda, fast):
    """Per-parameter posterior shift (mean_fast - mean_lsoda)/std_lsoda."""
    shift = {}
    for name, s in lsoda['stats'].items():
        if s is None or fast['stats'].get(name) is None:
            continue
        std = s['std']
        if std > 0:
            shift[name] = float((fast['stats'][name]['mean'] - s['mean']) / std)
    return shift


def max_abs_posterior_shift(shift):
    """Return the largest finite absolute posterior shift."""
    values = [abs(float(v)) for v in (shift or {}).values()
              if np.isfinite(v)]
    return max(values, default=float('inf'))


def posterior_distances(lsoda_sample, fast_sample, sample_params,
                        fast_sample_params=None):
    """Compute KS, Wasserstein and histogram-KL distances per parameter.

    These are descriptive diagnostics; they do not replace an effective
    sample-size/convergence assessment.  The KL estimate uses a shared,
    finite 32-bin histogram with a small pseudocount to remain defined when a
    chain has an empty bin.
    """
    from scipy.stats import ks_2samp, wasserstein_distance

    fast_sample_params = fast_sample_params or sample_params
    params_l, specials_l = split_sample_columns(sample_params)
    params_f, specials_f = split_sample_columns(fast_sample_params)
    params = [name for name in params_l if name in set(params_f)]
    ok_l = np.isfinite(lsoda_sample[:, specials_l['minuslogpost']])
    ok_f = np.isfinite(fast_sample[:, specials_f['minuslogpost']])
    out = {}
    indices_l = param_indices(sample_params)
    indices_f = param_indices(fast_sample_params)
    for name in params:
        x = np.asarray(lsoda_sample[ok_l, indices_l[name]], dtype=float)
        y = np.asarray(fast_sample[ok_f, indices_f[name]], dtype=float)
        x = x[np.isfinite(x)]
        y = y[np.isfinite(y)]
        if not x.size or not y.size:
            continue
        ks = ks_2samp(x, y, method='auto')
        lo = min(float(x.min()), float(y.min()))
        hi = max(float(x.max()), float(y.max()))
        if hi <= lo:
            kl = 0.0
        else:
            edges = np.linspace(lo, hi, 33)
            p = np.histogram(x, edges)[0].astype(float) + 1e-12
            q = np.histogram(y, edges)[0].astype(float) + 1e-12
            p /= p.sum()
            q /= q.sum()
            kl = float(np.sum(p * np.log(p / q)))
        out[name] = {'ks_stat': float(ks.statistic),
                     'ks_pvalue': float(ks.pvalue),
                     'wasserstein': float(wasserstein_distance(x, y)),
                     'kl_hist_32': kl}
    return out


def delta_logl_stats(delta):
    """Summary of per-point ``logL_fast - logL_lsoda`` (positive = fast better)."""
    delta = np.asarray(delta, dtype=float)
    delta = delta[np.isfinite(delta)]
    if delta.size == 0:
        return {'n': 0}
    abs_d = np.abs(delta)
    return {
        'n': int(delta.size),
        'max': float(delta.max()),
        'max_abs': float(abs_d.max()),
        'median': float(np.median(delta)),
        'mean': float(delta.mean()),
        'p95_abs': float(np.percentile(abs_d, 95)),
        'frac_abs_gt_0.5': float((abs_d > 0.5).mean()),
        'frac_abs_gt_1.0': float((abs_d > 1.0).mean()),
    }


def thin_idx(n, k):
    """Evenly spaced sample of ``k`` indices from ``range(n)`` (k >= n -> all)."""
    n = int(n)
    k = int(min(max(k, 1), n))
    if k >= n:
        return np.arange(n)
    return np.linspace(0, n - 1, k).astype(int)


# ---------------------------------------------------------------------------
# Cobaya-dependent parts (imported lazily so --help / unit tests work without
# a cobaya installation)
# ---------------------------------------------------------------------------

def build_info(yaml_path, engine, accuracy_mode, fast_threads, samples, seed,
               fallback=True, initial_point=None, proposal_scale=None):
    """Load a Cobaya run-info yaml and override the stiffGW theory engine."""
    import yaml
    with open(yaml_path, encoding='utf-8') as f:
        info = yaml.safe_load(f)
    info = copy.deepcopy(info)
    info.setdefault('likelihood', {})
    info.pop('debug', None)
    info.pop('output', None)
    info.pop('resume', None)
    base = os.path.dirname(os.path.abspath(yaml_path))
    for section in ('theory', 'likelihood'):
        for cfg in (info.get(section) or {}).values():
            if not isinstance(cfg, dict):
                continue
            pp = cfg.get('python_path')
            if isinstance(pp, str):
                cfg['python_path'] = (pp if os.path.isabs(pp)
                                      else os.path.join(base, pp))
            elif isinstance(pp, list):
                cfg['python_path'] = [p if os.path.isabs(p)
                                      else os.path.join(base, p) for p in pp]
            if section == 'likelihood':
                for key in cfg:
                    if key == 'path' or key == 'CC_file' or key.endswith('_file'):
                        val = cfg[key]
                        if isinstance(val, str) and not os.path.isabs(val):
                            cfg[key] = os.path.join(base, val)
    theories = info.setdefault('theory', {})
    # Always exercise the installable qualified class.  Keeping the legacy
    # short key (``stiffGW``) makes Cobaya resolve a module-relative class and
    # can skip ``stiffGW.yaml`` defaults (h/col_step/... become undefined).
    th = theories.pop('stiffGW', None)
    if th is None:
        th = theories.pop('stiffgwpy.cobaya.stiffGW.stiffGW', {})
    info['theory'] = {
        'stiffgwpy.cobaya.stiffGW.stiffGW': th,
        **theories,
    }
    th = info['theory']['stiffgwpy.cobaya.stiffGW.stiffGW']
    # The class is importable from the installed package; a source-relative
    # python_path would mask failures in wheel/qualified-name resolution.
    th.pop('python_path', None)
    th['engine'] = engine
    th['fallback'] = bool(fallback)
    if engine == 'fast':
        if accuracy_mode:
            th['accuracy_mode'] = accuracy_mode
        if fast_threads:
            th['fast_threads'] = int(fast_threads)
    mcmc = info.setdefault('sampler', {}).setdefault('mcmc', {})
    mcmc['max_samples'] = int(samples)
    mcmc['seed'] = int(seed)
    if proposal_scale is not None:
        if not np.isfinite(float(proposal_scale)) or float(proposal_scale) <= 0:
            raise ValueError('proposal_scale must be a positive finite number')
        mcmc['proposal_scale'] = float(proposal_scale)
    if initial_point:
        params = info.get('params', {})
        for name, value in initial_point.items():
            if name not in params or not isinstance(params[name], dict):
                raise ValueError('initial-point parameter is not defined in run-info: %s' % name)
            if not np.isfinite(float(value)):
                raise ValueError('initial-point value must be finite: %s=%r' % (name, value))
            # Cobaya uses scalar refs as deterministic starting points while
            # leaving the actual prior unchanged.  Derived parameters are not
            # sampled and therefore must not be included in this mapping.
            if 'derived' in params[name] or 'value' in params[name]:
                raise ValueError('initial-point parameter must be sampled/free: %s' % name)
            params[name]['ref'] = float(value)
    return info


def run_chain(info, label, out_dir, samples):
    """Run one Cobaya MCMC chain and return its samples plus wall time."""
    from cobaya.run import run
    chain_out = os.path.join(out_dir, 'chains', label)
    t0 = time.time()
    _, sampler = run(info, output=chain_out, force=True, resume=False)
    products = sampler.products()
    wall_s = time.time() - t0
    engine_stats = None
    # Cobaya keeps the initialized theory instances on the sampler model;
    # expose our telemetry when available.  Keep this best-effort so old
    # Cobaya releases still produce the normal chain record.
    model = getattr(sampler, 'model', None)
    theories = getattr(model, 'theory', {}) if model is not None else {}
    for theory in (theories.values() if hasattr(theories, 'values') else theories):
        if hasattr(theory, 'engine_stats'):
            engine_stats = theory.engine_stats
            break
    return {
        'label': label,
        'wall_s': wall_s,
        'per_sample_s': wall_s / max(samples, 1),
        'sample': products['sample'].to_numpy(),
        'sample_params': list(products['sample'].columns),
        'engine_stats': engine_stats,
    }


def _point_dict(model, row, indices):
    """Parameter dict for one chain row: only sampled params the model can set."""
    d = {name: row[i] for name, i in indices.items()}
    try:
        sampled = set(model.parameterization.sampled_params())
    except Exception:
        sampled = None
    if sampled:
        d = {k: v for k, v in d.items() if k in sampled}
    return d


def eval_minuslogpost_batch(model, rows, sample_params):
    """Minuslogpost of every point (inf on failure); also returns fail count."""
    indices = param_indices(sample_params)
    out = np.empty(len(rows))
    fails = 0
    for i, row in enumerate(rows):
        try:
            p = _point_dict(model, row, indices)
            ml = -float(model.logpost(p))  # minuslogpost (inf when the point fails)
            if not np.isfinite(ml):
                fails += 1
                ml = float('inf')
        except Exception:
            fails += 1
            ml = float('inf')
        out[i] = ml
    return out, fails


def compare_engines_at_points(yaml_path, fast_mode, fast_threads, points,
                              sample_params,
                              n_eval, engine_labels=('lsoda', 'fast')):
    """Evaluate both engines at thinned points; returns (delta, meta)."""
    from cobaya.model import get_model
    idx = thin_idx(len(points), n_eval)
    rows = np.asarray(points)[idx]
    results = {}
    for engine in engine_labels:
        info = build_info(yaml_path, engine, fast_mode, fast_threads,
                          samples=1, seed=1)
        model = get_model(info)
        t0 = time.time()
        ml, fails = eval_minuslogpost_batch(model, rows, sample_params)
        results[engine] = {'minuslogpost': ml, 'failures': fails,
                           'wall_s': time.time() - t0}
    ml_l = results['lsoda']['minuslogpost']
    ml_f = results['fast']['minuslogpost']
    both_ok = np.isfinite(ml_l) & np.isfinite(ml_f)
    delta = np.full(len(idx), np.nan)
    delta[both_ok] = ml_l[both_ok] - ml_f[both_ok]  # logL_fast - logL_lsoda
    meta = {
        'n_eval': int(len(idx)),
        'n_both_ok': int(both_ok.sum()),
        'lsoda_failures': results['lsoda']['failures'],
        'fast_failures': results['fast']['failures'],
        'lsoda_eval_wall_s': results['lsoda']['wall_s'],
        'fast_eval_wall_s': results['fast']['wall_s'],
    }
    return delta, meta


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _git_sha():
    try:
        return subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], stderr=subprocess.DEVNULL,
            text=True).strip()
    except Exception:
        return 'unknown'


def _yaml_sha256(path):
    """Return a content hash so the two engine runs are auditable twins."""
    try:
        h = hashlib.sha256()
        with open(path, 'rb') as fh:
            for block in iter(lambda: fh.read(1024 * 1024), b''):
                h.update(block)
        return h.hexdigest()
    except OSError:
        return 'unavailable'


def _fmt(x, nd=3):
    return '%.*g' % (nd, x) if isinstance(x, float) else str(x)


def parse_args(argv=None):
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--yaml', default=os.path.join(
        repo, 'stiffgwpy', 'cobaya', 'mcmc_compare.yaml'))
    p.add_argument('--out', default=os.path.join(repo, 'docs', 'mcmc'))
    p.add_argument('--samples', type=int, default=2000)
    p.add_argument('--seed', type=int, default=20260830)
    p.add_argument('--initial-point', default=None,
                   help='JSON file with explicit sampled-parameter starting point; '
                        'applied identically to both chains')
    p.add_argument('--fast-mode', choices=['', 'reference', 'production',
                                           'ultra-fast'], default='production')
    p.add_argument('--fast-threads', type=int, default=0)
    p.add_argument('--n-eval', type=int, default=50,
                   help='number of chain points at which both engines are '
                        'evaluated for Delta logL (LSODA eval is expensive)')
    p.add_argument('--min-effective-samples', type=float, default=2000.0,
                   help='ESS threshold required before labeling posterior '
                        'equivalence as certified (default: 2000)')
    p.add_argument('--max-posterior-shift', type=float, default=0.1,
                   help='maximum allowed absolute per-parameter posterior '
                        'shift in sigma units for certification (default: 0.1)')
    p.add_argument('--proposal-scale', type=float, default=None,
                   help='optional common Cobaya MCMC proposal scale override '
                        '(applied identically to both engines)')
    p.add_argument('--skip-lsoda-chain', action='store_true',
                   help='skip the long LSODA chain; compare fast chain '
                        'against per-point LSODA evaluations only')
    p.add_argument('--skip-pointwise', action='store_true',
                   help='skip same-point LSODA/fast Delta logL evaluation '
                        '(useful for chain pilots; not a certification run)')
    p.add_argument('--no-fallback', action='store_true',
                   help='disable LSODA fallback for a diagnostic pilot; '
                        'never eligible for posterior certification')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if not os.path.exists(args.yaml):
        sys.exit('yaml not found: %s' % args.yaml)
    try:
        import cobaya  # noqa: F401
    except ImportError:
        sys.exit('cobaya is not installed; run `pip install .[cobaya]` '
                 'on a machine with a working Cobaya environment.')

    os.makedirs(args.out, exist_ok=True)
    initial_point = None
    if args.initial_point:
        with open(args.initial_point, encoding='utf-8') as f:
            initial_point = json.load(f)
        if not isinstance(initial_point, dict) or not initial_point:
            sys.exit('--initial-point must contain a non-empty JSON object')
    if args.max_posterior_shift <= 0:
        sys.exit('--max-posterior-shift must be positive')
    engines = {}
    if not args.skip_lsoda_chain:
        info = build_info(args.yaml, 'lsoda', '', 0, args.samples, args.seed,
                          fallback=not args.no_fallback,
                          initial_point=initial_point,
                          proposal_scale=args.proposal_scale)
        engines['lsoda'] = run_chain(info, 'lsoda', args.out, args.samples)
        print('[mcmc_compare] lsoda chain: %.1f s (%.2f s/sample)'
              % (engines['lsoda']['wall_s'], engines['lsoda']['per_sample_s']))
    info = build_info(args.yaml, 'fast', args.fast_mode, args.fast_threads,
                      args.samples, args.seed,
                      fallback=not args.no_fallback,
                      initial_point=initial_point,
                      proposal_scale=args.proposal_scale)
    engines['fast'] = run_chain(info, 'fast', args.out, args.samples)
    print('[mcmc_compare] fast chain: %.1f s (%.3f s/sample)'
          % (engines['fast']['wall_s'], engines['fast']['per_sample_s']))

    # Chain-based posterior summaries.
    stats = {}
    for label, eng in engines.items():
        stats[label] = chain_stats(eng['sample'], eng['sample_params'])
        print('[mcmc_compare] %s: %d finite / %d samples, failure rate %.2g'
              % (label, stats[label]['n_finite'], stats[label]['n_total'],
                 stats[label]['failure_rate']))

    # Per-point Delta logL on the fast chain (or lsoda chain if available).
    # A pilot may explicitly skip this expensive LSODA diagnostic; such a
    # record is never eligible for posterior certification.
    if args.skip_pointwise:
        dl = {'n': 0, 'skipped': True}
        eval_meta = {'n_eval': 0, 'n_both_ok': 0,
                     'lsoda_failures': 0, 'fast_failures': 0,
                     'lsoda_eval_wall_s': 0.0, 'fast_eval_wall_s': 0.0,
                     'skipped': True}
        print('[mcmc_compare] pointwise Delta logL evaluation skipped '
              '(pilot mode; not certification evidence)')
    else:
        ref_engine = 'lsoda' if 'lsoda' in engines else 'fast'
        points = engines[ref_engine]['sample']
        print('[mcmc_compare] evaluating both engines at %d thinned points '
              '(LSODA eval, may take a while)...' % args.n_eval)
        delta, eval_meta = compare_engines_at_points(
            args.yaml, args.fast_mode, args.fast_threads, points,
            engines[ref_engine]['sample_params'],
            args.n_eval)
        dl = delta_logl_stats(delta)

    shift = posterior_shift(stats['lsoda'], stats['fast']) if 'lsoda' in stats \
        else {}
    max_shift = max_abs_posterior_shift(shift)
    distances = posterior_distances(
        engines['lsoda']['sample'], engines['fast']['sample'],
        engines['lsoda']['sample_params'],
        engines['fast']['sample_params']) if 'lsoda' in engines else {}
    speedup = {'mcmc_wall': None}
    if 'lsoda' in engines and engines['lsoda']['wall_s'] > 0:
        speedup['mcmc_wall'] = engines['lsoda']['wall_s'] / \
            engines['fast']['wall_s']
        speedup['per_sample'] = engines['lsoda']['per_sample_s'] / \
            engines['fast']['per_sample_s']
    if eval_meta['fast_eval_wall_s'] > 0 and \
            eval_meta['lsoda_eval_wall_s'] > 0:
        speedup['per_point_eval'] = eval_meta['lsoda_eval_wall_s'] / \
            eval_meta['fast_eval_wall_s']

    if 'lsoda' in stats:
        min_ess = min(stats['lsoda']['min_ess'], stats['fast']['min_ess'])
        posterior_certified = bool(
            min_ess >= args.min_effective_samples and
            stats['lsoda']['failure_rate'] == 0.0 and
            stats['fast']['failure_rate'] == 0.0 and
            not args.skip_pointwise and not args.no_fallback and
            eval_meta.get('n_both_ok', 0) > 0 and
            max_shift <= args.max_posterior_shift)
    else:
        min_ess = stats['fast']['min_ess']
        posterior_certified = False

    record = {
        'schema': 'mcmc_compare_v2',
        'date': time.strftime('%Y-%m-%dT%H:%M:%S'),
        'git_sha': _git_sha(),
        'yaml': os.path.basename(args.yaml),
        'yaml_sha256': _yaml_sha256(args.yaml),
        'samples': args.samples,
        'seed': args.seed,
        'initial_point': initial_point,
        'proposal_scale': args.proposal_scale,
        'fast_mode': args.fast_mode or None,
        'fast_threads': args.fast_threads,
        'fallback_enabled': not args.no_fallback,
        'n_eval': args.n_eval,
        'engines': stats,
        'engine_telemetry': {label: eng.get('engine_stats')
                             for label, eng in engines.items()},
        'posterior_shift': shift,
        'max_abs_posterior_shift': max_shift,
        'max_posterior_shift_threshold': args.max_posterior_shift,
        'posterior_distances': distances,
        'min_effective_samples': float(min_ess),
        'posterior_certified': posterior_certified,
        'delta_logL': dl,
        'eval': eval_meta,
        'speedup': speedup,
        # The comparison contract is intentionally explicit: both chains are
        # built from the same YAML/prior/initial-point policy, seed and MCMC
        # settings; only the theory engine knobs differ.
        'comparison_contract': {
            'same_yaml': True,
            'same_prior': True,
            'same_initial_point_policy': True,
            'same_sampler_settings': True,
            'same_seed': True,
            'engine_only_difference': True,
        },
    }
    out_path = os.path.join(args.out, 'mcmc_compare.jsonl')
    with open(out_path, 'a', encoding='utf-8') as f:
        f.write(json.dumps(record, allow_nan=True) + '\n')

    print('\n=== posterior shift (mean_fast - mean_lsoda)/std_lsoda ===')
    for name, s in sorted(shift.items()):
        print('  %-14s %s' % (name, _fmt(s)))
    print('=== Delta logL (fast - lsoda), n=%s ===' % dl.get('n'))
    for k in ('median', 'mean', 'max', 'max_abs', 'p95_abs',
              'frac_abs_gt_0.5', 'frac_abs_gt_1.0'):
        if k in dl:
            print('  %-14s %s' % (k, _fmt(dl[k])))
    if distances:
        print('=== posterior distances (KS / Wasserstein / histogram KL) ===')
        for name, d in sorted(distances.items()):
            print('  %-14s KS=%s W=%s KL=%s' %
                  (name, _fmt(d['ks_stat']), _fmt(d['wasserstein']),
                   _fmt(d['kl_hist_32'])))
    print('=== speedup ===')
    for k, v in speedup.items():
        print('  %-16s %s' % (k, _fmt(v)))
    print('=== engine telemetry ===')
    for label, eng in engines.items():
        telem = eng.get('engine_stats') or {}
        print('  %-8s fast_evals=%s fast_failures=%s failure_fraction=%s '
              'guard_rejections=%s physical_rejections=%s lsoda_evals=%s '
              'lsoda_fallbacks=%s fallback_fraction=%s' %
              (label, telem.get('fast_evals', 0),
               telem.get('fast_failures', 0),
               _fmt(telem.get('fast_failure_fraction', 0.0)),
               telem.get('fast_guard_rejections', 0),
               telem.get('fast_physical_rejections', 0),
               telem.get('lsoda_evals', 0),
               telem.get('lsoda_fallbacks', 0),
               _fmt(telem.get('fallback_fraction', 0.0))))
    print('=== posterior certification ===')
    print('  min ESS=%s threshold=%s max|shift|=%s threshold=%s certified=%s' %
          (_fmt(min_ess), _fmt(args.min_effective_samples),
           _fmt(max_shift), _fmt(args.max_posterior_shift),
           posterior_certified))
    print('\nrecord appended to %s' % out_path)


if __name__ == '__main__':
    main()
