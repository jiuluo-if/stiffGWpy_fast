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
    total = int(sample.shape[0])
    failure_rate = float(1.0 - ok.sum() / total) if total else 0.0
    stats = {}
    for name in params:
        col = sample[ok, indices[name]].astype(float)
        if col.size:
            stats[name] = {
                'mean': float(col.mean()),
                'std': float(col.std()),
                'median': float(np.percentile(col, 50)),
                'p16': float(np.percentile(col, 16)),
                'p84': float(np.percentile(col, 84)),
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
    return {'params': params, 'n_finite': int(ok.sum()), 'n_total': total,
            'failure_rate': failure_rate, 'stats': stats, 'map': map_row}


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

def build_info(yaml_path, engine, accuracy_mode, fast_threads, samples, seed):
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
    th = info.setdefault('theory', {}).setdefault('stiffGW', {})
    th['engine'] = engine
    th['fallback'] = True
    if engine == 'fast':
        if accuracy_mode:
            th['accuracy_mode'] = accuracy_mode
        if fast_threads:
            th['fast_threads'] = int(fast_threads)
    mcmc = info.setdefault('sampler', {}).setdefault('mcmc', {})
    mcmc['max_samples'] = int(samples)
    mcmc['seed'] = int(seed)
    return info


def run_chain(info, label, out_dir, samples):
    """Run one Cobaya MCMC chain and return its samples plus wall time."""
    from cobaya.run import run
    chain_out = os.path.join(out_dir, 'chains', label)
    t0 = time.time()
    _, sampler = run(info, output=chain_out, force=True, resume=False)
    products = sampler.products()
    wall_s = time.time() - t0
    return {
        'label': label,
        'wall_s': wall_s,
        'per_sample_s': wall_s / max(samples, 1),
        'sample': products['sample'].to_numpy(),
        'sample_params': list(products['sample'].columns),
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
    p.add_argument('--fast-mode', choices=['', 'reference', 'production',
                                           'ultra-fast'], default='production')
    p.add_argument('--fast-threads', type=int, default=0)
    p.add_argument('--n-eval', type=int, default=50,
                   help='number of chain points at which both engines are '
                        'evaluated for Delta logL (LSODA eval is expensive)')
    p.add_argument('--skip-lsoda-chain', action='store_true',
                   help='skip the long LSODA chain; compare fast chain '
                        'against per-point LSODA evaluations only')
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
    engines = {}
    if not args.skip_lsoda_chain:
        info = build_info(args.yaml, 'lsoda', '', 0, args.samples, args.seed)
        engines['lsoda'] = run_chain(info, 'lsoda', args.out, args.samples)
        print('[mcmc_compare] lsoda chain: %.1f s (%.2f s/sample)'
              % (engines['lsoda']['wall_s'], engines['lsoda']['per_sample_s']))
    info = build_info(args.yaml, 'fast', args.fast_mode, args.fast_threads,
                      args.samples, args.seed)
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

    record = {
        'schema': 'mcmc_compare_v1',
        'date': time.strftime('%Y-%m-%dT%H:%M:%S'),
        'git_sha': _git_sha(),
        'yaml': os.path.basename(args.yaml),
        'samples': args.samples,
        'seed': args.seed,
        'fast_mode': args.fast_mode or None,
        'fast_threads': args.fast_threads,
        'n_eval': args.n_eval,
        'engines': stats,
        'posterior_shift': shift,
        'delta_logL': dl,
        'eval': eval_meta,
        'speedup': speedup,
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
    print('=== speedup ===')
    for k, v in speedup.items():
        print('  %-16s %s' % (k, _fmt(v)))
    print('\nrecord appended to %s' % out_path)


if __name__ == '__main__':
    main()
