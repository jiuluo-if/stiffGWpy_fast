# -*- coding: utf-8 -*-
"""
plot_param_sweep.py -- audit phase 3 deliverables from the parameter sweep.

Reads ``docs/paramsweep/sweep_phase3.jsonl`` (written by ``param_sweep.py``)
and produces::

    docs/paramsweep/summary.json          totals + error distribution stats
    docs/paramsweep/worst_cases.json      worst-20 points (by DN_gw and dex)
    docs/paramsweep/error_distribution.png
    docs/paramsweep/parameter_error_map/  scatter of error vs each parameter

Usage::

    python scripts/plot_param_sweep.py --out docs/paramsweep
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib  # noqa: E402
from param_sweep import PARAMS  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

PERCENTILES = (50, 90, 95, 99, 99.9)


def pct(a, q):
    a = np.asarray([float(x) for x in a if x is not None and np.isfinite(x)])
    if not a.size:
        return None
    return float(np.percentile(a, q))


def load_records(out_dir):
    path = os.path.join(out_dir, 'sweep_phase3.jsonl')
    recs = []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if line:
                recs.append(json.loads(line))
    return recs


def spearman(x, y):
    """Simple rank correlation (no scipy dependency for the report)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 3 or np.all(x == x[0]):
        return None
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    dx = rx - rx.mean()
    dy = ry - ry.mean()
    denom = np.sqrt((dx * dx).sum() * (dy * dy).sum())
    return float((dx * dy).sum() / denom) if denom else None


def stats(a):
    a = np.asarray([float(x) for x in a if x is not None and np.isfinite(x)])
    if not a.size:
        return None
    return {'n': int(a.size), 'min': float(a.min()), 'max': float(a.max()),
            'mean': float(a.mean()), 'median': float(np.median(a)),
            'pct': {str(q): pct(a, q) for q in PERCENTILES}}


def main(argv=None):
    ap = argparse.ArgumentParser(description='phase-3 sweep deliverables')
    ap.add_argument('--out', default='docs/paramsweep')
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    recs = load_records(args.out)
    if not recs:
        print('no records found')
        return 1
    ok = [r for r in recs if r.get('status') == 'ok']
    dn_rel = np.array([r['DN_gw_last_rel'] for r in ok
                       if r.get('DN_gw_last_rel') is not None], float)
    dex = np.array([r['dex_max'] for r in ok if r.get('dex_max') is not None], float)
    lin = np.array([r['lin_max'] for r in ok if r.get('lin_max') is not None], float)
    status_counts = {}
    for r in recs:
        status_counts[r.get('status')] = status_counts.get(r.get('status'), 0) + 1

    summary = {
        'n_points': len(recs),
        'n_ok': len(ok),
        'status_counts': status_counts,
        'DN_gw_last_rel': stats(dn_rel),
        'dex_max': stats(dex),
        'lin_max': stats(lin),
        'lsoda_t_s': stats([r.get('t_s_lsoda') for r in recs
                            if r.get('t_s_lsoda') is not None]),
        'fast_t_s': stats([r.get('t_s_fast') for r in recs
                           if r.get('t_s_fast') is not None]),
        'iters_lsoda': stats([r.get('iters_lsoda') for r in ok
                              if r.get('iters_lsoda') is not None]),
        'iters_fast': stats([r.get('iters_fast') for r in ok
                             if r.get('iters_fast') is not None]),
    }

    # ---- worst cases ----
    worst_dn = sorted(ok, key=lambda r: r.get('DN_gw_last_rel') or -1,
                      reverse=True)[:20]
    worst_dex = sorted(ok, key=lambda r: r.get('dex_max') or -1, reverse=True)[:20]
    worst = {'by_DN_gw_last_rel': worst_dn, 'by_dex_max': worst_dex}
    with open(os.path.join(args.out, 'worst_cases.json'), 'w', encoding='utf-8') as fh:
        json.dump(worst, fh, ensure_ascii=False, indent=1)
    summary['worst_by_DN_gw_last_rel'] = [
        {'id': r['id'], 'DN_gw_last_rel': r.get('DN_gw_last_rel'),
         'dex_max': r.get('dex_max'), 'lin_max': r.get('lin_max'),
         'params': r['params']} for r in worst_dn]
    summary['worst_by_dex_max'] = [
        {'id': r['id'], 'DN_gw_last_rel': r.get('DN_gw_last_rel'),
         'dex_max': r.get('dex_max'), 'lin_max': r.get('lin_max'),
         'params': r['params']} for r in worst_dex]

    # ---- parameter correlations ----
    err_log = np.log10(np.maximum(dn_rel, 1e-16))
    corr = {}
    for name, lo, hi, mode in PARAMS:
        x = np.array([r['params'][name] for r in ok], float)
        corr[name] = spearman(x, err_log)
    summary['spearman_log10DNrel_vs_param'] = corr

    with open(os.path.join(args.out, 'summary.json'), 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)

    # ---- error distribution figure ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, data, label, xlab in (
            (axes[0], dn_rel, 'DN_gw last rel.', 'log10 rel. err'),
            (axes[1], dex, 'dex max', 'log10 dex'),
            (axes[2], lin, 'lin Omega max', 'log10 lin rel. err')):
        v = np.log10(np.maximum(np.asarray(data, float), 1e-16))
        ax.hist(v, bins=40, color='steelblue', alpha=0.75)
        ax.axvline(np.median(v), color='crimson', ls='--', label='median')
        ax.axvline(np.percentile(v, 95), color='darkorange', ls='--', label='p95')
        ax.set_xlabel(xlab)
        ax.set_ylabel('count')
        ax.set_title(label)
        ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'error_distribution.png'), dpi=140)
    plt.close(fig)

    # ---- parameter error maps ----
    map_dir = os.path.join(args.out, 'parameter_error_map')
    os.makedirs(map_dir, exist_ok=True)
    for name, lo, hi, mode in PARAMS:
        x = np.array([r['params'][name] for r in ok], float)
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.scatter(x, err_log, s=12, alpha=0.55, color='steelblue')
        ax.set_xlabel(name)
        ax.set_ylabel('log10 DN_gw last rel. err')
        ax.set_title('error vs %s (spearman %.2f)' %
                     (name, corr[name] if corr[name] is not None else float('nan')))
        if mode == 'log':
            ax.set_xscale('log')
        fig.tight_layout()
        fig.savefig(os.path.join(map_dir, 'map_%s.png' % name), dpi=130)
        plt.close(fig)

    print('summary: n=%d ok=%d status=%s' %
          (len(recs), len(ok), status_counts))
    print('DN_gw_last_rel: median %.3e  p95 %.3e  max %.3e' %
          (float(np.median(dn_rel)), float(np.percentile(dn_rel, 95)),
           float(dn_rel.max())))
    print('dex_max: median %.3e  p95 %.3e  max %.3e' %
          (float(np.median(dex)), float(np.percentile(dex, 95)),
           float(dex.max())))
    print('wrote %s' % os.path.join(args.out, 'summary.json'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
