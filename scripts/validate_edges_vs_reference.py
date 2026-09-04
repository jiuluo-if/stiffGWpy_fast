# -*- coding: utf-8 -*-
"""validate_edges_vs_reference.py -- matched production-vs-reference spots.

Extends the Layer-A-style matched certification (fast production, grid-
independent freq grid, z_tail=8, rtol=1e-9) to parameter-axis edges and
transition-sensitive interiors.  LSODA is never a truth anchor: the oracle is
the continuous-sigma DOP853 reference in ``stiffgwpy_fast.reference``.

Phases (checkpointed; reruns skip completed labels):
  python scripts/validate_edges_vs_reference.py --phase reference [--pool 6]
  python scripts/validate_edges_vs_reference.py --phase summary

Outputs (default out dir ``docs/paramsweep_z8b``):
  reference_points.jsonl     per-point matched fast+reference records
  validation_summary.json    aggregate stats + gates (Layer-A conventions)
  validation_summary.md      human-readable table
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

_spec = importlib.util.spec_from_file_location(
    "validate_fast_vs_reference",
    os.path.join(ROOT, "scripts", "validate_fast_vs_reference.py"))
VFR = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(VFR)

Z_TAIL = 8.0
RTOL = 1e-9
FREQ_RES = 1.0

# Parameter-axis edges (u ~ 0.02 / 0.98 of the Sobol physical box) and
# transition-sensitive interiors; params subsets follow SINGLE_POINTS style.
EDGE_POINTS = {
    "edge_r_lo":      dict(r=1.26e-6, cr=1, T_re=2e3, kappa10=1e-2),
    "edge_r_hi":      dict(r=7.94e-2, cr=1, T_re=2e3, kappa10=1e-2),
    "edge_tre_lo":    dict(r=1e-2, cr=1, T_re=1.26e1, kappa10=1e-2),
    "edge_tre_hi":    dict(r=1e-2, cr=1, T_re=7.94e5, kappa10=1e-2),
    # DN_re is only physically active on the free-tilt branch (cr=0): under
    # the consistency relation (cr=1) the model overwrites DN_re from (r,
    # T_re, kappa10), so DN_re-axis probes must run with cr=0, n_t=0.
    "edge_dnre_lo":   dict(r=1e-2, cr=0, n_t=0.0, T_re=2e3, DN_re=0.6,
                           kappa10=1e-2),
    "edge_dnre_hi":   dict(r=1e-2, cr=0, n_t=0.0, T_re=2e3, DN_re=29.4,
                           kappa10=1e-2),
    "edge_kap_lo":    dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.3e-6),
    "edge_kap_hi":    dict(r=1e-2, cr=1, T_re=2e3, kappa10=0.76),
    "edge_nt_red":    dict(r=1e-2, cr=0, n_t=-0.48, T_re=2e3, kappa10=1e-2),
    "edge_nt_blue":   dict(r=1e-2, cr=0, n_t=0.48, T_re=2e3, kappa10=1e-2),
    "interior_r05":   dict(r=5e-2, cr=1, T_re=2e3, kappa10=1e-2),
    "interior_tre300": dict(r=3e-2, cr=1, T_re=3e2, kappa10=1e-1),
    "interior_ntred": dict(r=1e-2, cr=0, n_t=-0.3, T_re=2e3, kappa10=1e-2),
    "interior_dnre20": dict(r=3e-2, cr=0, n_t=0.0, T_re=2e3, DN_re=20.0,
                            kappa10=1e-2),
    "interior_kap03": dict(r=2e-2, cr=1, T_re=1e3, kappa10=0.3),
    "interior_tre1000_r5e3": dict(r=5e-3, cr=1, T_re=1e3, kappa10=5e-2),
}


def env_meta():
    meta = {'python': sys.version.split()[0], 'platform': sys.platform,
            'numpy': np.__version__}
    for mod in ('numba', 'scipy'):
        try:
            meta[mod] = __import__(mod).__version__
        except Exception:
            pass
    try:
        out = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True,
                             text=True, timeout=10)
        if out.returncode == 0:
            meta['commit'] = out.stdout.strip()
    except Exception:
        pass
    return meta


def _run_case(case):
    kw, label = case
    return VFR.reference_point(kw, z_tail=Z_TAIL, rtol=RTOL,
                               freq_res=FREQ_RES, label=label, kind='edge')


def run_reference(args):
    import multiprocessing as mp
    path = os.path.join(args.out, 'reference_points.jsonl')
    os.makedirs(args.out, exist_ok=True)
    cases = [(dict(kw), name) for name, kw in EDGE_POINTS.items()]
    seen = VFR.done_ids(path)
    todo = [c for c in cases if c[1] not in seen]
    if not todo:
        print('all %d edge points already recorded' % len(cases))
        return
    print('edges: %d done, %d to solve (pool=%d)' % (
        len(cases) - len(todo), len(todo), args.pool), flush=True)
    t_start = time.perf_counter()
    if args.pool > 1 and len(todo) > 1:
        ctx = mp.get_context('spawn')
        with ctx.Pool(args.pool) as pool:
            for rec in pool.imap_unordered(_run_case, todo):
                VFR.append_line(path, rec)
                print('done %-18s status=%-10s DN_rel=%+.3e sig_rel=%.2e '
                      'ref_dt=%.0fs elapsed=%.0fs' % (
                          rec.get('label'), rec.get('status'),
                          rec.get('DN_gw_rel', float('nan')),
                          ((rec.get('signal') or {}).get('rel') or {}).get(
                              'max', float('nan')),
                          rec.get('ref_dt') or 0.0,
                          time.perf_counter() - t_start), flush=True)
    else:
        for case in todo:
            rec = _run_case(case)
            VFR.append_line(path, rec)
            print('done %-18s status=%s' % (rec.get('label'),
                                            rec.get('status')), flush=True)
    print('edges reference phase finished in %.0fs' % (
        time.perf_counter() - t_start), flush=True)


def run_summary(args):
    path = os.path.join(args.out, 'reference_points.jsonl')
    recs = [r for r in VFR.load_lines(path) if r.get('status') == 'ok'
            and r.get('DN_gw_ref')]
    all_recs = VFR.load_lines(path)
    if not recs:
        print('no ok matched edge records; run --phase reference first')
        return 1
    rels = np.abs([r['DN_gw_rel'] for r in recs])
    summary = {'n_points': len(recs), 'n_total': len(all_recs),
               'meta': env_meta(),
               'settings': {'z_tail': Z_TAIL, 'rtol': RTOL,
                            'freq_res': FREQ_RES,
                            'grid': 'grid_independent (matched fast/reference)',
                            'engine': 'production fast vs continuous-sigma DOP853'},
               'DN_gw_rel_abs': {
                   'max': float(rels.max()),
                   'p95': float(np.percentile(rels, 95)),
                   'p99': float(np.percentile(rels, 99)),
                   'median': float(np.median(rels))},
               'acceptance': {}}
    for band in ('signal', 'transition'):
        dex = [r[band]['dex']['max'] for r in recs
               if (r.get(band) or {}).get('dex')]
        rel = [r[band]['rel']['max'] for r in recs
               if (r.get(band) or {}).get('rel')]
        summary[band] = {'dex_max': max(dex), 'dex_p95': float(
            np.percentile(dex, 95)), 'rel_max': max(rel),
            'rel_p95': float(np.percentile(rel, 95))}
        summary['acceptance']['%s_rel_lt_1e-3' % band] = bool(max(rel) < 1e-3)
    summary['acceptance']['DN_gw_rel_lt_1e-4'] = bool(rels.max() < 1e-4)
    with open(os.path.join(args.out, 'validation_summary.json'), 'w',
              encoding='utf-8') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    L = []
    L.append('# Extended matched fast-vs-reference spots (parameter edges + interiors)')
    L.append('')
    L.append('Production fast engine vs continuous-sigma DOP853 reference on the SAME '
             'grid-independent frequency grid, SAME `z_tail=8`, rtol=1e-9, freq_res=1.0. '
             'Points cover parameter-axis edges (u ~ 0.02/0.98 of the sampled box) and '
             'transition-sensitive interiors.')
    L.append('')
    L.append('| label | r | n_t | cr | T_re | DN_re | kappa10 | sig rel max | DN rel | status |')
    L.append('| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |')
    for r in sorted(recs, key=lambda x: x['label']):
        p = r['params']
        L.append('| %s | %.3e | %s | %s | %.3e | %s | %.3e | %.3e | %+.3e | ok |'
                 % (r['label'], p.get('r', float('nan')),
                    ('%.2f' % p['n_t']) if p.get('n_t') is not None else '-',
                    p.get('cr', '-'), p.get('T_re', float('nan')),
                    ('%.1f' % p['DN_re']) if p.get('DN_re') is not None else '-',
                    p.get('kappa10', float('nan')),
                    r['signal']['rel']['max'], r['DN_gw_rel']))
    L.append('')
    L.append('Aggregates over %d ok points: signal rel max **%.3e** (gate <1e-3: %s), '
             'transition rel max **%.3e** (gate <1e-3: %s), DN rel abs median **%.3e** / '
             'max **%.3e** (gate <1e-4: %s).'
             % (len(recs), summary['signal']['rel_max'],
                'PASS' if summary['acceptance']['signal_rel_lt_1e-3'] else 'FAIL',
                summary['transition']['rel_max'],
                'PASS' if summary['acceptance']['transition_rel_lt_1e-3'] else 'FAIL',
                summary['DN_gw_rel_abs']['median'],
                summary['DN_gw_rel_abs']['max'],
                'PASS' if summary['acceptance']['DN_gw_rel_lt_1e-4'] else 'FAIL'))
    n_guard = sum(1 for r in all_recs
                  if r.get('status') == 'fast_rejected'
                  and r.get('reason') == 'shared_Neff_guard')
    L.append('')
    L.append('Non-ok records: %d total (%d explicit shared_Neff_guard, others exception/'
             'rejected recorded per point; never silent).'
             % (len(all_recs) - len(recs), n_guard))
    with open(os.path.join(args.out, 'validation_summary.md'), 'w',
              encoding='utf-8') as fh:
        fh.write('\n'.join(L) + '\n')
    print('summary saved:', os.path.join(args.out, 'validation_summary.json'))
    print('signal rel max %.3e | transition rel max %.3e | DN rel abs med %.3e / max %.3e'
          % (summary['signal']['rel_max'], summary['transition']['rel_max'],
             summary['DN_gw_rel_abs']['median'], summary['DN_gw_rel_abs']['max']))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase', choices=['reference', 'summary'], default='reference')
    ap.add_argument('--out', default=os.path.join('docs', 'paramsweep_z8b'))
    ap.add_argument('--pool', type=int, default=6)
    args = ap.parse_args()
    if not os.path.isabs(args.out):
        args.out = os.path.join(ROOT, args.out)
    if args.phase == 'reference':
        run_reference(args)
    else:
        return run_summary(args)
    return 0


if __name__ == '__main__':
    sys.exit(main())
