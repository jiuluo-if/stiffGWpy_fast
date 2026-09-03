# -*- coding: utf-8 -*-
"""validate_plain_grid_vs_reference.py -- plain-grid (fast tier) vs oracle.

Physics-first, LSODA-free.  The continuous-sigma DOP853 reference
(``stiffgwpy.reference``) is the truth anchor; the plain-grid tier
(``accuracy_mode='fast'``: h=0.02 / col_step=8 / no transition_refine /
phase_max=0 / construct frequency grid) is validated on the SAME z_tail
(z_tail=8, matched to the certified Layer A runs) at its OWN frequency nodes:
the reference spectrum is evaluated exactly on the plain-grid nodes, so the
residual isolates the plain-grid engine error (ODE stepping + tail + no-kink
handling) without any frequency-grid interpolation.

Phases (checkpointed; reruns skip completed labels):
  python scripts/validate_plain_grid_vs_reference.py --phase plain [--pool 3] [--workers 3]
  python scripts/validate_plain_grid_vs_reference.py --phase summary [--out docs/paramsweep_plain]

Outputs (default out dir ``docs/paramsweep_plain``):
  plain_points.jsonl        per-point fast+reference matched records
  validation_summary.json   aggregate stats + gates
  validation_summary.md     human-readable report
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy import fast_sgwb as FS  # noqa: E402
from stiffgwpy import global_param as gp  # noqa: E402
from stiffgwpy import reference as REF  # noqa: E402
from stiffgwpy.stiff_SGWB import LCDM_SG  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "validate_fast_vs_reference",
    os.path.join(ROOT, "scripts", "validate_fast_vs_reference.py"))
_VFR = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_VFR)

ln10 = math.log(10.0)
Z_TAIL = 8.0
RTOL = 1e-9
FREQ_RES = 1.0


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


def plain_point(kw, label, z_tail=Z_TAIL, rtol=RTOL, freq_res=FREQ_RES):
    """Plain-grid fast solve + reference on the plain-grid frequency nodes."""
    rec = {'label': label, 'kind': 'single', 'params': dict(kw),
           'z_tail': z_tail, 'rtol': rtol, 'freq_res': freq_res,
           'mode': 'fast(plain-grid)', 'status': 'ok'}
    cfg = FS.apply_accuracy_mode('fast')       # h=0.02, col_step=8, no refine
    FS.set_z_tail(z_tail)                       # matched tail (isolate engine)
    m = LCDM_SG(**kw)
    t0 = time.perf_counter()
    try:
        res = FS.SGWB_iter_fast(m, tol=cfg['tol'], freq_res=freq_res,
                                transition_refine=False, sigma_exact=False)
    except Exception as exc:
        rec['status'] = 'fast_exception'
        rec['error'] = '%s: %s' % (type(exc).__name__, exc)
        return rec
    rec['fast_dt'] = time.perf_counter() - t0
    if res is None:
        rec['status'] = 'fast_rejected'
        rec['reason'] = getattr(m, 'fast_failure_reason', 'unknown')
        return rec
    dn_total = float(m.cosmo_param['DN_eff'])
    dn_gw_fast = float(np.asarray(m.DN_gw)[-1])
    order = np.argsort(np.asarray(m.f, dtype=float))
    f_fast = np.asarray(m.f, dtype=float)[order]
    lo_fast = np.asarray(m.log10OmegaGW, dtype=float)[order]
    G = f_fast.astype(np.float64)
    if G.size < 8:
        rec['status'] = 'fast_exception'
        rec['error'] = 'plain grid too small: %d nodes' % int(G.size)
        return rec
    t0 = time.perf_counter()
    try:
        Ogw, Oj, Opgw, used = REF.spectrum_reference(
            m, G, dn_total, z_tail=z_tail, rtol=rtol,
            workers=int(os.environ.get('PLAIN_REF_WORKERS', '3')))
    except Exception as exc:
        rec['status'] = 'ref_exception'
        rec['error'] = '%s: %s' % (type(exc).__name__, exc)
        return rec
    rec['ref_dt'] = time.perf_counter() - t0
    g2, qerr, ierr = REF.integrate_spectrum(G, Ogw, Oj)
    Omega_nu = gp.Omega_nh2 / m.derived_param['h'] ** 2
    dn_gw_ref = float(gp.Neff0 * g2 / Omega_nu)
    lo_ref = np.log10(np.maximum(Ogw - Oj, 1e-300))
    rec['DN_gw_fast'] = dn_gw_fast
    rec['DN_gw_ref'] = dn_gw_ref
    rec['DN_gw_rel'] = (dn_gw_fast - dn_gw_ref) / dn_gw_ref
    rec['quadrature_error'] = float(qerr) if qerr is not None else None
    rec['interpolation_error'] = float(ierr) if ierr is not None else None
    rec['used_tail_frac'] = float(np.mean(used))
    rec['n_freq'] = int(G.size)
    for bname, band in (('all', (-30.0, 10.0)),
                        ('signal', _VFR.SIGNAL_BAND),
                        ('transition', _VFR.TRANSITION_BAND)):
        rec[bname] = _VFR.band_stats(lo_ref, lo_fast, f_fast, band)
    return rec


def _run_case(case):
    kw, label = case
    return plain_point(kw, label)


def run_plain(args):
    import multiprocessing as mp
    path = os.path.join(args.out, 'plain_points.jsonl')
    os.makedirs(args.out, exist_ok=True)
    cases = [(dict(kw), name) for name, kw in _VFR.SINGLE_POINTS.items()]
    seen = _VFR.done_ids(path)
    todo = [c for c in cases if c[1] not in seen]
    if not todo:
        print('all %d plain-grid points already recorded in %s' % (len(cases), path))
        return
    print('plain phase: %d done, %d to solve (pool=%d)' % (
        len(cases) - len(todo), len(todo), args.pool), flush=True)
    t_start = time.perf_counter()
    if args.pool > 1 and len(todo) > 1:
        ctx = mp.get_context('spawn')
        with ctx.Pool(args.pool) as pool:
            for rec in pool.imap_unordered(_run_case, todo):
                with open(path, 'a', encoding='utf-8') as fh:
                    fh.write(json.dumps(rec, ensure_ascii=False) + '\n')
                print('done %-12s status=%-14s ref_dt=%.0fs  elapsed=%.0fs'
                      % (rec.get('label'), rec.get('status'),
                         rec.get('ref_dt') or 0.0,
                         time.perf_counter() - t_start), flush=True)
    else:
        for case in todo:
            rec = _run_case(case)
            with open(path, 'a', encoding='utf-8') as fh:
                fh.write(json.dumps(rec, ensure_ascii=False) + '\n')
            print('done %-12s status=%-14s ref_dt=%.0fs' % (
                rec.get('label'), rec.get('status'),
                rec.get('ref_dt') or 0.0), flush=True)
    print('plain phase finished in %.0fs' % (time.perf_counter() - t_start))


def _stats_abs(values):
    arr = np.abs(np.asarray(values, dtype=float))
    if arr.size == 0:
        return {}
    return {'n': int(arr.size), 'median': float(np.median(arr)),
            'p95': float(np.percentile(arr, 95)),
            'p99': float(np.percentile(arr, 99)),
            'max': float(np.max(arr))}


def _merge_bands(rows):
    out = {}
    for bname in ('all', 'signal', 'transition'):
        dex = np.concatenate([r[bname]['dex']['max'] and [r[bname]['dex']['max']] or []
                              for r in rows if r.get(bname)])
        rel = np.concatenate([r[bname]['rel']['max'] and [r[bname]['rel']['max']] or []
                              for r in rows if r.get(bname)])
        out[bname] = {'dex': _stats_abs(dex), 'rel': _stats_abs(rel)}
    return out


def run_summary(args):
    path = os.path.join(args.out, 'plain_points.jsonl')
    rows = [r for r in _VFR.load_lines(path) if r.get('status') == 'ok']
    all_rows = _VFR.load_lines(path)
    if not rows:
        print('no ok plain-grid records; run --phase plain first')
        return 1
    sig = _stats_abs([r['signal']['rel']['max'] for r in rows])
    tra = _stats_abs([r['transition']['rel']['max'] for r in rows])
    dn = _stats_abs([r['DN_gw_rel'] for r in rows])
    gate_sig = sig['max'] < 1e-3
    gate_tra = tra['max'] < 1e-3
    gate_dn = dn['median'] < 1e-4
    summary = {
        'n_points': len(rows), 'n_total': len(all_rows),
        'meta': env_meta(),
        'settings': {'z_tail': Z_TAIL, 'rtol': RTOL, 'freq_res': FREQ_RES,
                     'mode': 'fast(plain-grid)', 'reference': 'continuous-sigma DOP853'},
        'acceptance': {
            'signal_rel_lt_1e-3': gate_sig,
            'transition_rel_lt_1e-3': gate_tra,
            'DN_gw_rel_lt_1e-4': gate_dn,
        },
        'signal': sig, 'transition': tra, 'DN_gw_rel_abs': dn,
        'runtime_s': {'fast_median': float(np.median([r['fast_dt'] for r in rows])),
                      'ref_median': float(np.median([r['ref_dt'] for r in rows]))},
    }
    with open(os.path.join(args.out, 'validation_summary.json'), 'w',
              encoding='utf-8') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    L = []
    L.append('# Plain-grid tier (fast) vs continuous-sigma reference — matched z8')
    L.append('')
    L.append('Plain-grid engine (`accuracy_mode="fast"`: h=0.02, col_step=8, no '
             'transition_refine, phase_max=0, construct grid) validated against the '
             'independent continuous-sigma DOP853 reference (rtol=1e-9) on the SAME '
             '`z_tail=8` at the plain-grid\'s OWN frequency nodes — the residual '
             'isolates plain-grid engine error (no frequency-grid interpolation).')
    L.append('')
    L.append('| label | n_freq | signal rel max | transition rel max | DN_gw rel | status |')
    L.append('| --- | ---: | ---: | ---: | ---: | --- |')
    for r in rows:
        L.append('| %s | %d | %.3e | %.3e | %.3e | %s |' % (
            r['label'], r['n_freq'], r['signal']['rel']['max'],
            r['transition']['rel']['max'], r['DN_gw_rel'], 'ok'))
    L.append('')
    L.append('Aggregates over %d points: signal-band rel max **%.3e** (%s <1e-3), '
             'transition-band rel max **%.3e** (%s <1e-3), integrated Delta_Neff rel '
             'abs median **%.3e** / p95 **%.3e** / max **%.3e** (<1e-4: %s).'
             % (len(rows), sig['max'], 'PASS' if gate_sig else 'FAIL',
                tra['max'], 'PASS' if gate_tra else 'FAIL',
                dn['median'], dn['p95'], dn['max'],
                'PASS' if gate_dn else 'FAIL'))
    L.append('')
    L.append('Reference runtime per point (workers=%s) median %.0f s; fast plain-grid '
             'median %.2f s.' % (os.environ.get('PLAIN_REF_WORKERS', '3'),
                                 summary['runtime_s']['ref_median'],
                                 summary['runtime_s']['fast_median']))
    with open(os.path.join(args.out, 'validation_summary.md'), 'w',
              encoding='utf-8') as fh:
        fh.write('\n'.join(L) + '\n')
    print('summary written to %s' % os.path.join(args.out, 'validation_summary.json'))
    print('signal rel max %.3e | transition rel max %.3e | DN med %.3e' % (
        sig['max'], tra['max'], dn['median']))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase', choices=['plain', 'summary'], default='plain')
    ap.add_argument('--out', default=os.path.join('docs', 'paramsweep_plain'))
    ap.add_argument('--pool', type=int, default=3)
    ap.add_argument('--workers', type=int, default=3)
    args = ap.parse_args()
    if not os.path.isabs(args.out):
        args.out = os.path.join(ROOT, args.out)
    if args.phase == 'plain':
        os.environ['PLAIN_REF_WORKERS'] = str(args.workers)
        run_plain(args)
    else:
        return run_summary(args)
    return 0


if __name__ == '__main__':
    sys.exit(main())
