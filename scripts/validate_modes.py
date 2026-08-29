# -*- coding: utf-8 -*-
"""
validate_modes.py -- audit: three recommended accuracy modes vs LSODA.

For each named accuracy mode in ``fast_sgwb.ACCURACY_MODES`` and a small set of
deterministic points, solves with the LSODA reference (same sigma grid, same
h/z_tail/freq_res, rtol=1e-8, outer tol=1e-7) and with the fast engine (mode
settings) and records the engine difference (Delta N_eff, DN_gw[-1], kappa_r,
spectrum dex, linear-Omega error), wall times and a rough RSS delta of the
fast call.  With ``--threads`` it also measures fast thread scaling for the
production mode.  Results are appended to ``docs/modes/validate_modes.jsonl``
and ``docs/modes/thread_scaling.json`` (checkpointed: existing mode/point pairs
are skipped unless ``--force``).

The ``reference`` preset runs LSODA at h=0.00125, which costs tens of minutes
per point; its same-grid engine difference is already measured in the phase-2
convergence study (docs/audit_phase2.md, h=0.00125 row), so by default the
reference mode is validated with the fast engine only (status ``fast_only``)
unless ``--include-reference-lsoda`` is given.

Usage::

    python scripts/validate_modes.py                 # default matrix
    python scripts/validate_modes.py --points default,p1
    python scripts/validate_modes.py --modes production --threads
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from convergence_study import metrics, run_fast, run_lsoda  # noqa: E402

from stiffgwpy import fast_sgwb as FS  # noqa: E402

try:
    import psutil
    _PROC = psutil.Process()
except Exception:  # pragma: no cover - optional measurement
    _PROC = None

POINTS = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'p1': dict(r=5e-3, cr=1, T_re=1e2, kappa10=1e-1),
    'p2': dict(r=1e-2, cr=0, n_t=0.2, DN_re=5.0, T_re=5e3, kappa10=1e-3),
}

LSODA_RTOL = 1e-8
LSODA_TOL = 1e-7


def _rss_mb():
    return _PROC.memory_info().rss / 1e6 if _PROC is not None else None


def one_case(mode_name, point_name, kw, include_lsoda=True):
    cfg = FS.ACCURACY_MODES[mode_name]
    rec = {'mode': mode_name, 'point': point_name, 'params': dict(kw),
           'cfg': dict(cfg)}
    ml = None
    if include_lsoda:
        t0 = time.perf_counter()
        try:
            ml = run_lsoda(kw, h=cfg['h'], z_tail=cfg['z_tail'],
                           freq_res=cfg['freq_res'], rtol=LSODA_RTOL,
                           tol=LSODA_TOL)
            rec['lsoda_status'] = 'ok'
        except Exception as exc:  # pragma: no cover - defensive
            rec['lsoda_status'] = 'failed'
            rec['lsoda_error'] = '%s: %s' % (type(exc).__name__, exc)
            ml = None
        rec['t_s_lsoda'] = time.perf_counter() - t0
    else:
        rec['lsoda_status'] = 'fast_only'
    rss0 = _rss_mb()
    t0 = time.perf_counter()
    try:
        mf = run_fast(kw, h=cfg['h'], col_step=cfg['col_step'],
                      z_tail=cfg['z_tail'], freq_res=cfg['freq_res'],
                      tol=cfg['tol'])
        mf = mf if getattr(mf, 'SGWB_converge', False) else None
        rec['fast_status'] = 'ok' if mf is not None else 'aborted'
    except Exception as exc:  # pragma: no cover - defensive
        rec['fast_status'] = 'failed'
        rec['fast_error'] = '%s: %s' % (type(exc).__name__, exc)
        mf = None
    rec['t_s_fast'] = time.perf_counter() - t0
    if _PROC is not None and rss0 is not None:
        rec['fast_rss_delta_mb'] = round(_rss_mb() - rss0, 1)
    if ml is not None and mf is not None:
        met = metrics(ml, mf)
        rec['DN_eff_rel'] = (abs(met['DN_eff_ref'] - met['DN_eff_test']) /
                             abs(met['DN_eff_ref']) if met['DN_eff_ref'] else None)
        rec['DN_gw_last_rel'] = met['DN_gw_last_rel']
        rec['kappa_r_rel'] = met['kappa_r_rel']
        rec['dex_max'] = met['spectrum']['dex_max']
        rec['dex_p50'] = met['spectrum']['dex_p50']
        rec['lin_max'] = met['spectrum']['lin_max']
        rec['lin_median'] = met['spectrum']['lin_median']
        if rec['t_s_fast']:
            rec['speedup'] = rec['t_s_lsoda'] / rec['t_s_fast']
    return rec


def thread_scaling(point_kw):
    cfg = FS.ACCURACY_MODES['production']
    run_fast(point_kw, h=cfg['h'], col_step=cfg['col_step'],
             z_tail=cfg['z_tail'], freq_res=cfg['freq_res'], tol=cfg['tol'])
    out = []
    for th in (1, 2, 4, 8, 16):
        if th > FS._MAX_THREADS:
            break
        FS.set_threads(th)
        times = []
        for _ in range(5):
            t0 = time.perf_counter()
            run_fast(point_kw, h=cfg['h'], col_step=cfg['col_step'],
                     z_tail=cfg['z_tail'], freq_res=cfg['freq_res'],
                     tol=cfg['tol'])
            times.append(time.perf_counter() - t0)
        out.append({'threads': th, 't_s_median': float(np.median(times))})
    FS.set_threads(min(FS.ACCURACY_MODES['production']['threads'],
                       FS._MAX_THREADS))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description='accuracy-mode validation vs LSODA')
    ap.add_argument('--out', default='docs/modes')
    ap.add_argument('--points', default=None,
                    help='comma list of %s' % ','.join(POINTS))
    ap.add_argument('--modes', default=None,
                    help='comma list of %s' % ','.join(FS.ACCURACY_MODES))
    ap.add_argument('--threads', action='store_true',
                    help='also measure fast thread scaling (production mode)')
    ap.add_argument('--include-reference-lsoda', action='store_true',
                    help='also run the expensive LSODA reference for the '
                         'reference mode (default: fast-only, cite phase 2)')
    ap.add_argument('--force', action='store_true',
                    help='re-run mode/point pairs already in the output file')
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    points = {k: v for k, v in POINTS.items()
              if args.points is None or k in args.points.split(',')}
    modes = [m for m in FS.ACCURACY_MODES
             if args.modes is None or m in args.modes.split(',')]
    out_path = os.path.join(args.out, 'validate_modes.jsonl')
    done = set()
    if os.path.exists(out_path):
        for line in open(out_path, encoding='utf-8'):
            try:
                r = json.loads(line)
                done.add((r['mode'], r['point']))
            except Exception:  # pragma: no cover - corrupt tail
                pass
    with open(out_path, 'a', encoding='utf-8') as fh:
        for pn, kw in points.items():
            for mn in modes:
                # The reference-mode LSODA at h=0.00125 is ~8x pricier; when
                # explicitly requested, run it only on the default point.
                if (mn == 'reference' and args.include_reference_lsoda
                        and pn != 'default'):
                    continue
                if (mn, pn) in done and not args.force:
                    continue
                include_lsoda = (mn != 'reference'
                                 or args.include_reference_lsoda)
                rec = one_case(mn, pn, kw, include_lsoda=include_lsoda)
                fh.write(json.dumps(rec, ensure_ascii=False) + '\n')
                fh.flush()
                print('[%s/%s] mode=%-11s point=%-8s lsoda=%-6s fast=%-7s '
                      't_lsoda=%6.1fs t_fast=%7.4fs speedup=%8.0fx' %
                      (len(done) + 1, len(points) * len(modes), mn, pn,
                       rec.get('lsoda_status', '?'),
                       rec.get('fast_status', '?'),
                       rec.get('t_s_lsoda', float('nan')),
                       rec.get('t_s_fast', float('nan')),
                       rec.get('speedup', float('nan'))), flush=True)
                done.add((mn, pn))
    if args.threads:
        ts = thread_scaling(POINTS['default'])
        with open(os.path.join(args.out, 'thread_scaling.json'), 'w',
                  encoding='utf-8') as fh:
            json.dump({'point': 'default', 'mode': 'production',
                       'threads': ts}, fh, ensure_ascii=False, indent=1)
        for row in ts:
            print('threads=%2d t_fast_median=%.4fs' %
                  (row['threads'], row['t_s_median']), flush=True)
    print('wrote %s' % out_path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
