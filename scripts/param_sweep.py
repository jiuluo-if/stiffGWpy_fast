# -*- coding: utf-8 -*-
"""
param_sweep.py -- audit phase 3: large parameter-space validation.

Deterministic sampling (Sobol + LHS + edge-oversampled LHS + hand-picked
extreme points, >= 1000 total); every point is solved by BOTH the LSODA
reference and the fast engine on the same sigma grid at default resolution
(h=0.01, z_tail=5, col_step=4, freq_res=1.0; LSODA rtol=1e-8, outer tol=1e-7).

Per point it records the convergence status, runtimes, final DN_eff / DN_gw /
kappa_r, the spectrum error (max dex, dex percentiles, max linear-Omega
error) and the DN_gw curve error.  Results are appended to
``docs/paramsweep/sweep_phase3.jsonl`` and the run is checkpointed: rerunning
skips point ids already present, so an interrupted sweep resumes in place.

Usage::

    python scripts/param_sweep.py --out docs/paramsweep --workers 4
    python scripts/param_sweep.py --out docs/paramsweep --limit 12   # smoke test
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from convergence_study import metrics, run_fast, run_lsoda  # noqa: E402
from scipy.stats import qmc  # noqa: E402

from stiffgwpy import fast_sgwb as FS  # noqa: E402
from stiffgwpy.stiff_SGWB import LCDM_SG  # noqa: E402

SEED = 20260830
N_SOBOL = 400
N_LHS = 350
N_EDGE = 200

# name, lo, hi, mode ('lin' | 'log' | 'disc')
PARAMS = [
    ('Omega_bh2', 0.018, 0.026, 'lin'),
    ('Omega_ch2', 0.09, 0.15, 'lin'),
    ('H0', 60.0, 76.0, 'lin'),
    ('DN_eff', 0.0, 2.0, 'lin'),
    ('A_s', 1e-9, 4e-9, 'log'),
    ('r', 1e-4, 1e-1, 'log'),
    ('n_t', -0.5, 0.5, 'lin'),
    ('cr', 0.0, 1.0, 'disc'),
    ('T_re', 10.0, 1e6, 'log'),
    ('DN_re', 0.0, 30.0, 'lin'),
    ('kappa10', 1e-4, 1.0, 'log'),
]
PARAM_NAMES = [p[0] for p in PARAMS]

BASE = dict(Omega_bh2=0.0223828, Omega_ch2=0.1201075, H0=67.32117,
            DN_eff=0.0, A_s=2.100549e-9, r=1e-2, n_t=0.0, cr=1.0,
            T_re=2e3, DN_re=0.0, kappa10=1e-2)

SWEEP_SETTINGS = dict(h=0.01, z_tail=5.0, col_step=4, freq_res=1.0,
                      lsoda_rtol=1e-8, outer_tol=1e-7)


def to_physical(u, lo, hi, mode):
    """Map a unit-square coordinate to a physical parameter value."""
    if mode == 'log':
        return lo * (hi / lo) ** u
    if mode == 'disc':
        return 1.0 if u >= 0.5 else 0.0
    return lo + (hi - lo) * u


def edge_push(u, strength=0.6):
    """Push unit-square coordinates toward the 0/1 boundaries (oversampling)."""
    return np.clip(0.5 + np.sign(u - 0.5) * np.abs(u - 0.5) ** strength, 0.0, 1.0)


def extreme_params():
    """Hand-picked extreme / boundary points (stiff-relevant corners + anchors)."""
    out = []
    for name, lo, hi, mode in PARAMS:
        for val in (lo, hi):
            p = dict(BASE)
            p[name] = float(val)
            if name == 'cr':
                p['cr'] = float(val)
            out.append(p)
    joints = [
        dict(r=1e-4, kappa10=1e-4, T_re=1e6, cr=0.0, n_t=0.0, DN_re=30.0),
        dict(r=1e-1, kappa10=1.0, T_re=10.0, cr=1.0),
        dict(r=1e-4, kappa10=1.0, T_re=10.0, cr=0.0, n_t=-0.5, DN_re=0.0),
        dict(r=1e-1, kappa10=1e-4, T_re=1e6, cr=1.0),
        dict(r=1e-2, kappa10=1.0, T_re=2e3, cr=1.0),
        dict(r=1e-2, kappa10=1e-4, T_re=2e3, cr=0.0, n_t=0.5, DN_re=30.0),
        dict(r=3.6e-2, kappa10=1.0, T_re=1e3, cr=0.0, n_t=0.0, DN_re=0.0),
        dict(r=1e-2, kappa10=1e-2, T_re=10.0, cr=1.0),
        dict(r=1e-2, kappa10=1e-2, T_re=2e3, cr=1.0),
        dict(r=1e-2, kappa10=1e-2, T_re=2e3, cr=0.0, n_t=0.0, DN_re=0.0),
        dict(r=1e-1, kappa10=1.0, T_re=1e6, cr=1.0),
        dict(r=1e-4, kappa10=1e-4, T_re=10.0, cr=0.0, n_t=0.0, DN_re=30.0),
        dict(r=1e-1, kappa10=1e-2, T_re=2e3, cr=0.0, n_t=-0.5, DN_re=30.0),
        dict(r=1e-4, kappa10=1e-2, T_re=2e3, cr=0.0, n_t=0.5, DN_re=0.0),
        dict(DN_eff=2.0, r=1e-1, cr=1.0, T_re=2e3, kappa10=1.0),
        dict(DN_eff=2.0, r=1e-4, cr=0.0, n_t=0.0, DN_re=30.0, T_re=10.0),
        dict(H0=60.0, r=1e-1, cr=1.0, T_re=1e6, kappa10=1e-4),
        dict(H0=76.0, r=1e-4, cr=0.0, n_t=-0.5, DN_re=0.0, T_re=1e6, kappa10=1.0),
        dict(Omega_bh2=0.018, Omega_ch2=0.15, H0=60.0, DN_eff=2.0, r=1e-1,
             T_re=10.0, kappa10=1.0, cr=1.0),
        dict(Omega_bh2=0.026, Omega_ch2=0.09, H0=76.0, r=1e-4, T_re=1e6,
             kappa10=1e-4, cr=1.0),
        dict(A_s=1e-9, r=1e-1, cr=1.0, T_re=2e3, kappa10=1.0),
        dict(A_s=4e-9, r=1e-4, cr=0.0, n_t=0.5, DN_re=30.0, T_re=10.0,
             kappa10=1e-4),
    ]
    for j in joints:
        p = dict(BASE)
        p.update(j)
        out.append(p)
    # cr=1 corner grid over the stiff-relevant pair (r, T_re)
    for r in (1e-4, 1e-3, 1e-2, 1e-1):
        for tre in (10.0, 1e2, 1e3, 1e4, 1e5, 1e6):
            p = dict(BASE, r=r, T_re=tre, cr=1.0)
            out.append(p)
    # cr=0 corner grid over (n_t, DN_re) at fixed stiff
    for nt in (-0.5, 0.0, 0.5):
        for dnre in (0.0, 5.0, 15.0, 30.0):
            p = dict(BASE, r=1e-2, cr=0.0, n_t=nt, DN_re=dnre, kappa10=1e-2)
            out.append(p)
    return out


def build_points():
    """Deterministic point list: Sobol + LHS + edge + extremes."""
    d = len(PARAMS)
    points = []
    sobol = qmc.Sobol(d=d, scramble=True, seed=SEED)
    for i, u in enumerate(sobol.random(N_SOBOL)):
        points.append({'id': 'sobol-%04d' % i, 'method': 'sobol',
                       'params': {name: float(to_physical(v, lo, hi, mode))
                                  for (name, lo, hi, mode), v in zip(PARAMS, u)}})
    lhs = qmc.LatinHypercube(d=d, seed=SEED + 1)
    for i, u in enumerate(lhs.random(N_LHS)):
        points.append({'id': 'lhs-%04d' % i, 'method': 'lhs',
                       'params': {name: float(to_physical(v, lo, hi, mode))
                                  for (name, lo, hi, mode), v in zip(PARAMS, u)}})
    edge = qmc.LatinHypercube(d=d, seed=SEED + 2)
    for i, u in enumerate(edge_push(edge.random(N_EDGE))):
        points.append({'id': 'edge-%04d' % i, 'method': 'edge',
                       'params': {name: float(to_physical(v, lo, hi, mode))
                                  for (name, lo, hi, mode), v in zip(PARAMS, u)}})
    for i, p in enumerate(extreme_params()):
        points.append({'id': 'extreme-%04d' % i, 'method': 'extreme',
                       'params': {name: float(p[name]) for name in PARAM_NAMES}})
    return points


def eval_point(pid, method, params):
    """Run LSODA then fast for one point and return the comparison record."""
    rec = {'id': pid, 'method': method, 'params': params,
           'settings': dict(SWEEP_SETTINGS)}
    kw = dict(params)
    lsoda = None
    fast = None
    t0 = time.perf_counter()
    try:
        ml = run_lsoda(kw, h=SWEEP_SETTINGS['h'], z_tail=SWEEP_SETTINGS['z_tail'],
                       freq_res=SWEEP_SETTINGS['freq_res'],
                       rtol=SWEEP_SETTINGS['lsoda_rtol'],
                       tol=SWEEP_SETTINGS['outer_tol'])
        lsoda = ml
    except Exception as exc:  # pragma: no cover - defensive
        rec['lsoda_error'] = '%s: %s' % (type(exc).__name__, exc)
    rec['t_s_lsoda'] = time.perf_counter() - t0
    t0 = time.perf_counter()
    try:
        mf = run_fast(kw, h=SWEEP_SETTINGS['h'], col_step=SWEEP_SETTINGS['col_step'],
                      z_tail=SWEEP_SETTINGS['z_tail'],
                      freq_res=SWEEP_SETTINGS['freq_res'],
                      tol=SWEEP_SETTINGS['outer_tol'])
        # run_fast always returns the model object; convergence is signalled
        # by SGWB_converge (the fast iteration may abort, e.g. N_eff > 5).
        fast = mf if getattr(mf, 'SGWB_converge', False) else None
    except Exception as exc:  # pragma: no cover - defensive
        rec['fast_error'] = '%s: %s' % (type(exc).__name__, exc)
    rec['t_s_fast'] = time.perf_counter() - t0
    if lsoda is None and fast is None:
        rec['status'] = 'both_failed'
        return rec
    if lsoda is None:
        rec['status'] = 'lsoda_failed'
        rec['DN_eff_fast'] = float(fast.cosmo_param['DN_eff'])
        rec['DN_gw_fast'] = float(fast.DN_gw[-1])
        rec['kappa_r_fast'] = float(fast.kappa_r)
        rec['iters_fast'] = int(getattr(fast, '_iters', -1))
        return rec
    if fast is None:
        rec['status'] = 'fast_failed'
        rec['DN_eff_lsoda'] = float(lsoda.cosmo_param['DN_eff'])
        rec['DN_gw_lsoda'] = float(lsoda.DN_gw[-1])
        rec['kappa_r_lsoda'] = float(lsoda.kappa_r)
        rec['iters_lsoda'] = int(getattr(lsoda, '_iters', -1))
        return rec
    met = metrics(lsoda, fast)
    rec['status'] = 'ok'
    rec['DN_eff_lsoda'] = met['DN_eff_ref']
    rec['DN_eff_fast'] = met['DN_eff_test']
    rec['DN_gw_lsoda'] = met['DN_gw_last_ref']
    rec['DN_gw_fast'] = met['DN_gw_last_test']
    rec['kappa_r_lsoda'] = float(lsoda.kappa_r)
    rec['kappa_r_fast'] = float(fast.kappa_r)
    rec['iters_lsoda'] = int(getattr(lsoda, '_iters', -1))
    rec['iters_fast'] = int(getattr(fast, '_iters', -1))
    rec['DN_gw_last_rel'] = met['DN_gw_last_rel']
    rec['DN_eff_rel'] = (abs(met['DN_eff_ref'] - met['DN_eff_test']) /
                         abs(met['DN_eff_ref']) if met['DN_eff_ref'] else None)
    rec['kappa_r_rel'] = met['kappa_r_rel']
    rec['dex_max'] = met['spectrum']['dex_max']
    rec['dex_p50'] = met['spectrum']['dex_p50']
    rec['dex_p95'] = met['spectrum']['dex_p95']
    rec['dex_p99'] = met['spectrum']['dex_p99']
    rec['lin_max'] = met['spectrum']['lin_max']
    rec['lin_median'] = met['spectrum']['lin_median']
    dg = met['DN_gw']
    rec['curve_dn_gw_max'] = dg['max'] if dg else None
    rec['curve_dn_gw_median'] = dg['median'] if dg else None
    return rec


def _worker(task):
    pid, method, params = task
    return eval_point(pid, method, params)


def main(argv=None):
    ap = argparse.ArgumentParser(description='phase-3 parameter-space sweep')
    ap.add_argument('--out', default='docs/paramsweep')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--limit', type=int, default=None,
                    help='max number of points to evaluate (smoke tests)')
    ap.add_argument('--offset', type=int, default=0)
    ap.add_argument('--no-warmup', action='store_true')
    ap.add_argument('--retry-failed', action='store_true',
                    help='re-evaluate points whose previous status was not ok')
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    points = build_points()
    with open(os.path.join(args.out, 'points.json'), 'w', encoding='utf-8') as fh:
        json.dump(points, fh, ensure_ascii=False, indent=1)
    out_path = os.path.join(args.out, 'sweep_phase3.jsonl')
    done = set()
    prev = {}
    if os.path.exists(out_path):
        with open(out_path, encoding='utf-8') as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                    done.add(r['id'])
                    prev[r['id']] = r
                except Exception:  # pragma: no cover - corrupt tail
                    pass
    if args.retry_failed:
        # Re-evaluate every point that previously failed (or never ran) so a
        # transient MemoryError / import hiccup is retried on the next call.
        todo = [p for p in points
                if p['id'] not in prev or prev[p['id']].get('status') != 'ok']
    else:
        todo = [p for p in points if p['id'] not in done]
    todo = todo[args.offset:]
    if args.limit is not None:
        todo = todo[:args.limit]
    if not todo:
        print('nothing to do (%d points already present)' % len(done))
        return 0
    if not args.no_warmup:
        # Compile all numba kernels once in the parent so pool workers only load cache.
        FS.set_col_step(SWEEP_SETTINGS['col_step'])
        FS.set_z_tail(SWEEP_SETTINGS['z_tail'])
        FS.set_h(SWEEP_SETTINGS['h'])
        mw = LCDM_SG(**BASE)
        FS.SGWB_iter_fast(mw, tol=SWEEP_SETTINGS['outer_tol'],
                          freq_res=SWEEP_SETTINGS['freq_res'])
        print('warmup done', flush=True)
    n = len(todo)
    t_start = time.perf_counter()
    with open(out_path, 'a', encoding='utf-8') as fh:
        # ProcessPoolExecutor workers are non-daemonic, so each can still spawn
        # the inner mp.Pool(4) used by run_SGWB (mp.Pool workers are daemonic
        # and would raise "daemonic processes are not allowed to have children").
        from concurrent.futures import ProcessPoolExecutor
        tasks = [(p['id'], p['method'], p['params']) for p in todo]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for i, rec in enumerate(pool.map(_worker, tasks, chunksize=1)):
                fh.write(json.dumps(rec, ensure_ascii=False) + '\n')
                fh.flush()
                if (i + 1) % 5 == 0 or i + 1 == n:
                    el = time.perf_counter() - t_start
                    rate = (i + 1) / el
                    eta = (n - i - 1) / rate if rate > 0 else float('nan')
                    print('[%d/%d] elapsed %.0fs eta %.0fs last=%s %s' %
                          (i + 1, n, el, eta, rec['id'], rec['status']),
                          flush=True)
    print('wrote %s (%d points)' % (out_path, n))
    return 0


if __name__ == '__main__':
    sys.exit(main())
