# -*- coding: utf-8 -*-
"""
validate_random.py -- randomized parameter-space precision gate (P1).

Samples deterministic stratified-random points over the full 11-parameter
model space (beyond the 12 hand-picked cases of ``validate_fast.py``) and
compares the original LSODA solver with the fast solver.  Exits non-zero when
any tolerance gate is exceeded.  This is the starting point of the
full-parameter-space certification requested by the audit; run with a larger
``--n`` for a real validation campaign.

Usage:
    python scripts/validate_random.py [--n 10] [--seed 20260829] [--out r.jsonl]
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stiffgwpy import fast_sgwb as FS
from stiffgwpy._metrics import dex_abs, rel_linear_omega, signal_mask
from stiffgwpy.stiff_SGWB import LCDM_SG

GATES = {
    'status_match': True,
    'DN_gw_last_rel': 1e-3,
    'kappa_r_rel': 1e-3,
    'log10OmegaGW_maxabs': 2e-2,
    'OmegaGW_rel_linear_max': 2e-2,
}

# (name, low, high, log_scale)
_RANGES = {
    'Omega_bh2': (0.018, 0.026, False),
    'Omega_ch2': (0.09, 0.15, False),
    'H0': (60.0, 74.0, False),
    'DN_eff': (0.0, 0.5, False),
    'A_s': (1e-9, 3e-9, True),
    'r': (1e-4, 1e-1, True),
    'cr': (0.0, 1.0, False),
    'T_re': (1e1, 1e4, True),
    'kappa10': (1e-3, 1.0, True),
    'n_t': (-0.1, 0.1, False),
    'DN_re': (0.0, 5.0, False),
}


def sample_points(n, seed):
    rng = np.random.default_rng(seed)
    points = []
    for j in range(n):
        kw = {}
        for name, (lo, hi, log_scale) in _RANGES.items():
            u = (j + rng.uniform()) / n       # stratified random in [0, 1)
            u = min(max(u, 1e-6), 1.0 - 1e-6)
            val = 10 ** (np.log10(lo) + u * (np.log10(hi) - np.log10(lo))) if log_scale \
                else lo + u * (hi - lo)
            if name == 'cr':
                kw[name] = 1 if val >= 0.5 else 0
            else:
                kw[name] = float(val)
        if kw['cr'] == 0:
            kw['n_t'] = float(kw['n_t'])
            kw['DN_re'] = float(kw['DN_re'])
        else:
            kw.pop('n_t')
            kw.pop('DN_re')
        points.append(kw)
    return points


def env_meta():
    meta = {
        'python': sys.version.split()[0],
        'platform': sys.platform,
        'cpu_count': os.cpu_count(),
        'threads': FS._THREADS,
        'col_step': FS._COL_STEP,
        'numpy': np.__version__,
    }
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


def run_point(idx, kw):
    rec = dict(idx=idx, kw={k: repr(v) for k, v in kw.items()})
    m = LCDM_SG(**kw)
    if m.derived_param['N_inf'] is None:
        rec['skip'] = 'N_inf None (invalid combo)'
        return rec
    mo = LCDM_SG(**kw)
    t0 = time.perf_counter()
    mo.SGWB_iter(engine='lsoda')
    rec['t_orig'] = time.perf_counter() - t0
    mf = LCDM_SG(**kw)
    FS.SGWB_iter_fast(mf)
    mf2 = LCDM_SG(**kw)
    t0 = time.perf_counter()
    FS.SGWB_iter_fast(mf2)
    rec['t_fast_warm_ms'] = (time.perf_counter() - t0) * 1e3
    rec['conv_orig'] = bool(mo.SGWB_converge)
    rec['conv_fast'] = bool(mf2.SGWB_converge)
    rec['status_match'] = bool(mo.SGWB_converge == mf2.SGWB_converge)
    if not (mo.SGWB_converge and mf2.SGWB_converge):
        rec['note'] = 'not both converged'
        return rec
    lo_old = np.asarray(mo.log10OmegaGW, float)
    lo_new = np.asarray(mf2.log10OmegaGW, float)
    if lo_old.size == lo_new.size:
        rec['log10OmegaGW_maxabs'] = float(dex_abs(lo_old, lo_new).max())
        mask = signal_mask(lo_new)
        rec['OmegaGW_rel_linear_max'] = (
            float(rel_linear_omega(lo_old[mask], lo_new[mask]).max()) if mask.any() else 0.0)
    else:
        rec['log10OmegaGW_maxabs'] = None
        rec['OmegaGW_rel_linear_max'] = None
    a = np.asarray(mo.DN_gw, float)
    b = np.asarray(mf2.DN_gw, float)
    rec['DN_gw_last_rel'] = (float(abs(a[-1] - b[-1]) / abs(a[-1]))
                             if a[-1] != 0 else None)
    rec['kappa_r_rel'] = float(abs(mo.kappa_r - mf2.kappa_r) / abs(mo.kappa_r))
    rec['DN_eff_final'] = (float(mo.cosmo_param['DN_eff']),
                           float(mf2.cosmo_param['DN_eff']))
    return rec


def main(argv=None):
    ap = argparse.ArgumentParser(description='randomized precision gate')
    ap.add_argument('--n', type=int, default=10, help='number of random points')
    ap.add_argument('--seed', type=int, default=20260829)
    ap.add_argument('--out', default=None, help='JSON output path')
    args = ap.parse_args(argv)
    points = sample_points(args.n, args.seed)
    out_path = args.out or os.path.join(os.getcwd(), 'validate_random_results.jsonl')
    results = []
    failures = []
    for i, kw in enumerate(points):
        t0 = time.perf_counter()
        r = run_point(i, kw)
        r['wall_s'] = time.perf_counter() - t0
        results.append(r)
        print(json.dumps(r, ensure_ascii=False), flush=True)
        if r.get('skip'):
            continue
        if not r.get('status_match'):
            failures.append('point %d: status_match' % i)
        for key, limit in GATES.items():
            if key == 'status_match':
                continue
            val = r.get(key)
            if val is None:
                failures.append('point %d: %s missing' % (i, key))
            elif val > limit:
                failures.append('point %d: %s = %.3e > %.3e' % (i, key, val, limit))
    payload = {'meta': env_meta(), 'n': args.n, 'seed': args.seed,
               'gates': GATES, 'results': results}
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(json.dumps(payload, ensure_ascii=False, indent=2) + '\n')
    print('wrote %s' % out_path)
    if failures:
        print('FAILED gates:')
        for f in failures:
            print('  - ' + f)
        return 1
    n_skip = sum(1 for r in results if r.get('skip'))
    print('all gates passed (%d point(s), %d skipped as invalid combos)' %
          (len(results), n_skip))
    return 0


if __name__ == '__main__':
    sys.exit(main())
