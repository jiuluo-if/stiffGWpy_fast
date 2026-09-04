# -*- coding: utf-8 -*-
"""
validate_fast.py -- precision gate for the fast solver (12 fixed cases).

Compares the original LSODA ``SGWB_iter()`` with ``fast_sgwb.SGWB_iter_fast()``
on a fixed parameter grid and exits non-zero when a tolerance gate is
exceeded, so the script can be used as a CI-style accuracy gate.

Error metrics (per the independent audit):
  * Omega_GW: absolute dex difference and relative difference of the *linear*
    Omega value, evaluated on the signal region (Omega >= 1e-30);
  * DN_gw[-1] / kappa_r: plain relative difference;
  * full-curve DN_gw / g2 / hubble: relative difference plus length match.

Usage:
    python scripts/validate_fast.py [--cases 0 3 11] [--out results.jsonl]
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
from stiffgwpy._metrics import dex_abs, rel_abs, rel_linear_omega, signal_mask
from stiffgwpy.stiff_SGWB import LCDM_SG

CASES = [
    dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    dict(r=3.6e-2, cr=1, T_re=2e3, kappa10=1e-2),
    dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-3),
    dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    dict(r=3.6e-2, cr=0, T_re=1e3, kappa10=1.0),
    dict(r=1e-2, cr=0, T_re=2e3, kappa10=1e-2),
    dict(r=1e-1, cr=1, T_re=1e1, kappa10=1.0),
    dict(r=1e-2, cr=0, T_re=1e4, kappa10=1e-1),
]

# Tolerance gates.  Measured values on the 12-case grid are roughly one order
# of magnitude below these limits (see docs/audit_final_report.md).
GATES = {
    'status_match': True,
    'DN_gw_last_rel': 1e-3,
    'kappa_r_rel': 1e-3,
    'hubble_maxrel': 1e-2,
    'log10OmegaGW_maxabs': 1e-2,
    'OmegaGW_rel_linear_max': 1e-2,
}


def env_meta():
    meta = {
        'python': sys.version.split()[0],
        'platform': sys.platform,
        'cpu_count': os.cpu_count(),
        'threads': FS._THREADS,
        'col_step': FS._COL_STEP,
        'numpy': np.__version__,
    }
    for mod, name in (('numba', 'numba'), ('scipy', 'scipy')):
        try:
            meta[name] = __import__(mod).__version__
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


def run_case(idx, kw):
    rec = dict(idx=idx, kw={k: repr(v) for k, v in kw.items()})
    m = LCDM_SG(**kw)
    if m.derived_param['N_inf'] is None:
        rec['skip'] = 'N_inf None (invalid combo)'
        return rec
    mo = LCDM_SG(**kw)
    t0 = time.perf_counter()
    mo.SGWB_iter(engine='lsoda')
    to = time.perf_counter() - t0
    rec['t_orig'] = to
    mf = LCDM_SG(**kw)
    t0 = time.perf_counter()
    FS.SGWB_iter_fast(mf)
    tf1 = time.perf_counter() - t0
    mf2 = LCDM_SG(**kw)
    t0 = time.perf_counter()
    FS.SGWB_iter_fast(mf2)
    tf2 = time.perf_counter() - t0
    rec['t_fast_first_ms'] = tf1 * 1e3
    rec['t_fast_warm_ms'] = tf2 * 1e3
    rec['speedup'] = to / tf2
    rec['conv_orig'] = bool(mo.SGWB_converge)
    rec['conv_fast'] = bool(mf2.SGWB_converge)
    rec['status_match'] = bool(mo.SGWB_converge == mf2.SGWB_converge)
    if not (mo.SGWB_converge and mf2.SGWB_converge):
        rec['note'] = 'not both converged'
        return rec
    for name in ('f', 'hubble'):
        a = np.asarray(getattr(mo, name), float)
        b = np.asarray(getattr(mf2, name), float)
        rec['len_' + name] = (int(a.size), int(b.size))
        if a.size == b.size and b.size:
            rec[name + '_maxrel'] = float(rel_abs(a, b).max())
        else:
            rec[name + '_maxrel'] = None
    lo_old = np.asarray(mo.log10OmegaGW, float)
    lo_new = np.asarray(mf2.log10OmegaGW, float)
    rec['len_log10OmegaGW'] = (int(lo_old.size), int(lo_new.size))
    if lo_old.size == lo_new.size:
        rec['log10OmegaGW_maxabs'] = float(dex_abs(lo_old, lo_new).max())
        mask = signal_mask(lo_new)
        rec['n_signal'] = int(mask.sum())
        rec['OmegaGW_rel_linear_max'] = (
            float(rel_linear_omega(lo_old[mask], lo_new[mask]).max()) if mask.any() else 0.0)
    else:
        rec['log10OmegaGW_maxabs'] = None
        rec['OmegaGW_rel_linear_max'] = None
    a = np.asarray(mo.DN_gw, float)
    b = np.asarray(mf2.DN_gw, float)
    rec['len_DN_gw'] = (int(a.size), int(b.size))
    fin = np.isfinite(a) & np.isfinite(b) & (np.abs(b) > 0)
    rec['DN_gw_fin_maxrel'] = float(rel_abs(a, b)[fin].max()) if fin.any() else None
    rec['DN_gw_last'] = (float(a[-1]), float(b[-1]))
    rec['DN_gw_last_rel'] = (float(abs(a[-1] - b[-1]) / abs(a[-1]))
                             if a[-1] != 0 else None)
    rec['kappa_r'] = (float(mo.kappa_r), float(mf2.kappa_r))
    rec['kappa_r_rel'] = float(abs(mo.kappa_r - mf2.kappa_r) / abs(mo.kappa_r))
    rec['DN_eff_final'] = (float(mo.cosmo_param['DN_eff']),
                           float(mf2.cosmo_param['DN_eff']))
    g = rel_abs(mo.g2, mf2.g2)
    g2f = np.isfinite(mo.g2) & np.isfinite(mf2.g2) & (np.abs(mf2.g2) > 0)
    rec['g2_fin_maxrel'] = float(g[g2f].max()) if g2f.any() else None
    return rec


def main(argv=None):
    ap = argparse.ArgumentParser(description='precision gate for the stiffgwpy fast solver')
    ap.add_argument('--cases', nargs='+', type=int, default=None,
                    help='case indices (default: all 12)')
    ap.add_argument('--out', default=None,
                    help='JSONL output path (default: validate_results.jsonl)')
    args = ap.parse_args(argv)
    which = args.cases if args.cases is not None else list(range(len(CASES)))
    out_path = args.out or os.path.join(os.getcwd(), 'validate_results.jsonl')
    results = []
    failures = []
    for i in which:
        t0 = time.perf_counter()
        r = run_case(i, CASES[i])
        r['wall_s'] = time.perf_counter() - t0
        results.append(r)
        print(json.dumps(r, ensure_ascii=False), flush=True)
        if r.get('skip'):
            continue
        if not r.get('status_match'):
            failures.append('case %d: status_match' % i)
        for key, limit in GATES.items():
            if key == 'status_match':
                continue
            val = r.get(key)
            if val is None:
                failures.append('case %d: %s missing' % (i, key))
            elif val > limit:
                failures.append('case %d: %s = %.3e > %.3e' % (i, key, val, limit))
    payload = {'meta': env_meta(), 'gates': GATES, 'results': results}
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(json.dumps(payload, ensure_ascii=False, indent=2) + '\n')
    print('wrote %s' % out_path)
    if failures:
        print('FAILED gates:')
        for f in failures:
            print('  - ' + f)
        return 1
    print('all gates passed (%d case(s))' % len(results))
    return 0


if __name__ == '__main__':
    sys.exit(main())
