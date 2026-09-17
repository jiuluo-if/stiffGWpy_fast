"""Standalone A/B for reusing an unchanged gen_fast grid structure."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import statistics
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
}


def digest(value):
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def make_cached_gen_fast(original):
    def cached(m, h=0.01, kink_split=False):
        d = m.derived_param
        n_re_abs = d['N_inf'] - d['N_re']
        cache = getattr(m, '_gen_fast_grid_cache', None)
        can_reuse = (
            cache is not None and cache['h'] == h and
            cache['kink_split'] == kink_split and
            len(m.Nv) == cache['size'] and
            abs(n_re_abs - cache['n_re_abs']) < 0.25 * h and
            abs(m.Nv[-2] - cache['second_last']) < 1e-14
        )
        if not can_reuse:
            result = original(m, h=h, kink_split=kink_split)
            m._gen_fast_grid_cache = {
                'h': h,
                'kink_split': kink_split,
                'size': len(m.Nv),
                'index_re': int(np.argmin(np.abs(m.Nv - n_re_abs))),
                'n_re_abs': n_re_abs,
                'second_last': float(m.Nv[-2]),
            }
            return result

        # The only changed node in the audited outer updates is the continuous
        # present-day anchor. Keep all other grid coordinates and the index
        # contract bit-for-bit unchanged.
        m.Nv[-1] = d['N_inf']
        m.N = m.Nv - m.Nv[-1]
        index_re = cache['index_re']
        p = m.cosmo_param
        Omh2 = d['Omega_mh2']
        Osh2 = d['Omega_sh2']
        Oerh2 = gp.Omega_ph2 * 7 / 8 * (4 / 11) ** (4 / 3) * p['DN_eff']
        Otrh2 = gp.Omega_orh2 + Oerh2
        Otreh2 = gp.Omega_ph2 * gp.rho_th[-1] + Oerh2
        OLh2 = (d['h'] ** 2 - Omh2 - gp.Omega_mnuh2 - gp.Omega_nh2 * 2 / 3
                - gp.Omega_ph2 - Oerh2 - Osh2)
        Sv = np.empty(len(m.Nv))
        f_hor = np.empty(len(m.Nv))
        delta_f = math.log(2 * math.pi / d['H_0'])
        FS.gen_kernel(
            m.Nv, Sv, f_hor, index_re, Omh2, Osh2, Oerh2, Otrh2, Otreh2,
            OLh2, gp.Omega_mnuh2, gp.Omega_ph2, gp.Omega_nh2, gp.nu_today,
            gp.N_fin, gp.N_max, FS._FD_X0, FS._FD_DX, FS._FD_C_RHO,
            FS._FD_C_P, FS._TH_X, FS._TH_C_RHO, FS._TH_C_RHOP, delta_f,
            FS.ln10)
        m.Nv = m.Nv
        m.N = m.Nv - m.Nv[-1]
        m.sigma = Sv
        m.f_hor = f_hor
        m.f_re = f_hor[index_re]
        return cache['size'], index_re
    return cached


def run(case, repeats, cached):
    original = FS.gen_fast
    if cached:
        FS.gen_fast = make_cached_gen_fast(original)
    rows = []
    try:
        for _ in range(repeats):
            model = LCDM_SG(**CASES[case])
            start = time.perf_counter()
            result = FS.SGWB_iter_fast(model, kink_split=True)
            elapsed = time.perf_counter() - start
            rows.append({
                'total_s': elapsed,
                'converged': result is model and bool(model.SGWB_converge),
                'digest_f': digest(model.f),
                'digest_spectrum': digest(model.log10OmegaGW),
                'digest_DN_gw': digest(model.DN_gw),
                'DN_gw_last': float(np.asarray(model.DN_gw)[-1]),
            })
    finally:
        FS.gen_fast = original
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=sorted(CASES), default='default')
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--json', required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode('fast')
    baseline = run(args.case, args.repeats, cached=False)
    candidate = run(args.case, args.repeats, cached=True)
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'case_id': args.case,
        'threads': int(FS._THREADS),
        'repeats': args.repeats,
        'baseline_median_s': statistics.median(row['total_s'] for row in baseline),
        'candidate_median_s': statistics.median(row['total_s'] for row in candidate),
        'baseline_p95_s': float(np.percentile([row['total_s'] for row in baseline], 95)),
        'candidate_p95_s': float(np.percentile([row['total_s'] for row in candidate], 95)),
        'candidate_over_baseline': statistics.median(row['total_s'] for row in candidate) /
        statistics.median(row['total_s'] for row in baseline),
        'all_converged': all(row['converged'] for row in baseline + candidate),
        'digest_equal': all(
            left[key] == right[key]
            for left, right in zip(baseline, candidate)
            for key in ('digest_f', 'digest_spectrum', 'digest_DN_gw')
        ),
        'dn_max_abs': max(abs(left['DN_gw_last'] - right['DN_gw_last'])
                          for left, right in zip(baseline, candidate)),
    }
    with open(args.json, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
