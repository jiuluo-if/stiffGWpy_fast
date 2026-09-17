"""Standalone 观测量敏感度代理诊断，不改变正式 outer reuse 判据。"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    'positive_tilt': dict(r=1e-2, cr=0, T_re=2e3, kappa10=1e-2, nt=0.2),
}


def _relative_max(a, b):
    return float(np.max(np.abs(a - b)) / max(1.0, float(np.max(np.abs(b)))))


def _run(name):
    saved = (FS.gen_fast, FS.prep_frequency_only, FS.solve_kernel,
             FS._OUTER_FULL_REUSE_SIGMA_TOL, FS._OUTER_FULL_REUSE_FHOR_TOL)
    sigma = []
    f_hor = []
    phases = []
    s2s = []
    spectra = []
    frequency_x = None

    def gen(*args, **kwargs):
        out = saved[0](*args, **kwargs)
        sigma.append(np.asarray(args[0].sigma, dtype=float).copy())
        return out

    def prep(*args, **kwargs):
        out = saved[1](*args, **kwargs)
        f_hor.append(np.asarray(out[1], dtype=float).copy())
        return out

    def solve(*args, **kwargs):
        nonlocal frequency_x
        out = saved[2](*args, **kwargs)
        ogw = np.asarray(args[16], dtype=float)
        oj = np.asarray(args[17], dtype=float)
        spectra.append(np.maximum(ogw[:, -1] - oj[:, -1], 1e-300).copy())
        frequency_x = np.log10(np.asarray(args[10], dtype=float).copy())
        phases.append(np.asarray(args[1], dtype=float).copy())
        s2s.append(np.asarray(args[3], dtype=float).copy())
        return out

    FS.gen_fast = gen
    FS.prep_frequency_only = prep
    FS.solve_kernel = solve
    try:
        model = LCDM_SG(**CASES[name])
        FS.apply_accuracy_mode('fast')
        FS.set_threads(2)
        # 禁用 shortcut 以捕获真实 first -> second outer update。
        FS._OUTER_FULL_REUSE_SIGMA_TOL = 0.0
        FS._OUTER_FULL_REUSE_FHOR_TOL = 0.0
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal')
    finally:
        (FS.gen_fast, FS.prep_frequency_only, FS.solve_kernel,
         FS._OUTER_FULL_REUSE_SIGMA_TOL, FS._OUTER_FULL_REUSE_FHOR_TOL) = saved

    if len(spectra) < 2:
        return {'point': name, 'kernel_calls': len(spectra), 'eligible': False}
    old_i = spectra[0]
    new_i = spectra[1]
    delta_i = np.log10(new_i) - np.log10(old_i)
    delta_phi = _relative_max(phases[0], phases[1]) if len(phases) >= 2 else None
    delta_s2 = _relative_max(s2s[0], s2s[1]) if len(s2s) >= 2 else None
    horizon_shift = (float(np.max(np.abs(f_hor[1] - f_hor[0])))
                     if len(f_hor) >= 2 else None)
    order = np.argsort(frequency_x)
    x_sorted = frequency_x[order]
    sensitivity = float(np.trapezoid(np.abs(new_i[order] - old_i[order]), x=x_sorted)
                        / max(np.trapezoid(new_i[order], x=x_sorted), 1e-300))
    return {
        'point': name,
        'kernel_calls': len(spectra),
        'eligible': True,
        'delta_log_omega_max': float(np.max(np.abs(delta_i))),
        'delta_phi_relative_max': delta_phi,
        'delta_s2_relative_max': delta_s2,
        'horizon_shift_max': horizon_shift,
        'frequency_weighted_dn_sensitivity_proxy': sensitivity,
        'dn_gw_first': float(model.DN_gw[0]),
        'dn_gw_last': float(np.asarray(model.DN_gw)[-1]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--points', nargs='+', choices=sorted(CASES), default=sorted(CASES))
    parser.add_argument('--out', default='docs/outer_observable_proxy_round.json')
    args = parser.parse_args()
    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': {'numba_threads': 2, 'blas_threads': 1, 'workers': 1},
        'forced_reuse': False,
        'rows': [_run(name) for name in args.points],
        'semantics': ('First-to-second outer update observable proxies on the same '
                      'forced-reuse solve; diagnostic only, not an error bound.'),
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
