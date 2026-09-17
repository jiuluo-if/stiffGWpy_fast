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
    snapshots = []
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
        snapshot = tuple(value.copy() if isinstance(value, np.ndarray) else value
                         for value in args)
        snapshots.append(snapshot)
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

    def run_input_variant(kind):
        variant = list(snapshots[0])
        if kind == 'phi':
            variant[1] = variant[1] + (snapshots[1][1] - snapshots[0][1])
            variant[2] = variant[2] + (snapshots[1][2] - snapshots[0][2])
        elif kind == 's2':
            variant[3] = variant[3] + (snapshots[1][3] - snapshots[0][3])
            variant[4] = variant[4] + (snapshots[1][4] - snapshots[0][4])
        else:
            raise ValueError(kind)
        variant[11] = 1
        for index in (16, 17, 18):
            variant[index] = np.zeros_like(snapshots[0][index])
        variant[22] = np.full_like(snapshots[0][22], -1.0)
        saved[2](*tuple(variant))
        return np.maximum(variant[16][:, -1] - variant[17][:, -1], 1e-300)

    phi_response = run_input_variant('phi')
    s2_response = run_input_variant('s2')
    phi_pred = float(np.max(np.abs(np.log10(phi_response) - np.log10(old_i))))
    s2_pred = float(np.max(np.abs(np.log10(s2_response) - np.log10(old_i))))
    return {
        'point': name,
        'kernel_calls': len(spectra),
        'outer_kernel_calls': len(spectra),
        'response_kernel_calls': 2,
        'eligible': True,
        'delta_log_omega_max': float(np.max(np.abs(delta_i))),
        'delta_phi_relative_max': delta_phi,
        'delta_s2_relative_max': delta_s2,
        'horizon_shift_max': horizon_shift,
        'frequency_weighted_dn_sensitivity_proxy': sensitivity,
        'linear_response_predicted_dlogomega_phi': phi_pred,
        'linear_response_predicted_dlogomega_s2': s2_pred,
        'linear_response_actual_dlogomega': float(np.max(np.abs(delta_i))),
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
