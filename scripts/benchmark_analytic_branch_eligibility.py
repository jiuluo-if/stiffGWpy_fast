# -*- coding: utf-8 -*-
"""统计实际 kernel 中可进入常系数 radiation/stiff branch 的 segment 比例。"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import exact_background as EB  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
}


def run_case(name):
    model = LCDM_SG(**CASES[name])
    FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal')
    nv = model.Nv.astype(float)
    phi, _, _, _, *_ = EB.fast_phi_s2_split(
        model, nv, model.cosmo_param['DN_eff'], sigma_nodes=model.sigma,
    )
    _, _, j0s, z0s, _ = FS.prep_frequency_only(model, nv, model.f)
    total = radiation_exact = stiff_exact = eligible = 0
    for j0, z0 in zip(j0s, z0s):
        k = int(j0)
        z = float(z0)
        while k < len(nv) - 1 and z < FS._Z_TAIL:
            z_end = float(z0 + phi[k + 1] - phi[k])
            sigma_left = float(model.sigma[k])
            sigma_right = float(model.sigma[k + 1])
            total += 1
            if abs(sigma_left - 4 / 3) < 1e-6 and abs(sigma_right - 4 / 3) < 1e-6:
                radiation_exact += 1
            if abs(sigma_left - 2) < 1e-6 and abs(sigma_right - 2) < 1e-6:
                stiff_exact += 1
            if (abs(sigma_left - 4 / 3) < 1e-4 and abs(sigma_right - 4 / 3) < 1e-4) \
                    or (abs(sigma_left - 2) < 1e-4 and abs(sigma_right - 2) < 1e-4):
                eligible += 1
            k += 1
            z = z_end
    return {
        'case': name,
        'n_modes': int(len(j0s)),
        'segments_before_tail': total,
        'radiation_exact_segments': radiation_exact,
        'stiff_exact_segments': stiff_exact,
        'eligible_segments_sigma_1e-4': eligible,
        'eligible_fraction': eligible / max(total, 1),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode('fast')
    FS.set_threads(2)
    payload = {
        'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'experiment': 'analytic_branch_eligibility_standalone',
        'resources': {'numba_threads': 2, 'blas_threads': 1, 'workers': 1},
        'z_tail': FS._Z_TAIL,
        'rows': [run_case(name) for name in CASES],
    }
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
