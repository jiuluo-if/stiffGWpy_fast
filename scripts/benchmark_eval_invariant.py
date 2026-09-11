# -*- coding: utf-8 -*-
"""Check that native likelihood nodes do not alter the bolometric DN result."""
from __future__ import annotations

import json
import os
import sys

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:
    from _resource_budget import apply_environment, telemetry

apply_environment()
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'low_r': dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
}


def solve(kw, eval_freqs=None):
    model = LCDM_SG(**kw)
    config = FS.FastSolverConfig(h=0.005, col_step=8, z_tail=5.0,
                                 phase_max=0.25, freq_grid='goal', threads=1)
    FS.SGWB_iter_fast(model, kink_split=True, config=config,
                      eval_freqs=eval_freqs)
    return model


def main():
    FS.apply_accuracy_mode('fast')
    rows = []
    eval_freqs = np.array([-1.73, -0.42, 0.37, 1.21])
    for name, kw in CASES.items():
        base = solve(kw)
        with_eval = solve(kw, eval_freqs)
        delta = abs(float(with_eval.DN_gw[-1]) - float(base.DN_gw[-1]))
        rows.append({
            'case': name,
            'base_n_freq': int(base.f.size),
            'eval_n_freq': int(with_eval.f.size),
            'base_DN_gw': float(base.DN_gw[-1]),
            'eval_DN_gw': float(with_eval.DN_gw[-1]),
            'abs_delta': delta,
            'rel_delta': delta / abs(float(base.DN_gw[-1])),
            'failure_base': getattr(base, 'fast_failure_reason', None),
            'failure_eval': getattr(with_eval, 'fast_failure_reason', None),
        })
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(threads=1),
        'eval_freqs': eval_freqs.tolist(),
        'rows': rows,
    }
    with open('docs/benchmark_eval_invariant_head.json', 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
