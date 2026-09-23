"""Round 47 fresh canonical stage profile; diagnostic only."""
from __future__ import annotations

import json
import os
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import (  # noqa: E402
    apply_environment,
    limit_affinity,
    telemetry,
)

apply_environment()
import psutil  # noqa: E402

from scripts import profile_fast_breakdown as profile  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
}


def main() -> None:
    process = psutil.Process()
    selected = limit_affinity(process, 2)
    os.environ['FAST_THREADS'] = '2'
    original = profile.FS.SGWB_iter_fast

    def canonical(*args, **kwargs):
        kwargs['freq_grid'] = 'goal'
        kwargs['frequency_quadrature'] = 'pchip'
        return original(*args, **kwargs)

    profile.FS.SGWB_iter_fast = canonical
    rows = {}
    try:
        for name, case in CASES.items():
            records = []
            for _ in range(25):
                _, totals, model = profile.run_once(case, kink_split=True)
                row = {
                    key: float(values[0])
                    for key, values in totals.items()
                    if values and key not in ('outer_iteration', 'j0_cv', 'j0_span',
                                              'steps_per_channel')
                }
                row['outer_iterations'] = int(len(model.DN_gw))
                row['kernel_calls'] = int(len(totals['tensor_solve_kernel']))
                row['steps_per_channel'] = float(totals['steps_per_channel'][-1])
                row['j0_span'] = float(totals['j0_span'][-1])
                records.append(row)
            names = sorted(k for k in records[0]
                           if k not in ('outer_iterations', 'kernel_calls',
                                        'steps_per_channel', 'j0_span'))
            rows[name] = {
                'median_s': {f'{k}_s': statistics.median(r[k] for r in records)
                             for k in names},
                'p95_s': {f'{k}_s': float(__import__('numpy').percentile(
                    [r[k] for r in records], 95)) for k in names},
                'outer_iterations': sorted({r['outer_iterations'] for r in records}),
                'kernel_calls': sorted({r['kernel_calls'] for r in records}),
                'steps_per_channel': sorted({r['steps_per_channel'] for r in records}),
                'j0_span': sorted({r['j0_span'] for r in records}),
                'records': records,
            }
    finally:
        profile.FS.SGWB_iter_fast = original
    payload = {
        'experiment': 'round47_fresh_canonical_profile',
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(process, threads=2),
        'affinity_selected': selected,
        'threads': 2,
        'numba_threading_layer': os.environ.get('NUMBA_THREADING_LAYER'),
        'goal': True,
        'frequency_quadrature': 'pchip',
        'kink_split': True,
        'repeats': 25,
        'rows': rows,
    }
    out = Path('docs/profile_fast_breakdown_round47_20260923.json')
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
