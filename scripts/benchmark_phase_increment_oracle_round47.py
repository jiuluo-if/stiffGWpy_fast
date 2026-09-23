"""Independent continuous-sigma oracle screen for the Round 47 twin."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts._resource_budget import apply_environment  # noqa: E402

apply_environment()

from scripts.benchmark_phase_increment_kernel_twin import (  # noqa: E402
    solve_kernel_phase_increment,
)
from scripts.benchmark_same_grid_reference import CASES  # noqa: E402
from scripts.validate_fast_vs_reference import reference_point  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


def main():
    FS.set_threads(2)
    original = FS.solve_kernel

    def candidate(*args, **kwargs):
        # validate_fast_vs_reference uses the non-kink production call shape.
        if len(args) == 23:
            args = args + (-1, 0.0, 0.0)
        return solve_kernel_phase_increment(*args, 32)

    FS.solve_kernel = candidate
    try:
        records = []
        for name in ('default', 'highT', 'sobol_010'):
            rec = reference_point(
                dict(CASES[name]), z_tail=5.0, rtol=1e-9,
                freq_res=1.0, label=name, kind='round47_candidate')
            records.append(rec)
    finally:
        FS.solve_kernel = original
    payload = {
        'experiment': 'round47_phase_increment_independent_oracle',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'settings': {'threads': 2, 'z_tail': 5.0, 'rtol': 1e-9,
                     'freq_res': 1.0, 'reanchor': 32},
        'records': records,
    }
    out = ROOT / 'docs/phase_increment_oracle_round47_20260923.json'
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + '\n',
                   encoding='utf-8')
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
