"""Full outer self-consistency matrix for the phase recurrence prototype.

This is a reference-only candidate audit.  It monkeypatches the in-memory
``solve_kernel`` binding for one model at a time and restores it immediately;
the formal fast source and public API are not changed.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import (  # noqa: E402
    solve_kernel_recurrence,
)
from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _run(point, candidate, threads):
    FS.apply_accuracy_mode('fast')
    FS.set_threads(threads)
    model = LCDM_SG(**CASES[point])
    original = FS.solve_kernel
    try:
        if candidate:
            FS.solve_kernel = solve_kernel_recurrence
        result = FS.SGWB_iter_fast(
            model, kink_split=True, freq_grid='goal',
            frequency_quadrature='pchip')
    finally:
        FS.solve_kernel = original
    if result is None:
        reason = getattr(model, 'fast_failure_reason', 'unknown')
        status = 'physical_guard' if reason == 'shared_Neff_guard' else 'failure'
        return {'status': status, 'failure': reason}
    return {
        'status': 'ok',
        'n_freq': int(len(model.f)),
        'spectrum': np.asarray(model.log10OmegaGW, dtype=float),
        'DN_gw': float(np.asarray(model.DN_gw, dtype=float)[-1]),
    }


def _stats(values):
    values = np.asarray(values, dtype=float)
    return {
        'p50': float(np.percentile(values, 50)),
        'p95': float(np.percentile(values, 95)),
        'max': float(np.max(values)),
        'n': int(values.size),
    }


def run_case(point, threads):
    base = _run(point, False, threads)
    candidate = _run(point, True, threads)
    status_match = base['status'] == candidate['status']
    row = {
        'point': point,
        'base_status': base['status'],
        'candidate_status': candidate['status'],
        'status_match': status_match,
        'base_failure': base.get('failure'),
        'candidate_failure': candidate.get('failure'),
    }
    if base['status'] == 'ok' and candidate['status'] == 'ok':
        if base['n_freq'] != candidate['n_freq']:
            row.update({
                'status_match': False,
                'failure': 'frequency_count_mismatch',
                'base_n_freq': base['n_freq'],
                'candidate_n_freq': candidate['n_freq'],
            })
            return row
        delta = np.abs(candidate['spectrum'] - base['spectrum'])
        row.update({
            'n_freq': base['n_freq'],
            'spectrum_delta_dex': _stats(delta),
            'DN_rel': abs(candidate['DN_gw'] - base['DN_gw'])
            / max(abs(base['DN_gw']), 1e-300),
            'base_DN_gw': base['DN_gw'],
            'candidate_DN_gw': candidate['DN_gw'],
        })
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', choices=sorted(CASES),
                        default=sorted(CASES))
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--out', default='docs/phase_recurrence_outer_matrix_20260917.json')
    args = parser.parse_args(argv)
    rows = [run_case(point, args.threads) for point in args.points]
    false_safe = sum(
        row['base_status'] == 'physical_guard'
        and row['candidate_status'] == 'ok' for row in rows)
    accepted = [row for row in rows if row.get('spectrum_delta_dex')]
    payload = {
        'schema_version': 1,
        'experiment': 'standalone_phase_exp_recurrence_outer_matrix',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resource': telemetry(workers=1, threads=args.threads,
                              concurrent_processes=1),
        'settings': {'points': args.points, 'threads': args.threads,
                     'frequency_quadrature': 'pchip', 'freq_grid': 'goal',
                     'kink_split': True},
        'summary': {
            'total': len(rows),
            'base_ok': sum(row['base_status'] == 'ok' for row in rows),
            'candidate_ok': sum(row['candidate_status'] == 'ok' for row in rows),
            'status_mismatch': sum(not row['status_match'] for row in rows),
            'false_safe_count': false_safe,
            'max_spectrum_delta_dex': max(
                (row['spectrum_delta_dex']['max'] for row in accepted),
                default=None),
            'max_DN_rel': max((row['DN_rel'] for row in accepted), default=None),
        },
        'rows': rows,
    }
    output = ROOT / args.out
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + '\n',
                      encoding='utf-8')
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
