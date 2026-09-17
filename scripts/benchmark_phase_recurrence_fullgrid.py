"""Full native-grid oracle audit for the phase exponential recurrence prototype.

The formal solver is not changed.  For each named point this compares the
production kernel and the standalone recurrence twin on the complete native
goal grid, then compares both outputs with the independent Prüfer/DOP853
oracle at the same fixed ``DN_eff``.
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
    _make_args,
    _prepared,
    solve_kernel_recurrence,
)
from scripts.benchmark_prufer_oracle import (  # noqa: E402
    CASES as ORACLE_CASES,
)
from scripts.benchmark_prufer_oracle import (  # noqa: E402
    spectrum_prufer,
)
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402


def _stats(values):
    values = np.asarray(values, dtype=float)
    return {
        'p50': float(np.percentile(values, 50)),
        'p95': float(np.percentile(values, 95)),
        'max': float(np.max(values)),
        'n': int(values.size),
    }


def _integrated_dn(model, values):
    omega_nu = gp.Omega_nh2 / model.derived_param['h'] ** 2
    return float(gp.Neff0 * gp.ln10
                 * FS.integrate_frequency_pchip(model.f, values)
                 / omega_nu)


def run_case(point, threads, rtol):
    model, common = _prepared(point, threads, cases=ORACLE_CASES)
    baseline = _make_args(common)
    candidate = _make_args(common)
    FS.solve_kernel(*baseline)
    solve_kernel_recurrence(*candidate)

    base_omega = baseline[16][:, -1] - baseline[17][:, -1]
    cand_omega = candidate[16][:, -1] - candidate[17][:, -1]
    scale = max(float(np.max(np.abs(base_omega))), 1e-300)
    mask = np.abs(base_omega) > 1e-12 * scale
    recurrence_delta = np.abs(cand_omega[mask] - base_omega[mask]) \
        / np.maximum(np.abs(base_omega[mask]), 1e-300)

    dn_eff = float(model.cosmo_param['DN_eff'])
    oracle = spectrum_prufer(model, np.asarray(model.f, dtype=float), dn_eff,
                             z_tail=5.0, rtol=rtol)
    oracle_omega = oracle['Ogw'] - oracle['Oj']
    base_oracle = np.abs(base_omega - oracle_omega) \
        / np.maximum(np.abs(oracle_omega), 1e-300)
    cand_oracle = np.abs(cand_omega - oracle_omega) \
        / np.maximum(np.abs(oracle_omega), 1e-300)
    base_dn = _integrated_dn(model, base_omega)
    cand_dn = _integrated_dn(model, cand_omega)
    oracle_dn = float(oracle['DN_gw'])
    return {
        'point': point,
        'n_freq': int(len(model.f)),
        'dn_eff': dn_eff,
        'dn_gw': {
            'production': base_dn,
            'recurrence': cand_dn,
            'prufer': oracle_dn,
            'recurrence_vs_production_relative': abs(cand_dn - base_dn)
            / max(abs(base_dn), 1e-300),
            'production_vs_prufer_relative': abs(base_dn - oracle_dn)
            / max(abs(oracle_dn), 1e-300),
            'recurrence_vs_prufer_relative': abs(cand_dn - oracle_dn)
            / max(abs(oracle_dn), 1e-300),
        },
        'spectrum_relative': {
            'recurrence_vs_production': _stats(recurrence_delta),
            'production_vs_prufer': _stats(base_oracle),
            'recurrence_vs_prufer': _stats(cand_oracle),
        },
        'oracle': {
            'used_tail_fraction': float(oracle['used_tail_fraction']),
            'quadrature_error': float(oracle['quadrature_error']),
            'interpolation_error': oracle['interpolation_error'],
        },
        'status': 'ok',
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+',
                        choices=sorted(ORACLE_CASES),
                        default=['default', 'lowT', 'highT', 'stiff',
                                 'high_kappa'])
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--out', default='docs/phase_recurrence_fullgrid_20260917.json')
    args = parser.parse_args(argv)
    records = []
    for point in args.points:
        try:
            records.append(run_case(point, args.threads, args.rtol))
        except RuntimeError as exc:
            message = str(exc)
            status = ('physical_guard'
                      if 'shared_Neff_guard' in message else
                      'fast_preparation_failure')
            records.append({
                'point': point,
                'status': status,
                'message': message,
            })
    payload = {
        'schema_version': 1,
        'experiment': 'standalone_phase_exp_recurrence_full_native_grid_oracle',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resource': telemetry(workers=1, threads=args.threads,
                              concurrent_processes=1),
        'settings': {'points': args.points, 'threads': args.threads,
                     'rtol': args.rtol, 'freq_grid': 'goal',
                     'frequency_quadrature': 'pchip'},
        'records': records,
    }
    output = ROOT / args.out
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + '\n',
                      encoding='utf-8')
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
