"""Record per-frequency frozen-tail and adiabaticity diagnostics.

This is a science diagnostic, not a production solver or a tail correction.
The fast ``DN_eff`` and native frequency set are frozen while independent
reference modes are solved at several hand-off depths.
"""

import argparse
import hashlib
import json
import os
import sys

try:
    from scripts._resource_budget import apply_environment, nested_thread_budget, telemetry
except ImportError:
    from _resource_budget import apply_environment, nested_thread_budget, telemetry

apply_environment()

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benchmark_same_grid_reference import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

Z_TAILS = (5.0, 6.0, 7.0, 8.0, 10.0)
SCHEMA_VERSION = 1


def _commit():
    return os.popen('git rev-parse HEAD').read().strip()


def _reference_version():
    with open(REF.__file__, 'rb') as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def _cache_key(point, parameters, frequencies, dn_eff, z_tail, rtol):
    payload = {
        'commit': _commit(),
        'reference_version': _reference_version(),
        'point': point,
        'parameters': parameters,
        'frequencies': [float(value) for value in frequencies],
        'DN_eff': float(dn_eff),
        'z_tail': float(z_tail),
        'rtol': float(rtol),
        'solver_config': {'reference': 'DOP853', 'observable': 'mode tail diagnostics'},
    }
    encoded = json.dumps(payload, sort_keys=True,
                         separators=(',', ':')).encode('utf-8')
    return hashlib.sha256(encoded).hexdigest()


def _new_checkpoint():
    return {'schema_version': SCHEMA_VERSION, 'generated_commit': _commit(),
            'reference_version': _reference_version(), 'records': {}}


def _load_checkpoint(path):
    with open(path, encoding='utf-8') as handle:
        checkpoint = json.load(handle)
    expected = _new_checkpoint()
    if any(checkpoint.get(name) != expected[name]
           for name in ('schema_version', 'generated_commit', 'reference_version')):
        raise ValueError('tail diagnostic checkpoint does not match current code')
    return checkpoint


def _save_checkpoint(path, checkpoint):
    temporary = path + '.tmp'
    with open(temporary, 'w', encoding='utf-8') as handle:
        json.dump(checkpoint, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    os.replace(temporary, path)


def _frequencies(point, freq_count):
    FS.apply_accuracy_mode('fast')
    FS.set_z_tail(5.0)
    model = LCDM_SG(**CASES[point])
    FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                      frequency_quadrature='simpson')
    frequencies = np.sort(np.asarray(model.f, dtype=float))
    if freq_count is None:
        return model, frequencies
    if not 3 <= freq_count <= frequencies.size:
        raise ValueError('freq_count must be between 3 and the native grid size')
    picks = np.linspace(0, frequencies.size - 1, freq_count, dtype=int)
    return model, frequencies[picks]


def run_point(point, freq_count, rtol, workers, checkpoint, checkpoint_path,
              max_points=None):
    model, frequencies = _frequencies(point, freq_count)
    dn_eff = float(model.cosmo_param['DN_eff'])
    omega_nu = float(gp.Omega_nh2 / model.derived_param['h'] ** 2)
    rows = []
    records = checkpoint.setdefault('records', {})
    for freq in frequencies:
        for z_tail in Z_TAILS:
            key = _cache_key(point, CASES[point], [freq], dn_eff, z_tail, rtol)
            cached = records.get(key)
            if cached is not None:
                rows.append(dict(cached, cache_hit=True))
                continue
            solution = REF.solve_reference_mode(
                model, float(freq), dn_eff, z_tail=z_tail, rtol=rtol)
            row = {
                'point': point,
                'frequency': float(freq),
                'z_tail': float(z_tail),
                'Ogw_today': float(solution['Ogw_today']),
                'Oj_today': float(solution['Oj_today']),
                'Opgw_today': float(solution['Opgw_today']),
                'DN_weight': float(
                    gp.Neff0 * np.log(10.0)
                    * (solution['Ogw_today'] - solution['Oj_today'])
                    / omega_nu),
                'used_tail': bool(solution['used_tail']),
                'event_N': solution['event_N'],
                'phase_handoff': float(solution['phase_handoff']),
                'amplitude_handoff': float(solution['amplitude_handoff']),
                'omega_handoff': float(solution['omega_handoff']),
                'omega_prime_over_omega2': float(solution['omega_prime_over_omega2']),
                'omega_second_over_omega3': float(solution['omega_second_over_omega3']),
                'eps_handoff': solution['eps_handoff'],
                'matching_error_rel': solution['matching_error_rel'],
                'n_steps': int(solution['n_steps']),
                'cache_key': key,
                'cache_hit': False,
            }
            records[key] = row
            _save_checkpoint(checkpoint_path, checkpoint)
            rows.append(row)
            if max_points is not None and len(rows) >= max_points:
                return rows, frequencies, dn_eff
    return rows, frequencies, dn_eff


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', default=['default', 'lowT', 'highT', 'stiff'],
                        choices=sorted(CASES))
    parser.add_argument('--freq-count', type=int, default=8)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--max-points', type=int, default=None)
    parser.add_argument('--out', default='docs/tail_diagnostics_head.json')
    parser.add_argument('--checkpoint', default=None)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args(argv)
    inner_threads = nested_thread_budget(args.workers, args.threads)
    os.environ['NUMBA_NUM_THREADS'] = str(inner_threads)
    os.environ['FAST_THREADS'] = str(inner_threads)
    FS.set_threads(inner_threads)
    checkpoint_path = args.checkpoint or args.out + '.checkpoint.json'
    checkpoint = _load_checkpoint(checkpoint_path) if args.resume else _new_checkpoint()
    all_rows = []
    point_meta = []
    remaining = args.max_points
    for point in args.points:
        rows, frequencies, dn_eff = run_point(
            point, args.freq_count, args.rtol, args.workers, checkpoint,
            checkpoint_path, remaining)
        all_rows.extend(rows)
        point_meta.append({'point': point, 'frequencies': frequencies.tolist(),
                           'DN_eff': dn_eff})
        if remaining is not None:
            remaining -= len(rows)
            if remaining <= 0:
                break
    payload = {
        'schema_version': SCHEMA_VERSION,
        'generated_commit': _commit(),
        'reference_version': _reference_version(),
        'resources': telemetry(workers=args.workers, threads=inner_threads,
                               concurrent_processes=1),
        'rtol': args.rtol,
        'z_tails': list(Z_TAILS),
        'points': point_meta,
        'rows': all_rows,
        'checkpoint': checkpoint_path,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    print('wrote %s' % args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
