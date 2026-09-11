"""Measure independent-reference sensitivity to the frozen-tail hand-off.

The fast spectrum and its ``DN_eff`` are held fixed.  Only the continuous-
sigma DOP853 reference ``z_tail`` changes, so the result is an oracle-floor
diagnostic rather than a fast-solver benchmark.
"""

import argparse
import hashlib
import json
import os
import sys
import time

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benchmark_same_grid_reference import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

Z_TAILS = (5.0, 6.0, 7.0, 8.0, 10.0)
SCHEMA_VERSION = 2


def _commit() -> str:
    return os.popen('git rev-parse HEAD').read().strip()


def _reference_version() -> str:
    with open(REF.__file__, 'rb') as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def _cache_key(point, parameters, frequencies, dn_eff, z_tail, rtol):
    payload = {
        'commit': _commit(),
        'parameters': parameters,
        'frequencies': [float(value) for value in frequencies],
        'DN_eff': float(dn_eff),
        'z_tail': float(z_tail),
        'rtol': float(rtol),
        'reference_version': _reference_version(),
        'point': point,
    }
    encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True,
                         separators=(',', ':')).encode('utf-8')
    return hashlib.sha256(encoded).hexdigest()


def _new_checkpoint():
    return {
        'schema_version': SCHEMA_VERSION,
        'generated_commit': _commit(),
        'reference_version': _reference_version(),
        'records': {},
    }


def _load_checkpoint(path):
    with open(path, encoding='utf-8') as handle:
        checkpoint = json.load(handle)
    expected = _new_checkpoint()
    if (checkpoint.get('schema_version') != expected['schema_version'] or
            checkpoint.get('generated_commit') != expected['generated_commit'] or
            checkpoint.get('reference_version') != expected['reference_version']):
        raise ValueError('checkpoint 与当前 commit/schema/reference version 不一致')
    checkpoint.setdefault('records', {})
    return checkpoint


def _save_checkpoint(path, checkpoint):
    temporary = path + '.tmp'
    with open(temporary, 'w', encoding='utf-8') as handle:
        json.dump(checkpoint, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    os.replace(temporary, path)


def run_point(point, rtol=1e-10, workers=1, freq_count=None,
              freq_min=None, freq_max=None, checkpoint=None,
              checkpoint_path=None):
    """Run one fixed-background tail sweep and return an auditable record."""
    kw = CASES[point]
    FS.apply_accuracy_mode('fast')
    FS.set_z_tail(5.0)
    fast = LCDM_SG(**kw)
    t0 = time.perf_counter()
    FS.SGWB_iter_fast(fast, kink_split=True, freq_grid='goal',
                      frequency_quadrature='pchip')
    fast_runtime = time.perf_counter() - t0
    freqs = np.asarray(fast.f, dtype=float)
    if freq_min is not None or freq_max is not None:
        lo = -np.inf if freq_min is None else float(freq_min)
        hi = np.inf if freq_max is None else float(freq_max)
        freqs = freqs[(freqs >= lo) & (freqs <= hi)]
    if freq_count is not None:
        if freq_count < 3 or freq_count > freqs.size:
            raise ValueError('freq_count must be between 3 and the native grid size')
        order = np.argsort(freqs)
        picks = np.linspace(0, freqs.size - 1, freq_count, dtype=int)
        freqs = freqs[order[picks]]
    dn_eff = float(fast.cosmo_param['DN_eff'])

    rows = []
    for z_tail in Z_TAILS:
        key = _cache_key(point, kw, freqs, dn_eff, z_tail, rtol)
        point_records = checkpoint.setdefault('records', {}).setdefault(point, {}) if checkpoint is not None else {}
        cached = point_records.get(str(z_tail))
        if cached is not None and cached.get('cache_key') == key:
            row = dict(cached)
            row['cache_hit'] = True
            rows.append(row)
            continue
        ref_model = LCDM_SG(**kw)
        t0 = time.perf_counter()
        result = REF.run_reference(
            ref_model, dn_eff=dn_eff, freq_res=1.0, z_tail=z_tail,
            rtol=rtol, freq_subset=freqs, self_consistent=False,
            workers=workers)
        runtime = time.perf_counter() - t0
        rows.append({
            'z_tail': z_tail,
            'DN_gw': float(result['DN_gw']),
            'DN_eff': float(result['DN_eff']),
            'runtime_s': runtime,
            'used_tail_fraction': float(np.mean(result['used_tail'])),
            'quadrature_error': float(result['quadrature_error']),
            'interpolation_error': (
                None if result['interpolation_error'] is None
                else float(result['interpolation_error'])),
            'n_freq': int(result['n_freq']),
            'cache_key': key,
            'cache_hit': False,
        })
        if checkpoint is not None:
            point_records[str(z_tail)] = rows[-1]
            _save_checkpoint(checkpoint_path, checkpoint)

    summary = REF.summarize_tail_convergence(rows)
    summary.update({
        'point': point,
        'parameters': kw,
        'fast_DN_eff': dn_eff,
        'fast_DN_gw': float(fast.DN_gw[-1]),
        'fast_runtime_s': fast_runtime,
        'rtol': rtol,
        'n_freq': int(freqs.size),
        'frequency_grid': ('same native fast goal grid' if freq_count is None
                           and freq_min is None and freq_max is None
                           else 'representative subset of native grid'),
        'frequency_subset_requested': freq_count,
        'frequency_band_requested': [freq_min, freq_max],
        'oracle_semantics': 'fixed DN_eff; only reference z_tail varies',
        'acceptance_note': (
            'The observed deepest-tail difference is the reported systematic '
            'bound; fitted decay is descriptive and not used to shrink it.'),
        'cache_hits': int(sum(row.get('cache_hit', False) for row in rows)),
        'cache_keys': [row['cache_key'] for row in rows],
    })
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', default=['default'],
                        choices=sorted(CASES))
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--workers', type=int, default=1,
                        help='reference workers; keep at 1 unless explicitly raised')
    parser.add_argument('--freq-count', type=int, default=None,
                        help='optional representative subset size for a quick diagnostic')
    parser.add_argument('--freq-min', type=float, default=None)
    parser.add_argument('--freq-max', type=float, default=None)
    parser.add_argument('--out', default='docs/oracle_tail_convergence_head.json')
    parser.add_argument('--checkpoint', default=None,
                        help='逐个 z_tail 原子写入的 checkpoint 路径')
    parser.add_argument('--resume', action='store_true',
                        help='从匹配当前 commit/schema 的 checkpoint 恢复')
    args = parser.parse_args(argv)

    checkpoint_path = args.checkpoint or (args.out + '.checkpoint.json')
    checkpoint = _load_checkpoint(checkpoint_path) if args.resume else _new_checkpoint()
    records = [run_point(point, rtol=args.rtol, workers=args.workers,
                         freq_count=args.freq_count, freq_min=args.freq_min,
                         freq_max=args.freq_max, checkpoint=checkpoint,
                         checkpoint_path=checkpoint_path)
               for point in args.points]
    payload = {
        'schema_version': SCHEMA_VERSION,
        'generated_commit': _commit(),
        'reference_version': _reference_version(),
        'resources': telemetry(workers=args.workers, threads=2),
        'z_tails': list(Z_TAILS),
        'checkpoint': checkpoint_path,
        'points': records,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    print('wrote %s' % args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
