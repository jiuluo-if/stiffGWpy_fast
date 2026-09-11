# -*- coding: utf-8 -*-
"""A/B phase-cap candidates against the formal fast profile."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

try:
    from scripts._resource_budget import apply_environment, limit_affinity, telemetry
except ImportError:
    from _resource_budget import apply_environment, limit_affinity, telemetry

apply_environment()
import numpy as np  # noqa: E402
import psutil  # noqa: E402

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


def digest(value):
    return hashlib.sha256(np.ascontiguousarray(value, dtype=np.float64).tobytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase-max', type=float, nargs='+', default=[0.25, 0.35, 0.5])
    ap.add_argument('--reps', type=int, default=25)
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--out', default='docs/benchmark_phase_candidate_head.json')
    args = ap.parse_args()
    process = psutil.Process()
    limit_affinity(process, args.threads)
    os.environ['FAST_THREADS'] = str(args.threads)
    FS.apply_accuracy_mode('fast')
    FS.set_threads(args.threads)
    rows = []
    baseline = {}
    for phase_max in args.phase_max:
        for name, kw in CASES.items():
            config = FS.FastSolverConfig(h=0.005, col_step=8, z_tail=5.0,
                                         phase_max=phase_max, freq_grid='goal',
                                         threads=args.threads)
            samples = []
            reference = None
            for _ in range(args.reps):
                model = LCDM_SG(**kw)
                start = time.perf_counter()
                result = FS.SGWB_iter_fast(model, config=config, kink_split=True)
                samples.append((time.perf_counter() - start) * 1e3)
                if result is not model:
                    reference = {'failure': getattr(model, 'fast_failure_reason', None)}
                    continue
                state = {
                    'DN_gw': float(model.DN_gw[-1]),
                    'spectrum': np.asarray(model.log10OmegaGW),
                    'digest_spectrum': digest(model.log10OmegaGW),
                    'digest_DN_gw': digest(model.DN_gw),
                    'n_freq': int(model.f.size),
                }
                if reference is None:
                    reference = state
                if phase_max == args.phase_max[0]:
                    baseline[name] = state
            finite = [x for x in samples if np.isfinite(x)]
            row = {
                'case': name,
                'phase_max': phase_max,
                'reps': args.reps,
                'warm_median_ms': float(np.median(finite)),
                'warm_p95_ms': float(np.percentile(finite, 95)),
                'n_freq': reference.get('n_freq') if reference else None,
                'DN_gw': reference.get('DN_gw') if reference else None,
                'digest_spectrum': reference.get('digest_spectrum') if reference else None,
                'digest_DN_gw': reference.get('digest_DN_gw') if reference else None,
                'failure': reference.get('failure') if reference else 'missing',
            }
            if name in baseline and reference and 'spectrum' in reference:
                base = baseline[name]
                delta = np.abs(10.0 ** (reference['spectrum'] - base['spectrum']) - 1.0)
                row['DN_rel_vs_formal'] = abs(
                    reference['DN_gw'] - base['DN_gw']) / max(abs(base['DN_gw']), 1e-300)
                row['spectrum_max_rel_vs_formal'] = float(np.max(delta))
            else:
                row['DN_rel_vs_formal'] = None
                row['spectrum_max_rel_vs_formal'] = None
            rows.append(row)
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(process, threads=args.threads),
        'rows': rows,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
