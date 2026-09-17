"""Fresh attribution of outer, assembly, and preparation work."""
from __future__ import annotations

import json
import os
import statistics
import sys
import time
from collections import defaultdict

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import exact_background as EB  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _median(values):
    return float(statistics.median(values)) if values else 0.0


def run_case(name, repeats=5):
    originals = {
        'gen': FS.gen_fast,
        'correct': FS._correct_kink_background,
        'goal': FA.goal_oriented_freqs,
        'prep': FS.prep_frequency_only,
        'primitive': EB.fast_phi_s2_split,
        'solve': FS.solve_kernel,
        'integrate': FS.int_SGWB_W,
        'pchip': FS._pchip_integrals_vectorized,
        'fine': FS.pchip_fine,
    }
    rows = []
    try:
        # JIT and import warmup are excluded from the reported repeats.
        warm = LCDM_SG(**CASES[name])
        FS.SGWB_iter_fast(warm, kink_split=True)

        def timed(name_key, fn):
            def wrapper(*args, **kwargs):
                start = time.perf_counter()
                result = fn(*args, **kwargs)
                elapsed = (time.perf_counter() - start) * 1e3
                rows[-1][name_key].append(elapsed)
                return result
            return wrapper

        FS.gen_fast = timed('gen_ms', originals['gen'])
        FS._correct_kink_background = timed('correct_ms', originals['correct'])
        FA.goal_oriented_freqs = timed('goal_ms', originals['goal'])
        FS.prep_frequency_only = timed('prep_ms', originals['prep'])
        EB.fast_phi_s2_split = timed('primitive_ms', originals['primitive'])
        FS.int_SGWB_W = timed('integrate_ms', originals['integrate'])
        FS._pchip_integrals_vectorized = timed('pchip_ms', originals['pchip'])
        FS.pchip_fine = timed('fine_ms', originals['fine'])

        def solve(*args, **kwargs):
            assemble = int(args[11])
            start = time.perf_counter()
            result = originals['solve'](*args, **kwargs)
            elapsed = (time.perf_counter() - start) * 1e3
            rows[-1]['solve_ms'].append(elapsed)
            rows[-1]['solve_assemble_ms' if assemble else 'solve_probe_ms'].append(elapsed)
            rows[-1]['assemble_flags'].append(assemble)
            return result

        FS.solve_kernel = solve
        for _ in range(repeats):
            row = defaultdict(list)
            rows.append(row)
            model = LCDM_SG(**CASES[name])
            start = time.perf_counter()
            FS.SGWB_iter_fast(model, kink_split=True)
            row['total_ms'].append((time.perf_counter() - start) * 1e3)
            row['iterations'] = len(row['solve_ms'])
            row['converged'] = bool(getattr(model, 'SGWB_converge', False))
            row['failure'] = getattr(model, 'fast_failure_reason', None)
    finally:
        FS.gen_fast = originals['gen']
        FS._correct_kink_background = originals['correct']
        FA.goal_oriented_freqs = originals['goal']
        FS.prep_frequency_only = originals['prep']
        EB.fast_phi_s2_split = originals['primitive']
        FS.solve_kernel = originals['solve']
        FS.int_SGWB_W = originals['integrate']
        FS._pchip_integrals_vectorized = originals['pchip']
        FS.pchip_fine = originals['fine']

    summary = {
        'case': name,
        'repeats': repeats,
        'total_median_ms': _median([r['total_ms'][0] for r in rows]),
        'stages_median_ms': {
            key: _median([value for row in rows for value in row[key]])
            for key in ('gen_ms', 'correct_ms', 'goal_ms', 'prep_ms',
                        'primitive_ms', 'solve_ms', 'solve_probe_ms',
                        'solve_assemble_ms', 'integrate_ms', 'pchip_ms',
                        'fine_ms')
        },
        'assemble_flags': [flag for row in rows for flag in row['assemble_flags']],
        'converged': all(row['converged'] for row in rows),
        'failures': [row['failure'] for row in rows if row['failure'] is not None],
    }
    return summary


def main():
    FS.apply_accuracy_mode('fast')
    FS.set_threads(2)
    rows = [run_case(name) for name in ('default', 'highT', 'stiff', 'high_kappa', 'lowT')]
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'experiment': 'outer_assembly_attribution',
        'production_path_changed': False,
        'resources': telemetry(workers=1, threads=2),
        'rows': rows,
    }
    out = 'docs/outer_attribution_round18_20260917.json'
    with open(out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
