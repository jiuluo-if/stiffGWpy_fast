"""Measure H2 interpolation vectorization on the transition-refined solver path.

中文：在转折细化求解路径上对比逐点与批量 H2 插值，并核验所有输出数组完全一致。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import sys
import time

import numpy as np

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:
    from _resource_budget import apply_environment, telemetry

apply_environment()

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import exact_background as EB  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

OUTPUT_FIELDS = (
    'f', 'DN_gw', 'g2', 'w2', 'N', 'Nv', 'sigma', 'f_hor',
    'Ogw_today', 'Oj_today', 'Opgw_today', 'log10OmegaGW',
)


def _scalar_h2_reference(N, model, DN_eff):
    """Pre-vectorization H2 loop retained as the exact A/B baseline.

    中文：保留向量化前的逐点公式作为 A/B 基线；所有结果必须与批量实现逐位一致。
    """
    N = np.asarray(N, dtype=float)
    d = model.derived_param
    Omh2 = d['Omega_mh2']
    Osh2 = d['Omega_sh2']
    Oerh2 = gp.Omega_ph2 * 7 / 8 * (4 / 11) ** (4 / 3) * DN_eff
    Otrh2 = gp.Omega_orh2 + Oerh2
    Otreh2 = gp.Omega_ph2 * gp.rho_th[-1] + Oerh2
    OLh2 = (d['h'] ** 2 - Omh2 - gp.Omega_mnuh2 - gp.Omega_nh2 * 2 / 3
            - gp.Omega_ph2 - Oerh2 - Osh2)
    N_inf = d['N_inf']
    fd_rho, _ = EB._fd_from_ref()
    out = np.empty_like(N)
    for i in range(N.size):
        Ni = float(N[i])
        eN = np.exp(N_inf - Ni)
        e3N = eN * eN * eN
        nu = gp.nu_today / eN
        if nu > 100.0:
            out[i] = (Omh2 + gp.Omega_mnuh2
                      + (gp.Omega_ph2 + 2.0 / 3.0 * gp.Omega_nh2 + Oerh2) * eN
                      + Osh2 * e3N + OLh2 / e3N)
        elif nu >= 0.1:
            rho_nu = float(fd_rho(np.log10(nu)))
            out[i] = (Omh2
                      + (gp.Omega_ph2 + (2.0 / 3.0 + rho_nu / 3.0) * gp.Omega_nh2 + Oerh2) * eN
                      + Osh2 * e3N + OLh2 / e3N)
        elif Ni > N_inf - gp.N_fin:
            out[i] = Omh2 + Otrh2 * eN + Osh2 * e3N + OLh2 / e3N
        elif Ni >= N_inf - gp.N_max:
            Nl = N_inf - Ni
            out[i] = (Omh2 + (gp.Omega_ph2 * float(gp.spl_rho(Nl)) + Oerh2) * eN
                      + Osh2 * e3N + OLh2 / e3N)
        else:
            out[i] = Omh2 + Otreh2 * eN + Osh2 * e3N + OLh2 / e3N
    return out


def _digest(model):
    payload = b''.join(
        np.ascontiguousarray(np.asarray(getattr(model, name), dtype=np.float64)).tobytes()
        for name in OUTPUT_FIELDS
    )
    return hashlib.sha256(payload).hexdigest()


def _run_once(point, config, h2_function):
    original = EB.H2_vec
    EB.H2_vec = h2_function
    try:
        model = LCDM_SG(**CASES[point])
        started = time.perf_counter()
        result = FS.SGWB_iter_fast(
            model,
            config=config,
            kink_split=False,
            transition_refine=True,
            freq_grid='goal',
            freq_res=1.0,
            tol=1e-6,
            frequency_quadrature='pchip',
        )
        elapsed = time.perf_counter() - started
    finally:
        EB.H2_vec = original
    if result is None or not model.SGWB_converge:
        reason = getattr(model, 'fast_failure_reason', None) or 'not converged'
        raise RuntimeError('transition-refined solve failed for %s: %s' % (point, reason))
    return {
        'seconds': elapsed,
        'DN_eff': float(model.cosmo_param['DN_eff']),
        'nodes': int(model.Nv.size),
        'digest': _digest(model),
        'outputs': {
            name: np.asarray(getattr(model, name), dtype=float).copy()
            for name in OUTPUT_FIELDS
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+',
                        choices=('edge_dnre_lo', 'edge_dnre_hi'),
                        default=['edge_dnre_lo', 'edge_dnre_hi'])
    parser.add_argument('--warmup', type=int, default=2)
    parser.add_argument('--repeats', type=int, default=10)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--out', default=None)
    args = parser.parse_args(argv)
    if args.warmup < 0 or args.repeats < 1:
        parser.error('--warmup must be >= 0 and --repeats must be >= 1')

    config = FS.resolve_config(
        'fast', h=0.005, col_step=8, z_tail=5.0, phase_max=0.25,
        freq_grid='goal', kink_split=False, threads=args.threads,
    )
    methods = {
        'scalar_reference': _scalar_h2_reference,
        'vectorized': EB.H2_vec,
    }
    points = {}
    for point in args.points:
        measured = {name: [] for name in methods}
        digests = {name: [] for name in methods}
        latest = {}
        for method, function in methods.items():
            for _ in range(args.warmup):
                _run_once(point, config, function)
        for repeat in range(args.repeats):
            order = tuple(methods) if repeat % 2 == 0 else tuple(reversed(methods))
            for method in order:
                result = _run_once(point, config, methods[method])
                measured[method].append(result['seconds'])
                digests[method].append(result['digest'])
                latest[method] = result

        reference = latest['scalar_reference']
        candidate = latest['vectorized']
        exact = (
            reference['DN_eff'] == candidate['DN_eff']
            and reference['nodes'] == candidate['nodes']
            and all(np.array_equal(reference['outputs'][name], candidate['outputs'][name])
                    for name in OUTPUT_FIELDS)
        )
        if not exact:
            raise RuntimeError('vectorized H2 changed solver outputs for %s' % point)
        if len(set(digests['scalar_reference'] + digests['vectorized'])) != 1:
            raise RuntimeError('solver digest changed across repeats for %s' % point)
        baseline = statistics.median(measured['scalar_reference'])
        optimized = statistics.median(measured['vectorized'])
        points[point] = {
            'scalar_reference_runtime_s': measured['scalar_reference'],
            'scalar_reference_median_s': baseline,
            'vectorized_runtime_s': measured['vectorized'],
            'vectorized_median_s': optimized,
            'speedup': baseline / optimized,
            'nodes': candidate['nodes'],
            'DN_eff': candidate['DN_eff'],
            'output_digest': candidate['digest'],
            'outputs_exact': exact,
        }
        print('%-14s scalar=%7.3f ms vector=%7.3f ms speedup=%.2fx exact=%s' % (
            point, baseline * 1e3, optimized * 1e3, baseline / optimized, exact))

    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(workers=1, threads=args.threads),
        'configuration': {
            'profile': 'fast',
            'path': 'transition-refine',
            'h': config.h,
            'col_step': config.col_step,
            'z_tail': config.z_tail,
            'phase_max': config.phase_max,
            'freq_grid': config.freq_grid,
            'freq_res': 1.0,
            'tol': 1e-6,
            'frequency_quadrature': 'pchip',
            'warmup': args.warmup,
            'repeats': args.repeats,
        },
        'points': points,
        'semantics': ('full transition-refined SGWB runtime A/B; the scalar '
                      'reference preserves the pre-vectorization H2 loop and '
                      'all solver output arrays must match exactly'),
    }
    if args.out:
        with open(args.out, 'w', encoding='utf-8') as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
        print('wrote', args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
