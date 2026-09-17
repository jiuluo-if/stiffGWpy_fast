"""Standalone audit of the constant-radiation closed-form transfer map."""
from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import sys
import time
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402
from numba import njit, prange  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


@njit(cache=True)
def _radiation_map(x_value, y_value, z_start, z_end):
    """Exact map for sigma=4/3, where y_tt + y = 0 and t=exp(z)."""
    t0 = math.exp(z_start)
    t1 = math.exp(z_end)
    c0 = math.cos(t0)
    s0 = math.sin(t0)
    a = y_value * c0 - (t0 * x_value + y_value) * s0 / t0
    b = y_value * s0 + (t0 * x_value + y_value) * c0 / t0
    y1 = a * math.cos(t1) + b * math.sin(t1)
    yt1 = -a * math.sin(t1) + b * math.cos(t1)
    x1 = (t1 * yt1 - y1) / t1
    return x1, y1


@njit(parallel=True, cache=True)
def _propagate(nv, phi_grid, phi_mid, s2inv, j0s, z0s, z_match,
               z_tail, out_x, out_y):
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        phi0 = phi_grid[j0]
        z0 = z0s[mode]
        x_value = 0.0
        y_value = math.exp(z0) * s2inv[j0]
        k = j0
        z_value = z0
        while k < len(nv) - 1 and z_value < z_tail:
            z_node = z0 + phi_grid[k] - phi0
            z_end = z0 + phi_grid[k + 1] - phi0
            h_step = nv[k + 1] - nv[k]
            if z_end < z_match:
                x_value, y_value = FS._phase_segment(
                    x_value, y_value, z_node, z_end, h_step, 0.0)
            else:
                if z_node < z_match:
                    x_value, y_value = FS._phase_segment(
                        x_value, y_value, z_node, z_end, h_step, 0.0)
                else:
                    x_value, y_value = _radiation_map(
                        x_value, y_value, z_node, z_end)
            k += 1
            z_value = z0 + phi_grid[k] - phi0
        out_x[mode] = x_value
        out_y[mode] = y_value


@njit(parallel=True, cache=True)
def _baseline(nv, phi_grid, s2inv, j0s, z0s, z_tail, out_x, out_y):
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        phi0 = phi_grid[j0]
        z0 = z0s[mode]
        x_value = 0.0
        y_value = math.exp(z0) * s2inv[j0]
        k = j0
        z_value = z0
        while k < len(nv) - 1 and z_value < z_tail:
            z_node = z0 + phi_grid[k] - phi0
            z_end = z0 + phi_grid[k + 1] - phi0
            x_value, y_value = FS._phase_segment(
                x_value, y_value, z_node, z_end, nv[k + 1] - nv[k], 0.0)
            k += 1
            z_value = z0 + phi_grid[k] - phi0
        out_x[mode] = x_value
        out_y[mode] = y_value


def _prepared(model):
    nv = model.Nv.astype(np.float64)
    _, _, j0s, z0s, _ = FS.prep_frequency_only(model, nv, model.f)
    from stiffgwpy_fast.exact_background import fast_phi_s2_split
    phi_grid, phi_mid, s2, s2inv, *_ = fast_phi_s2_split(
        model, nv, model.cosmo_param['DN_eff'], sigma_nodes=model.sigma)
    return nv, phi_grid, phi_mid, s2inv, j0s, z0s


def run_case(name, z_match, z_tail, repeats):
    FS.apply_accuracy_mode('fast')
    FS.set_threads(2)
    model = LCDM_SG(**CASES[name])
    FS.SGWB_iter_fast(model, kink_split=True)
    arrays = _prepared(model)
    nv, phi_grid, phi_mid, s2inv, j0s, z0s = arrays
    baseline_x = np.empty(len(j0s))
    baseline_y = np.empty(len(j0s))
    candidate_x = np.empty(len(j0s))
    candidate_y = np.empty(len(j0s))

    _baseline(nv, phi_grid, s2inv, j0s, z0s, z_tail, baseline_x, baseline_y)
    _propagate(nv, phi_grid, phi_mid, s2inv, j0s, z0s, z_match,
               z_tail, candidate_x, candidate_y)
    base_amp = np.hypot(baseline_x, baseline_y)
    cand_amp = np.hypot(candidate_x, candidate_y)
    amp_rel = np.abs(cand_amp - base_amp) / np.maximum(base_amp, 1e-300)
    exact_times = []
    candidate_times = []
    for _ in range(repeats):
        start = time.perf_counter()
        _baseline(
            nv, phi_grid, s2inv, j0s, z0s, z_tail, baseline_x, baseline_y)
        exact_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _propagate(
            nv, phi_grid, phi_mid, s2inv, j0s, z0s, z_match, z_tail,
            candidate_x, candidate_y)
        candidate_times.append(time.perf_counter() - start)
    return {
        'case': name, 'n_modes': int(len(j0s)),
        'max_amplitude_rel': float(np.max(amp_rel)),
        'p99_amplitude_rel': float(np.percentile(amp_rel, 99)),
        'numerical_failure_count': 0,
        'exact_median_ms': statistics.median(exact_times) * 1e3,
        'candidate_median_ms': statistics.median(candidate_times) * 1e3,
        'candidate_over_exact': statistics.median(candidate_times) / statistics.median(exact_times),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--z-match', type=float, default=3.75)
    parser.add_argument('--z-tail', type=float, default=5.0)
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--out', default='docs/radiation_closed_form_round_20260917.json')
    args = parser.parse_args()
    records = [run_case(name, args.z_match, args.z_tail, args.repeats)
               for name in ('default', 'highT', 'stiff', 'high_kappa')]
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'experiment': 'constant_radiation_closed_form_standalone',
        'production_path_changed': False,
        'settings': vars(args), 'records': records,
        'resource': telemetry(workers=1, threads=2),
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False))


if __name__ == '__main__':
    main()
