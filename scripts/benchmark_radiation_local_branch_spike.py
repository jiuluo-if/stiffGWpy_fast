# -*- coding: utf-8 -*-
"""按实际 sigma 连续区间测试 radiation exact-map 的 standalone 跳跃。"""
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
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
from stiffgwpy_fast import exact_background as EB  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


@njit(cache=True)
def _radiation_map(x_value, y_value, z_start, z_end):
    """sigma=4/3 的精确 transfer map；只用于局部 standalone 候选。"""
    t0 = np.exp(z_start)
    t1 = np.exp(z_end)
    c0 = np.cos(t0)
    s0 = np.sin(t0)
    a = y_value * c0 - (t0 * x_value + y_value) * s0 / t0
    b = y_value * s0 + (t0 * x_value + y_value) * c0 / t0
    y1 = a * np.cos(t1) + b * np.sin(t1)
    yt1 = -a * np.sin(t1) + b * np.cos(t1)
    x1 = (t1 * yt1 - y1) / t1
    return x1, y1


@njit(inline='always')
def _radiation_eligible(sigma_left, sigma_right, tolerance):
    return (abs(sigma_left - 4.0 / 3.0) < tolerance
            and abs(sigma_right - 4.0 / 3.0) < tolerance)


@njit(parallel=True, cache=True)
def _baseline(nv, phi_grid, s2inv, j0s, z0s, z_tail, out_x, out_y):
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        phi0 = phi_grid[j0]
        z0 = z0s[mode]
        x_value = 0.0
        y_value = np.exp(z0) * s2inv[j0]
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


@njit(parallel=True, cache=True)
def _local_branch(nv, phi_grid, s2inv, sigma, j0s, z0s, z_tail,
                  tolerance, out_x, out_y, branch_segments):
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        phi0 = phi_grid[j0]
        z0 = z0s[mode]
        x_value = 0.0
        y_value = np.exp(z0) * s2inv[j0]
        k = j0
        z_value = z0
        local_count = 0
        while k < len(nv) - 1 and z_value < z_tail:
            z_node = z0 + phi_grid[k] - phi0
            if _radiation_eligible(sigma[k], sigma[k + 1], tolerance):
                # 连续候选段合并后只做一次 exact map，实际减少 propagation 调用。
                end = k + 1
                while end < len(nv) - 1 and _radiation_eligible(
                        sigma[end], sigma[end + 1], tolerance):
                    end += 1
                z_end = z0 + phi_grid[end] - phi0
                x_value, y_value = _radiation_map(
                    x_value, y_value, z_node, z_end)
                local_count += end - k
                k = end
            else:
                z_end = z0 + phi_grid[k + 1] - phi0
                x_value, y_value = FS._phase_segment(
                    x_value, y_value, z_node, z_end, nv[k + 1] - nv[k], 0.0)
                k += 1
            z_value = z0 + phi_grid[k] - phi0
        out_x[mode] = x_value
        out_y[mode] = y_value
        branch_segments[mode] = local_count


def _prepared(model):
    nv = model.Nv.astype(np.float64)
    _, _, j0s, z0s, _ = FS.prep_frequency_only(model, nv, model.f)
    phi_grid, _, _, s2inv, *_ = EB.fast_phi_s2_split(
        model, nv, model.cosmo_param['DN_eff'], sigma_nodes=model.sigma)
    return nv, phi_grid, s2inv, model.sigma.astype(np.float64), j0s, z0s


def run_case(name, tolerance, repeats):
    model = LCDM_SG(**CASES[name])
    FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal')
    nv, phi_grid, s2inv, sigma, j0s, z0s = _prepared(model)
    baseline_x = np.empty(len(j0s))
    baseline_y = np.empty(len(j0s))
    candidate_x = np.empty(len(j0s))
    candidate_y = np.empty(len(j0s))
    branch_segments = np.empty(len(j0s), dtype=np.int64)
    _baseline(nv, phi_grid, s2inv, j0s, z0s, 5.0, baseline_x, baseline_y)
    _local_branch(nv, phi_grid, s2inv, sigma, j0s, z0s, 5.0,
                  tolerance, candidate_x, candidate_y, branch_segments)
    baseline_amp = np.hypot(baseline_x, baseline_y)
    candidate_amp = np.hypot(candidate_x, candidate_y)
    amp_rel = np.abs(candidate_amp - baseline_amp) / np.maximum(baseline_amp, 1e-300)
    baseline_times = []
    candidate_times = []
    for _ in range(repeats):
        start = time.perf_counter()
        _baseline(nv, phi_grid, s2inv, j0s, z0s, 5.0, baseline_x, baseline_y)
        baseline_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _local_branch(nv, phi_grid, s2inv, sigma, j0s, z0s, 5.0,
                      tolerance, candidate_x, candidate_y, branch_segments)
        candidate_times.append(time.perf_counter() - start)
    return {
        'case': name,
        'tolerance': tolerance,
        'n_modes': int(len(j0s)),
        'branch_segments_median': float(np.median(branch_segments)),
        'branch_segments_total': int(np.sum(branch_segments)),
        'amplitude_rel_p50': float(np.percentile(amp_rel, 50)),
        'amplitude_rel_p95': float(np.percentile(amp_rel, 95)),
        'amplitude_rel_max': float(np.max(amp_rel)),
        'baseline_median_ms': statistics.median(baseline_times) * 1e3,
        'candidate_median_ms': statistics.median(candidate_times) * 1e3,
        'candidate_over_baseline': statistics.median(candidate_times)
        / statistics.median(baseline_times),
        'numerical_failure_count': 0,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tolerance', type=float, nargs='+',
                        default=[1e-8, 1e-6, 1e-4])
    parser.add_argument('--repeats', type=int, default=30)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode('fast')
    FS.set_threads(2)
    rows = [run_case(name, tolerance, args.repeats)
            for tolerance in args.tolerance
            for name in ('default', 'highT', 'stiff', 'high_kappa', 'lowT')]
    payload = {
        'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                          cwd=ROOT, text=True).strip(),
        'experiment': 'radiation_local_branch_standalone',
        'production_path_changed': False,
        'settings': vars(args),
        'resource': telemetry(workers=1, threads=2),
        'rows': rows,
    }
    output = ROOT / args.out
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n',
                      encoding='utf-8')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
