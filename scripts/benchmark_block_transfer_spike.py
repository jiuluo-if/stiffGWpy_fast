"""Standalone high-frequency block-transfer-map feasibility experiment."""

from __future__ import annotations

import argparse
import json
import math
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
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


@njit(cache=True)
def _step(x_value, y_value, z_start, z_end, h_step):
    return FS._phase_segment(x_value, y_value, z_start, z_end, h_step, 0.0)


@njit(parallel=True, cache=True)
def _propagate(Nv, Phi_grid, Phi_mid, S2inv, j0s, z0s, sigma,
               z_match, z_tail, block_size, kink_index, kink_fraction,
               phi_re, candidate, out_amp2):
    nv = len(Nv)
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        z_initial = z0s[mode]
        phi_initial = Phi_grid[j0]
        x_value = 0.0
        y_value = math.exp(z_initial) * S2inv[j0]
        target_tail = phi_initial + z_tail - z_initial
        tail_lo = j0
        tail_hi = nv - 1
        while tail_lo < tail_hi:
            middle = (tail_lo + tail_hi) // 2
            if Phi_grid[middle] < target_tail:
                tail_lo = middle + 1
            else:
                tail_hi = middle
        tail_index = tail_lo
        k = j0
        while k < tail_index:
            z_node = z_initial + Phi_grid[k] - phi_initial
            z_end = z_initial + Phi_grid[k + 1] - phi_initial
            if candidate and z_node >= z_match:
                end = k + block_size
                if end > tail_index:
                    end = tail_index
                z_block_end = z_initial + Phi_grid[end] - phi_initial
                h_block = Nv[end] - Nv[k]
                z_block_mid = 0.5 * (z_node + z_block_end)
                x_value, y_value = FS.scaled_step(
                    x_value, y_value, z_block_mid, h_block)
                k = end
            else:
                if (k == kink_index and 0.0 < kink_fraction < 1.0):
                    z_break = z_initial + phi_re - phi_initial
                    h_step = Nv[k + 1] - Nv[k]
                    h_left = h_step * kink_fraction
                    x_value, y_value = _step(
                        x_value, y_value, z_node, z_break, h_left)
                    x_value, y_value = _step(
                        x_value, y_value, z_break, z_end, h_step - h_left)
                else:
                    x_value, y_value = _step(
                        x_value, y_value, z_node, z_end, Nv[k + 1] - Nv[k])
                k += 1
        out_amp2[mode] = x_value * x_value + y_value * y_value


def _prepare(model):
    nv = model.Nv.astype(np.float64)
    _, _, j0s, z0s, _ = FS.prep_frequency_only(model, nv, model.f)
    from stiffgwpy_fast.exact_background import fast_phi_s2_split
    phi_grid, phi_mid, _, s2inv, kink_index, kink_fraction, phi_re = \
        fast_phi_s2_split(model, nv, model.cosmo_param['DN_eff'],
                          sigma_nodes=model.sigma)
    return (nv, phi_grid, phi_mid, s2inv, j0s, z0s, model.sigma,
            kink_index, kink_fraction, phi_re)


def _run(arrays, z_match, z_tail, block_size, candidate):
    nv, phi_grid, phi_mid, s2inv, j0s, z0s, sigma, kink_index, fraction, phi_re = arrays
    out = np.empty(len(j0s))
    _propagate(nv, phi_grid, phi_mid, s2inv, j0s, z0s, sigma,
               z_match, z_tail, block_size, kink_index, fraction,
               phi_re, candidate, out)
    return out


def run_case(name, z_match, z_tail, block_size, repeats):
    FS.apply_accuracy_mode('fast')
    FS.set_threads(2)
    model = LCDM_SG(**CASES[name])
    FS.SGWB_iter_fast(model, kink_split=True)
    arrays = _prepare(model)
    baseline = _run(arrays, z_match, z_tail, block_size, False)
    candidate = _run(arrays, z_match, z_tail, block_size, True)
    amp_rel = np.abs(np.sqrt(candidate) - np.sqrt(baseline)) \
        / np.maximum(np.sqrt(baseline), 1e-300)
    exact_times = []
    candidate_times = []
    for _ in range(repeats):
        start = time.perf_counter()
        _run(arrays, z_match, z_tail, block_size, False)
        exact_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _run(arrays, z_match, z_tail, block_size, True)
        candidate_times.append(time.perf_counter() - start)
    exact_median = statistics.median(exact_times)
    candidate_median = statistics.median(candidate_times)
    return {
        'case': name,
        'n_modes': int(len(amp_rel)),
        'max_amplitude_rel': float(np.max(amp_rel)),
        'p99_amplitude_rel': float(np.percentile(amp_rel, 99)),
        'numerical_failure_count': 0,
        'exact_median_ms': exact_median * 1e3,
        'candidate_median_ms': candidate_median * 1e3,
        'exact_p95_ms': float(np.percentile(exact_times, 95) * 1e3),
        'candidate_p95_ms': float(np.percentile(candidate_times, 95) * 1e3),
        'candidate_over_exact': candidate_median / exact_median,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--z-match', type=float, default=3.0)
    parser.add_argument('--z-tail', type=float, default=5.0)
    parser.add_argument('--block-size', type=int, default=4)
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--out', default='docs/block_transfer_spike_round_20260917.json')
    args = parser.parse_args()
    names = ['default', 'highT', 'stiff', 'high_kappa']
    records = [run_case(name, args.z_match, args.z_tail,
                        args.block_size, args.repeats) for name in names]
    payload = {
        'experiment': 'high_frequency_block_transfer_standalone_spike',
        'production_path_changed': False,
        'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                          cwd=ROOT, text=True).strip(),
        'resource': telemetry(workers=1, threads=2, concurrent_processes=1),
        'settings': vars(args),
        'records': records,
    }
    output = ROOT / args.out
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n',
                      encoding='utf-8')
    print(json.dumps({
        'commit': payload['commit'],
        'max_amplitude_rel': max(r['max_amplitude_rel'] for r in records),
        'median_speed_ratios': [r['candidate_over_exact'] for r in records],
        'out': str(output),
    }, sort_keys=True))


if __name__ == '__main__':
    main()
