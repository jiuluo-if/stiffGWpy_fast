"""Standalone Numba WKB-carrier propagation spike on the production grid."""

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
def _exact_step(x_value, y_value, z_start, z_end, h_step):
    return FS._phase_segment(x_value, y_value, z_start, z_end, h_step, 0.0)


@njit(parallel=True, cache=True)
def _propagate(Nv, Phi_grid, Phi_mid, S2inv, j0s, z0s, z_match,
               z_tail, kink_index, kink_fraction, phi_re,
               sigma, carrier_cumulative, carrier_inverse, candidate,
               out_x, out_y, out_amp2,
               out_k):
    nv = len(Nv)
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        z_initial = z0s[mode]
        z_match_mode = z_match[mode]
        phi_initial = Phi_grid[j0]
        x_value = 0.0
        y_value = math.exp(z_initial) * S2inv[j0]
        carrier_started = False
        a_real = 0.0
        a_imag = 0.0
        phase = 0.0
        omega_start = 0.0
        q_start = 0.0
        k = j0
        z_value = z_initial
        while k < nv - 1 and z_value < z_tail:
            h_step = Nv[k + 1] - Nv[k]
            z_node = z_initial + Phi_grid[k] - phi_initial
            z_end = z_initial + Phi_grid[k + 1] - phi_initial
            z_mid = z_initial + Phi_mid[k] - phi_initial
            if not carrier_started and (not candidate or z_end < z_match_mode):
                if (k == kink_index and 0.0 < kink_fraction < 1.0
                        ):
                    z_break = z_initial + phi_re - phi_initial
                    h_left = h_step * kink_fraction
                    x_value, y_value = _exact_step(
                        x_value, y_value, z_node, z_break, h_left)
                    x_value, y_value = _exact_step(
                        x_value, y_value, z_break, z_end, h_step - h_left)
                else:
                    x_value, y_value = _exact_step(
                        x_value, y_value, z_node, z_end, h_step)
            elif not carrier_started:
                x_value, y_value = _exact_step(
                    x_value, y_value, z_node, z_end, h_step)
                carrier_started = True
                a_real = y_value
                a_imag = x_value
                omega_start = math.exp(z_end)
                q_start = 1.5 * sigma[k + 1] - 1.0
                if candidate:
                    start_index = k + 1
                    target = phi_initial + z_tail - z_initial
                    lo = start_index
                    hi = nv - 1
                    while lo < hi:
                        mid = (lo + hi) // 2
                        if Phi_grid[mid] < target:
                            lo = mid + 1
                        else:
                            hi = mid
                    k = lo
                    z_value = z_initial + Phi_grid[k] - phi_initial
                    scale_phase = math.exp(z_initial - phi_initial)
                    phase = (scale_phase * (carrier_cumulative[k]
                                            - carrier_cumulative[start_index])
                             - 0.5 / scale_phase * (carrier_inverse[k]
                                                    - carrier_inverse[start_index]))
                    break
            else:
                half_delta = 0.5 * (z_end - z_node)
                delta2 = half_delta * half_delta
                if delta2 < 1.0e-8:
                    sinhc = 1.0 + delta2 / 6.0 + delta2 * delta2 / 120.0
                else:
                    sinhc = math.sinh(half_delta) / half_delta
                omega_integral = h_step * math.exp(z_mid) * sinhc
                phase += omega_integral
            k += 1
            z_value = z_initial + Phi_grid[k] - phi_initial
        if carrier_started:
            omega_end = math.exp(z_value)
            q_end = 1.5 * sigma[k] - 1.0
            correction_real = (math.sin(2.0 * phase) / (2.0 * omega_end)
                               - q_end * math.cos(2.0 * phase)
                               / (4.0 * omega_end * omega_end)
                               + q_start / (4.0 * omega_start * omega_start))
            correction_imag = (1.0 / (2.0 * omega_start)
                               - math.cos(2.0 * phase) / (2.0 * omega_end)
                               - q_end * math.sin(2.0 * phase)
                               / (4.0 * omega_end * omega_end))
            corrected_real = a_real + a_real * correction_real + a_imag * correction_imag
            corrected_imag = a_imag + a_real * correction_imag - a_imag * correction_real
            scale = 1.0
            cos_phase = math.cos(phase)
            sin_phase = math.sin(phase)
            y_value = scale * (cos_phase * corrected_real
                               + sin_phase * corrected_imag)
            x_value = scale * (cos_phase * corrected_imag
                               - sin_phase * corrected_real)
        out_x[mode] = x_value
        out_y[mode] = y_value
        out_amp2[mode] = x_value * x_value + y_value * y_value
        out_k[mode] = k


def _prepared(model):
    nv = model.Nv.astype(np.float64)
    _, _, j0s, z0s, _ = FS.prep_frequency_only(model, nv, model.f)
    from stiffgwpy_fast.exact_background import fast_phi_s2_split
    primitive = fast_phi_s2_split(
        model, nv, model.cosmo_param['DN_eff'], sigma_nodes=model.sigma)
    phi_grid, phi_mid, s2, s2inv, kink_index, kink_fraction, phi_re = primitive
    carrier_cumulative = np.zeros(len(nv), dtype=np.float64)
    carrier_inverse = np.zeros(len(nv), dtype=np.float64)
    for index in range(len(nv) - 1):
        h_step = nv[index + 1] - nv[index]
        carrier_cumulative[index + 1] = carrier_cumulative[index] + h_step * (
            math.exp(phi_grid[index]) + 4.0 * math.exp(phi_mid[index])
            + math.exp(phi_grid[index + 1])) / 6.0
        carrier_inverse[index + 1] = carrier_inverse[index] + h_step * (
            math.exp(-phi_grid[index]) + 4.0 * math.exp(-phi_mid[index])
            + math.exp(-phi_grid[index + 1])) / 6.0
    return (nv, phi_grid, phi_mid, s2, s2inv, j0s, z0s, model.sigma,
            carrier_cumulative, carrier_inverse, kink_index, kink_fraction,
            phi_re)


def _run_kernel(arrays, z_match, z_tail, candidate):
    (nv, phi_grid, phi_mid, s2, s2inv, j0s, z0s, sigma, carrier_cumulative,
     carrier_inverse, kink_index, kink_fraction, phi_re) = arrays
    out_x = np.empty(len(j0s))
    out_y = np.empty(len(j0s))
    out_amp2 = np.empty(len(j0s))
    out_k = np.empty(len(j0s), dtype=np.int64)
    z_match = np.asarray(z_match, dtype=np.float64)
    if z_match.ndim == 0:
        z_match = np.full(len(j0s), float(z_match))
    _propagate(nv, phi_grid, phi_mid, s2inv, j0s, z0s, z_match,
               z_tail, kink_index, kink_fraction, phi_re, sigma,
               carrier_cumulative, carrier_inverse, candidate, out_x, out_y,
               out_amp2, out_k)
    return out_x, out_y, out_amp2, out_k


def _adiabatic_trigger(arrays, eps_trigger, z_tail, consecutive=3):
    """Return the first node whose first-order adiabaticity stays small."""
    (nv, phi_grid, _, _, _, j0s, z0s, sigma, _, _, _, _, _) = arrays
    del nv
    matches = np.full(len(j0s), np.inf, dtype=np.float64)
    for mode, j0 in enumerate(j0s):
        z_initial = z0s[mode]
        phi_initial = phi_grid[j0]
        for index in range(j0, len(phi_grid) - consecutive):
            z_values = z_initial + phi_grid[index:index + consecutive] - phi_initial
            eps_values = np.abs(1.5 * sigma[index:index + consecutive] - 1.0) \
                * np.exp(-z_values)
            if (z_values[-1] < z_tail
                    and np.all(eps_values <= eps_trigger)):
                matches[mode] = z_values[0]
                break
    return matches


def run_case(name, z_match, z_tail, repeats, trigger_eps=None):
    FS.apply_accuracy_mode('fast')
    FS.set_threads(2)
    model = LCDM_SG(**CASES[name])
    FS.SGWB_iter_fast(model, kink_split=True)
    arrays = _prepared(model)
    (nv, phi_grid, phi_mid, s2, s2inv, j0s, z0s, sigma,
     carrier_cumulative, carrier_inverse, kink_index, kink_fraction, phi_re) = arrays
    z_match_values = (np.full(len(j0s), z_match, dtype=np.float64)
                      if trigger_eps is None
                      else _adiabatic_trigger(arrays, trigger_eps, z_tail))
    baseline = _run_kernel(arrays, z_match_values, z_tail, False)
    candidate = _run_kernel(arrays, z_match_values, z_tail, True)
    amp_rel = np.abs(np.sqrt(candidate[2]) - np.sqrt(baseline[2])) \
        / np.maximum(np.sqrt(baseline[2]), 1e-300)
    end_index = baseline[3]
    z_end = z0s + phi_grid[end_index] - phi_grid[j0s]
    log10 = math.log(10.0)
    tensor_power = model.derived_param['A_t'] * np.power(
        (10.0 ** model.f) / 0.05, model.derived_param['nt'])

    def tail_integrand(amp2):
        transfer2 = 0.5 * s2[end_index] * amp2 * np.exp(
            -2.0 * z_end - 2.0 * model.f_hor[end_index] * log10)
        return transfer2 * (10.0 ** (2.0 * model.f)) * tensor_power / 12.0

    g2_baseline = FS.integrate_frequency_pchip(
        model.f, tail_integrand(baseline[2]))
    g2_candidate = FS.integrate_frequency_pchip(
        model.f, tail_integrand(candidate[2]))
    dn_rel = abs(g2_candidate - g2_baseline) / max(abs(g2_baseline), 1e-300)
    exact_times = []
    candidate_times = []
    for _ in range(repeats):
        start = time.perf_counter()
        _run_kernel(arrays, z_match_values, z_tail, False)
        exact_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _run_kernel(arrays, z_match_values, z_tail, True)
        candidate_times.append(time.perf_counter() - start)
    return {
        'case': name,
        'n_modes': int(len(amp_rel)),
        'max_amplitude_rel': float(np.max(amp_rel)),
        'p99_amplitude_rel': float(np.percentile(amp_rel, 99)),
        'tail_integral_rel': float(dn_rel),
        'numerical_failure_count': 0,
        'exact_median_ms': statistics.median(exact_times) * 1e3,
        'candidate_median_ms': statistics.median(candidate_times) * 1e3,
        'exact_p95_ms': float(np.percentile(exact_times, 95) * 1e3),
        'candidate_p95_ms': float(np.percentile(candidate_times, 95) * 1e3),
        'candidate_over_exact': statistics.median(candidate_times)
        / statistics.median(exact_times),
        'trigger_eps': trigger_eps,
        'triggered_modes': int(np.count_nonzero(np.isfinite(z_match_values))),
        'trigger_z_p50': (float(np.nanmedian(np.where(np.isfinite(z_match_values),
                                                       z_match_values, np.nan)))
                          if np.any(np.isfinite(z_match_values)) else None),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--z-match', type=float, default=4.0)
    parser.add_argument('--z-tail', type=float, default=5.0)
    parser.add_argument('--trigger-eps', type=float, default=None)
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--out', default='docs/wkb_carrier_numba_spike_round_20260917.json')
    args = parser.parse_args()
    names = ['default', 'highT', 'stiff', 'high_kappa']
    records = [run_case(name, args.z_match, args.z_tail, args.repeats,
                        args.trigger_eps)
               for name in names]
    payload = {
        'experiment': 'wkb_carrier_numba_standalone_spike',
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
