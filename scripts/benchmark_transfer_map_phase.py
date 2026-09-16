"""Standalone linear-z Magnus transfer-map probe for the fast propagation path."""

import argparse
import json
import math
import os
import sys
import time

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402
from scipy.linalg import expm  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.exact_background import fast_phi_s2_split  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _matrix(z):
    w = math.exp(z)
    return np.array([[-1.0, -w], [w, 1.0]], dtype=float)


def magnus_linear_step(xh, yh, z0, z1, h):
    """Apply the second-order Magnus map for linearly interpolated z."""
    a0 = _matrix(z0)
    a1 = _matrix(z1)
    commutator = a1 @ a0 - a0 @ a1
    omega = 0.5 * h * (a0 + a1) + (h * h / 12.0) * commutator
    state = expm(omega) @ np.array([xh, yh], dtype=float)
    return float(state[0]), float(state[1])


def _phase_segment(xh, yh, z0, z1, h, phase_max):
    zmid = 0.5 * (z0 + z1)
    n_sub = 1 if phase_max <= 0.0 or zmid <= 0.0 else max(
        1, int(math.ceil(h * math.exp(zmid) / phase_max)))
    for index in range(n_sub):
        left = z0 + (z1 - z0) * index / n_sub
        right = z0 + (z1 - z0) * (index + 1) / n_sub
        xh, yh = magnus_linear_step(xh, yh, left, right, h / n_sub)
    return xh, yh


def propagate(model, freq, dn_eff, z_tail, candidate, phase_max=0.25):
    nv = np.asarray(model.Nv, dtype=float)
    sigma = np.asarray(model.sigma, dtype=float)
    phi, _, s2, s2inv, kink_index, kink_fraction, phi_re = fast_phi_s2_split(
        model, nv, dn_eff, sigma_nodes=sigma)
    j0s = np.empty(1, dtype=np.int64)
    z0s = np.empty(1, dtype=float)
    fp_minus = np.empty(nv.size, dtype=float)
    FS.prep_frequency_kernel(model.f_hor, np.asarray([freq]), FS.ln10,
                              j0s, z0s, fp_minus)
    j0 = int(j0s[0])
    z0 = float(z0s[0])
    phi0 = phi[j0]
    xh, yh = 0.0, math.exp(z0) * s2inv[j0]
    k = j0
    zz = z0
    last_z = zz
    while k < nv.size - 1 and zz < z_tail:
        z_node = z0 + phi[k] - phi0
        z_end = z0 + phi[k + 1] - phi0
        h_step = float(nv[k + 1] - nv[k])
        if (k == kink_index and 0.0 < kink_fraction < 1.0):
            z_break = z0 + phi_re - phi0
            left_h = h_step * kink_fraction
            if candidate:
                xh, yh = _phase_segment(xh, yh, z_node, z_break,
                                         left_h, phase_max)
                xh, yh = _phase_segment(xh, yh, z_break, z_end,
                                         h_step - left_h, phase_max)
            else:
                xh, yh = FS._phase_segment(xh, yh, z_node, z_break,
                                            left_h, phase_max)
                xh, yh = FS._phase_segment(xh, yh, z_break, z_end,
                                            h_step - left_h, phase_max)
        elif candidate:
            xh, yh = _phase_segment(xh, yh, z_node, z_end,
                                     h_step, phase_max)
        else:
            xh, yh = FS._phase_segment(xh, yh, z_node, z_end,
                                        h_step, phase_max)
        k += 1
        zz = z0 + phi[k] - phi0
        if zz < z_tail:
            last_z = zz
    if zz < z_tail:
        k = nv.size - 1
    else:
        k = max(j0, k - 1)
    return {
        'amplitude_handoff': float(math.sqrt(0.5 * (xh * xh + yh * yh))),
        'phase_handoff': float(math.atan2(xh, yh)),
        'last_z': float(last_z),
        'k': int(k),
    }


def _relative(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', choices=sorted(CASES),
                        default=['default', 'lowT', 'highT', 'stiff', 'high_kappa'])
    parser.add_argument('--freq-count', type=int, default=8)
    parser.add_argument('--z-tail', type=float, default=5.0)
    parser.add_argument('--rtol', type=float, default=1e-9)
    parser.add_argument('--runtime-repeats', type=int, default=25)
    parser.add_argument('--out', default='docs/transfer_map_phase_round.json')
    args = parser.parse_args(argv)
    FS.apply_accuracy_mode('fast')
    rows = []
    probes = []
    for point in args.points:
        model = LCDM_SG(**CASES[point])
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                          frequency_quadrature='pchip')
        freqs = np.sort(np.asarray(model.f, dtype=float))
        picks = np.linspace(0, freqs.size - 1, args.freq_count, dtype=int)
        dn_eff = float(model.cosmo_param['DN_eff'])
        for freq in freqs[picks]:
            ref = REF.solve_reference_mode(model, float(freq), dn_eff,
                                           z_tail=args.z_tail, rtol=args.rtol)
            base = propagate(model, float(freq), dn_eff, args.z_tail, False)
            cand = propagate(model, float(freq), dn_eff, args.z_tail, True)
            rows.append({
                'point': point, 'frequency': float(freq),
                'used_tail': bool(ref['used_tail']),
                'baseline_amplitude_rel_to_oracle': _relative(
                    base['amplitude_handoff'], ref['amplitude_handoff']),
                'candidate_amplitude_rel_to_oracle': _relative(
                    cand['amplitude_handoff'], ref['amplitude_handoff']),
                'candidate_vs_baseline_amplitude_rel': _relative(
                    cand['amplitude_handoff'], base['amplitude_handoff']),
                'baseline_phase_abs_to_oracle': abs(
                    base['phase_handoff'] - ref['phase_handoff']),
                'candidate_phase_abs_to_oracle': abs(
                    cand['phase_handoff'] - ref['phase_handoff']),
                'candidate_vs_baseline_phase_abs': abs(
                    cand['phase_handoff'] - base['phase_handoff']),
                'baseline': base, 'candidate': cand,
            })
            if ref['used_tail'] and not any(item[0] == point for item in probes):
                probes.append((point, model, float(freq), dn_eff))
    tail_rows = [row for row in rows if row['used_tail']]
    summary = {
        'tail_mode_count': len(tail_rows),
        'max_candidate_vs_baseline_amplitude_rel': float(max(
            (row['candidate_vs_baseline_amplitude_rel'] for row in tail_rows),
            default=0.0)),
        'max_candidate_vs_baseline_phase_abs': float(max(
            (row['candidate_vs_baseline_phase_abs'] for row in tail_rows),
            default=0.0)),
        'candidate_amplitude_oracle_better_count': int(sum(
            row['candidate_amplitude_rel_to_oracle']
            < row['baseline_amplitude_rel_to_oracle'] for row in tail_rows)),
    }
    runtime = []
    for point, model, freq, dn_eff in probes:
        for candidate in (False, True):
            for _ in range(2):
                propagate(model, freq, dn_eff, args.z_tail, candidate)
            samples = []
            for _ in range(args.runtime_repeats):
                start = time.perf_counter()
                propagate(model, freq, dn_eff, args.z_tail, candidate)
                samples.append((time.perf_counter() - start) * 1000.0)
            runtime.append({
                'point': point, 'candidate': candidate,
                'repeats': args.runtime_repeats,
                'median_ms': float(np.median(samples)),
                'p95_ms': float(np.percentile(samples, 95)),
            })
    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'z_tail': args.z_tail, 'rtol': args.rtol,
        'resources': telemetry(workers=1, threads=2), 'summary': summary,
        'runtime_ab': runtime,
        'rows': rows,
        'semantics': 'linear-z second-order Magnus map versus current midpoint map and Cartesian DOP853',
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
