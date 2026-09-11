"""Compare a standalone Prüfer amplitude-phase oracle with Cartesian DOP853."""

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
from scipy.integrate import solve_ivp  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benchmark_same_grid_reference import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _prufer_derivatives(z, sigma, phase):
    """Return derivatives for x=A sin(phase), y=A cos(phase)."""
    return (1.5 * sigma - 1.0,
            1.5 * sigma - 2.0 + math.cos(2.0 * phase),
            -math.exp(z) - math.sin(2.0 * phase))


def _phase_averaged_power(amplitude):
    return 0.5 * amplitude * amplitude


def _prufer_rhs(N, state, bg):
    z, log_amplitude, phase = state
    _, sigma = REF._H2_and_sigma(bg, N)
    dz, dlog_amplitude, dphase = _prufer_derivatives(z, sigma, phase)
    return dz, dlog_amplitude, dphase


def solve_prufer(model, freq, dn_eff, z_tail, rtol=1e-10):
    bg = REF._Background(model, dn_eff)
    n_start = REF._find_start_N(bg, freq, -12.0)
    z0 = (freq - REF.f_hor_abs(bg, freq, n_start)) * REF.ln10
    result = solve_ivp(
        _prufer_rhs, (n_start, bg.N_inf), (z0, z0, 0.0),
        method='DOP853', rtol=rtol, atol=[1e-12, 1e-12, 1e-12],
        events=[REF._make_tail_event(z_tail)], args=(bg,))
    if not result.t_events[0].size:
        return {'used_tail': False, 'nfev': int(result.nfev)}
    event_N = float(result.t_events[0][0])
    zf, log_amplitude, phase = result.y_events[0][0]
    amplitude = math.exp(float(log_amplitude))
    return {
        'used_tail': True,
        'event_N': event_N,
        'z_handoff': float(zf),
        'amplitude_handoff': amplitude / math.sqrt(2.0),
        'phase_handoff': float(phase),
        'phase_averaged_power': _phase_averaged_power(amplitude),
        'nfev': int(result.nfev),
    }


def compare_mode(model, freq, dn_eff, z_tail, rtol=1e-10):
    start = time.perf_counter()
    cartesian = REF.solve_reference_mode(model, freq, dn_eff,
                                         z_tail=z_tail, rtol=rtol)
    cartesian_s = time.perf_counter() - start
    start = time.perf_counter()
    prufer = solve_prufer(model, freq, dn_eff, z_tail, rtol=rtol)
    prufer_s = time.perf_counter() - start
    if not cartesian['used_tail'] or not prufer['used_tail']:
        return {'frequency': float(freq), 'z_tail': float(z_tail),
                'used_tail': False}
    phase_delta = (prufer['phase_handoff'] - cartesian['phase_handoff'] + math.pi) \
        % (2.0 * math.pi) - math.pi
    amplitude_rel = abs(prufer['amplitude_handoff']
                        - cartesian['amplitude_handoff']) \
        / max(abs(cartesian['amplitude_handoff']), 1e-300)
    reference_power = cartesian['amplitude_handoff'] ** 2
    power_rel = abs(prufer['phase_averaged_power'] - reference_power) \
        / max(reference_power, 1e-300)
    return {
        'frequency': float(freq), 'z_tail': float(z_tail), 'used_tail': True,
        'amplitude_relative_error': float(amplitude_rel),
        'phase_delta': float(phase_delta),
        'phase_averaged_power_relative_error': float(power_rel),
        'cartesian_nfev': int(cartesian['n_steps']),
        'prufer_nfev': int(prufer['nfev']),
        'cartesian_runtime_s': cartesian_s,
        'prufer_runtime_s': prufer_s,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--point', choices=sorted(CASES), default='default')
    parser.add_argument('--freq-count', type=int, default=8)
    parser.add_argument('--z-tail', nargs='+', type=float, default=[5.0, 7.0])
    parser.add_argument('--out', default='docs/oracle_prufer_default.json')
    args = parser.parse_args(argv)
    if args.freq_count < 3:
        raise ValueError('freq-count must be at least 3')

    FS.apply_accuracy_mode('fast')
    model = LCDM_SG(**CASES[args.point])
    FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                      frequency_quadrature='simpson')
    native = np.sort(np.asarray(model.f, dtype=float))
    picks = np.linspace(0, native.size - 1, args.freq_count, dtype=int)
    freqs = native[picks]
    dn_eff = float(model.cosmo_param['DN_eff'])
    rows = [compare_mode(model, float(freq), dn_eff, z_tail)
            for freq in freqs for z_tail in args.z_tail]
    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'point': args.point,
        'frequencies': freqs.tolist(),
        'z_tails': [float(value) for value in args.z_tail],
        'DN_eff': dn_eff,
        'resources': telemetry(workers=1, threads=2),
        'rows': rows,
        'semantics': 'standalone Prüfer amplitude-phase versus Cartesian DOP853',
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
