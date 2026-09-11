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
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _prufer_derivatives(z, sigma, phase):
    """Return derivatives for x=A sin(phase), y=A cos(phase)."""
    return (1.5 * sigma - 1.0,
            1.5 * sigma - 2.0 + math.cos(2.0 * phase),
            -math.exp(z) - math.sin(2.0 * phase))


def _phase_averaged_power(amplitude):
    return 0.5 * amplitude * amplitude


def _phase_averaged_today(coefficient_squared, z_handoff, event_N, n_inf,
                          today_f_hor, freq, tensor_power):
    transfer_squared = coefficient_squared \
        * math.exp(-2.0 * z_handoff + 2.0 * event_N - 2.0 * n_inf)
    x_squared = transfer_squared * 10.0 ** (2.0 * (freq - today_f_hor))
    oj = -transfer_squared * tensor_power / 3.0
    opgw = x_squared * tensor_power / 36.0
    return {
        'Ogw_today': float(3.0 * opgw + oj),
        'Oj_today': float(oj),
        'Opgw_today': float(opgw),
    }


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
    tensor_power = model.derived_param['A_t'] \
        * (10.0 ** freq / gp.f_piv) ** model.derived_param['nt']
    if not result.t_events[0].size:
        zf, log_amplitude, phase = result.y[:, -1]
        amplitude = math.exp(float(log_amplitude))
        x = amplitude * math.sin(float(phase))
        y = amplitude * math.cos(float(phase))
        th = y / math.exp(float(zf))
        return {
            'Ogw_today': float((x * x + y * y) / 24.0 * tensor_power
                               + x * th / 3.0 * tensor_power),
            'Oj_today': float(x * th / 3.0 * tensor_power),
            'Opgw_today': float((-5.0 * x * x + 7.0 * y * y)
                                / 72.0 * tensor_power),
            'used_tail': False,
            'nfev': int(result.nfev),
        }
    event_N = float(result.t_events[0][0])
    zf, log_amplitude, phase = result.y_events[0][0]
    amplitude = math.exp(float(log_amplitude))
    f_today = REF.f_hor_abs(bg, freq, bg.N_inf)
    today = _phase_averaged_today(
        _phase_averaged_power(amplitude), float(zf), event_N, bg.N_inf,
        f_today, freq, tensor_power)
    return dict(today, **{
        'used_tail': True,
        'event_N': event_N,
        'z_handoff': float(zf),
        'amplitude_handoff': amplitude / math.sqrt(2.0),
        'phase_handoff': float(phase),
        'phase_averaged_power': _phase_averaged_power(amplitude),
        'nfev': int(result.nfev),
    })


def spectrum_prufer(model, freqs, dn_eff, z_tail, rtol=1e-10):
    freqs = np.asarray(freqs, dtype=float)
    rows = [solve_prufer(model, float(freq), dn_eff, z_tail, rtol=rtol)
            for freq in freqs]
    ogw = np.asarray([row['Ogw_today'] for row in rows])
    oj = np.asarray([row['Oj_today'] for row in rows])
    opgw = np.asarray([row['Opgw_today'] for row in rows])
    g2, quadrature_error, interpolation_error = REF.integrate_spectrum(
        freqs, ogw, oj)
    omega_nu = gp.Omega_nh2 / model.derived_param['h'] ** 2
    return {
        'Ogw': ogw,
        'Oj': oj,
        'Opgw': opgw,
        'DN_gw': float(gp.Neff0 * g2 / omega_nu),
        'quadrature_error': float(quadrature_error),
        'interpolation_error': interpolation_error,
        'used_tail_fraction': float(np.mean([row['used_tail'] for row in rows])),
        'rows': rows,
    }


def compare_spectrum(model, freqs, dn_eff, z_tail, rtol=1e-10):
    prufer = spectrum_prufer(model, freqs, dn_eff, z_tail, rtol=rtol)
    ogw, oj, opgw, used_tail = REF.spectrum_reference(
        model, freqs, dn_eff, z_tail=z_tail, rtol=rtol, workers=1)
    omega_nu = gp.Omega_nh2 / model.derived_param['h'] ** 2
    g2, _, _ = REF.integrate_spectrum(freqs, ogw, oj)
    cartesian_dn = float(gp.Neff0 * g2 / omega_nu)
    return {
        'z_tail': float(z_tail),
        'prufer_DN_gw': prufer['DN_gw'],
        'cartesian_DN_gw': cartesian_dn,
        'DN_relative_error': abs(prufer['DN_gw'] - cartesian_dn)
        / max(abs(cartesian_dn), 1e-300),
        'Ogw_max_relative_error': float(np.max(
            np.abs(prufer['Ogw'] - ogw) / np.maximum(np.abs(ogw), 1e-300))),
        'Oj_max_relative_error': float(np.max(
            np.abs(prufer['Oj'] - oj) / np.maximum(np.abs(oj), 1e-300))),
        'Opgw_max_relative_error': float(np.max(
            np.abs(prufer['Opgw'] - opgw) / np.maximum(np.abs(opgw), 1e-300))),
        'used_tail_fraction_prufer': prufer['used_tail_fraction'],
        'used_tail_fraction_cartesian': float(np.mean(used_tail)),
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
    opgw_rel = abs(prufer['Opgw_today'] - cartesian['Opgw_today']) \
        / max(abs(cartesian['Opgw_today']), 1e-300)
    oj_rel = abs(prufer['Oj_today'] - cartesian['Oj_today']) \
        / max(abs(cartesian['Oj_today']), 1e-300)
    ogw_rel = abs(prufer['Ogw_today'] - cartesian['Ogw_today']) \
        / max(abs(cartesian['Ogw_today']), 1e-300)
    return {
        'frequency': float(freq), 'z_tail': float(z_tail), 'used_tail': True,
        'amplitude_relative_error': float(amplitude_rel),
        'phase_delta': float(phase_delta),
        'phase_averaged_power_relative_error': float(power_rel),
        'Opgw_relative_error': float(opgw_rel),
        'Oj_relative_error': float(oj_rel),
        'Ogw_relative_error': float(ogw_rel),
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
    spectrum = [compare_spectrum(model, freqs, dn_eff, z_tail)
                for z_tail in args.z_tail]
    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'point': args.point,
        'frequencies': freqs.tolist(),
        'z_tails': [float(value) for value in args.z_tail],
        'DN_eff': dn_eff,
        'resources': telemetry(workers=1, threads=2),
        'rows': rows,
        'spectrum': spectrum,
        'semantics': 'standalone Prüfer amplitude-phase versus Cartesian DOP853',
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
