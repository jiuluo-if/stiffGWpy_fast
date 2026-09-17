"""Standalone WKB-carrier spike for the tensor propagation tail.

The probe keeps Cartesian DOP853 propagation up to ``z_match`` and replaces
the remaining oscillatory interval with a scalar carrier.  It is research
only: no production API, default, fallback, or guard semantics are changed.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, telemetry  # noqa: E402
from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

apply_environment()


def _carrier_rhs(n_value, state, bg):
    z, envelope, phase, cosine_integral, sine_integral = state
    _, sigma = REF._H2_and_sigma(bg, n_value)
    return (1.5 * sigma - 1.0, 1.5 * sigma - 2.0, math.exp(z),
            math.cos(2.0 * phase), math.sin(2.0 * phase))


def _event(level):
    def event(n_value, state, bg):
        del n_value, bg
        return state[0] - level

    event.terminal = True
    event.direction = 1
    return event


def _initial_state(model, freq, dn_eff):
    bg = REF._Background(model, dn_eff)
    n_start = REF._find_start_N(bg, freq, -12.0)
    z0 = (freq - REF.f_hor_abs(bg, freq, n_start)) * REF.ln10
    return bg, n_start, z0


def _exact_to_level(model, freq, dn_eff, level, rtol):
    bg, n_start, z0 = _initial_state(model, freq, dn_eff)
    result = solve_ivp(
        REF._tensor_orig,
        (n_start, bg.N_inf),
        (z0, 0.0, math.exp(z0)),
        method='DOP853',
        rtol=rtol,
        atol=[1e-12, 1e-22, 1e-22],
        events=[_event(level)],
        args=(bg,),
    )
    if not result.t_events[0].size:
        raise RuntimeError('Cartesian reference did not reach z_match')
    n_match = float(result.t_events[0][0])
    return bg, n_match, result.y_events[0][0], int(result.nfev)


def _exact_to_tail(model, freq, dn_eff, z_tail, rtol):
    bg, n_start, z0 = _initial_state(model, freq, dn_eff)
    result = solve_ivp(
        REF._tensor_orig,
        (n_start, bg.N_inf),
        (z0, 0.0, math.exp(z0)),
        method='DOP853',
        rtol=rtol,
        atol=[1e-12, 1e-22, 1e-22],
        events=[_event(z_tail)],
        args=(bg,),
    )
    if not result.t_events[0].size:
        raise RuntimeError('Cartesian reference did not reach z_tail')
    return float(result.t_events[0][0]), result.y_events[0][0], int(result.nfev)


def _carrier_to_tail(bg, n_match, state_match, z_tail, rtol):
    z_match, x_match, y_match = [float(value) for value in state_match]
    carrier = solve_ivp(
        _carrier_rhs,
        (n_match, bg.N_inf),
        (z_match, 0.0, 0.0, 0.0, 0.0),
        method='DOP853',
        rtol=rtol,
        atol=[1e-12, 1e-12, 1e-12, 1e-12, 1e-12],
        events=[_event(z_tail)],
        args=(bg,),
    )
    if not carrier.t_events[0].size:
        raise RuntimeError('carrier did not reach z_tail')
    n_tail = float(carrier.t_events[0][0])
    _, envelope, phase, cosine_integral, sine_integral = carrier.y_events[0][0]
    w_match = complex(y_match, x_match)
    carrier_correction = complex(cosine_integral, sine_integral)
    a_tail = w_match + w_match.conjugate() * carrier_correction
    w_tail = math.exp(envelope) * complex(math.cos(-phase), math.sin(-phase)) \
        * a_tail
    return n_tail, float(w_tail.imag), float(w_tail.real), int(carrier.nfev)


def _metrics(candidate, exact):
    x_c, y_c = candidate
    x_e, y_e = exact
    amp_c = math.hypot(x_c, y_c)
    amp_e = math.hypot(x_e, y_e)
    return {
        'amplitude_rel': abs(amp_c - amp_e) / max(amp_e, 1e-300),
        'phase_abs': abs(math.atan2(x_c, y_c) - math.atan2(x_e, y_e)),
        'state_rel_l2': math.hypot(x_c - x_e, y_c - y_e)
                         / max(amp_e, 1e-300),
    }


def run_mode(model, freq, dn_eff, z_match, z_tail, rtol):
    bg, n_match, state_match, nfev_match = _exact_to_level(
        model, freq, dn_eff, z_match, rtol)
    n_tail, exact_state, nfev_exact = _exact_to_tail(
        model, freq, dn_eff, z_tail, rtol)
    n_tail_candidate, x_candidate, y_candidate, nfev_carrier = _carrier_to_tail(
        bg, n_match, state_match, z_tail, rtol)
    return {
        'freq': float(freq),
        'z_match': float(z_match),
        'n_match': n_match,
        'n_tail': n_tail,
        'n_tail_candidate': n_tail_candidate,
        'metrics': _metrics((x_candidate, y_candidate),
                            (float(exact_state[1]), float(exact_state[2]))),
        'nfev_exact': nfev_exact,
        'nfev_match_plus_carrier': nfev_match + nfev_carrier,
        'numerical_failure': False,
    }


def runtime_sample(model_kwargs, frequencies, z_match, z_tail, repeats, rtol):
    model = LCDM_SG(**model_kwargs)
    for freq in frequencies:
        run_mode(model, freq, 0.0, z_match, z_tail, rtol)
    exact_times = []
    candidate_times = []
    candidate_model = LCDM_SG(**model_kwargs)

    def candidate_once():
        for freq in frequencies:
            bg, n_match, state_match, _ = _exact_to_level(
                candidate_model, freq, 0.0, z_match, rtol)
            _carrier_to_tail(bg, n_match, state_match, z_tail, rtol)

    for _ in range(repeats):
        start = time.perf_counter()
        for freq in frequencies:
            _exact_to_tail(model, freq, 0.0, z_tail, rtol)
        exact_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        candidate_once()
        candidate_times.append(time.perf_counter() - start)
    return {
        'exact_median_ms': statistics.median(exact_times) * 1e3,
        'candidate_median_ms': statistics.median(candidate_times) * 1e3,
        'exact_p95_ms': float(np.percentile(exact_times, 95) * 1e3),
        'candidate_p95_ms': float(np.percentile(candidate_times, 95) * 1e3),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--z-match', type=float, default=1.0)
    parser.add_argument('--z-tail', type=float, default=5.0)
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--rtol', type=float, default=1e-9)
    parser.add_argument('--out', default='docs/wkb_carrier_spike_round_20260917.json')
    args = parser.parse_args()
    frequencies = [-4.0, -2.0, 0.0, 1.0, 2.0]
    names = ['default', 'highT', 'stiff', 'high_kappa', 'edge_r_hi', 'sobol_000']
    records = []
    for name in names:
        model = LCDM_SG(**CASES[name])
        rows = [run_mode(model, freq, 0.0, args.z_match, args.z_tail,
                         args.rtol) for freq in frequencies]
        runtime = runtime_sample(CASES[name], frequencies, args.z_match,
                                 args.z_tail, args.repeats, args.rtol)
        records.append({'case': name, 'rows': rows, 'runtime': runtime})
    payload = {
        'experiment': 'wkb_carrier_standalone_spike',
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
    metrics = [row['metrics'] for record in records for row in record['rows']]
    print(json.dumps({
        'commit': payload['commit'],
        'modes': len(metrics),
        'max_amplitude_rel': max(row['amplitude_rel'] for row in metrics),
        'max_phase_abs': max(row['phase_abs'] for row in metrics),
        'max_state_rel_l2': max(row['state_rel_l2'] for row in metrics),
        'out': str(output),
    }, sort_keys=True))


if __name__ == '__main__':
    main()
