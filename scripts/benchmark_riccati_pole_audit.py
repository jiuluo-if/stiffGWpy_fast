"""Standalone audit of the raw tensor ratio Riccati reduction.

This is a throwaway research probe.  It never changes the production solver.
The ratio r=x/y has poles when y crosses zero, so a pole is reported as an
explicit fail-closed diagnostic rather than being continued or regularized.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_same_grid_reference import CASES as REFERENCE_CASES  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = dict(REFERENCE_CASES)
CASES.update({
    'edge_r_lo': dict(r=1.26e-6, cr=1, T_re=2e3, kappa10=1e-2),
    'edge_r_hi': dict(r=7.94e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'edge_tre_hi': dict(r=1e-2, cr=1, T_re=7.94e5, kappa10=1e-2),
    'sobol_000': dict(r=3.4845505440304265, cr=1.0,
                      n_t=-0.12252988945692778,
                      T_re=34989.571710260374,
                      DN_re=18.787082182243466,
                      kappa10=5.28549585835966e-06),
})


def _rhs(n_value, state, bg):
    z, ratio, log_abs_y = state
    _, sigma = REF._H2_and_sigma(bg, n_value)
    z_prime = 1.5 * sigma - 1.0
    omega = math.exp(z)
    # Exact quotient reduction of REF._tensor_orig, valid only while y != 0.
    ratio_prime = -2.0 * ratio - omega * (1.0 + ratio * ratio)
    log_abs_y_prime = -1.0 + 1.5 * sigma + omega * ratio
    return z_prime, ratio_prime, log_abs_y_prime


def _tail_event(z_tail):
    def event(n_value, state, bg):
        del n_value, bg
        return state[0] - z_tail

    event.terminal = True
    event.direction = 1
    return event


def _pole_event(sign, ratio_limit):
    def event(n_value, state, bg):
        del n_value, bg
        return sign * state[1] - ratio_limit

    event.terminal = True
    event.direction = 1
    return event


def audit_mode(model, freq, dn_eff, z_tail, ratio_limit, rtol):
    bg = REF._Background(model, dn_eff)
    n_start = REF._find_start_N(bg, freq, -12.0)
    z0 = (freq - REF.f_hor_abs(bg, freq, n_start)) * REF.ln10
    result = solve_ivp(
        _rhs,
        (n_start, bg.N_inf),
        (z0, 0.0, z0),
        method='DOP853',
        rtol=rtol,
        atol=[1e-12, 1e-12, 1e-12],
        events=[_tail_event(z_tail), _pole_event(1.0, ratio_limit),
                _pole_event(-1.0, ratio_limit)],
        args=(bg,),
    )
    tail = bool(result.t_events[0].size)
    pole_events = [event for event in result.t_events[1:] if event.size]
    pole = bool(pole_events)
    if pole:
        pole_n = float(next(event[0] for event in result.t_events[1:] if event.size))
        terminal = 'pole'
    elif tail:
        pole_n = None
        terminal = 'tail'
    else:
        pole_n = None
        terminal = 'domain_end'
    cartesian = None
    if pole_n is not None:
        cart_result = solve_ivp(
            REF._tensor_orig,
            (n_start, pole_n),
            (z0, 0.0, math.exp(z0)),
            method='DOP853',
            rtol=rtol,
            atol=[1e-12, 1e-22, 1e-22],
            args=(bg,),
        )
        z_cart, x_cart, y_cart = cart_result.y[:, -1]
        cartesian = {
            'success': bool(cart_result.success),
            'nfev': int(cart_result.nfev),
            'z': float(z_cart),
            'x': float(x_cart),
            'y': float(y_cart),
            'abs_y_over_abs_x': float(abs(y_cart) / max(abs(x_cart), 1e-300)),
        }
    return {
        'freq': float(freq),
        'terminal': terminal,
        'pole': pole,
        'tail_reached': tail,
        'pole_N': pole_n,
        'nfev': int(result.nfev),
        'success': bool(result.success),
        'message': str(result.message),
        'cartesian_at_pole': cartesian,
    }


def run_case(name, kwargs, frequencies, args):
    model = LCDM_SG(**kwargs)
    rows = [audit_mode(model, freq, args.dn_eff, args.z_tail,
                       args.ratio_limit, args.rtol)
            for freq in frequencies]
    poles = sum(row['pole'] for row in rows)
    return {
        'case': name,
        'parameters': kwargs,
        'frequencies': [float(value) for value in frequencies],
        'rows': rows,
        'pole_count': int(poles),
        'tail_count': int(sum(row['tail_reached'] for row in rows)),
        'numerical_failure_count': int(sum(not row['success'] for row in rows)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cases', nargs='+', default=['default', 'highT', 'stiff',
                                                        'high_kappa', 'edge_r_lo',
                                                        'edge_r_hi', 'edge_tre_hi',
                                                        'sobol_000'])
    parser.add_argument('--freqs', nargs='+', type=float,
                        default=np.linspace(-4.0, 2.0, 8).tolist())
    parser.add_argument('--dn-eff', type=float, default=0.0)
    parser.add_argument('--z-tail', type=float, default=5.0)
    parser.add_argument('--ratio-limit', type=float, default=1.0e6)
    parser.add_argument('--rtol', type=float, default=1.0e-10)
    parser.add_argument('--out', default='docs/riccati_pole_audit_round_20260917.json')
    args = parser.parse_args()

    unknown = sorted(set(args.cases) - set(CASES))
    if unknown:
        raise ValueError('unknown cases: %s' % ', '.join(unknown))
    records = [run_case(name, CASES[name], args.freqs, args) for name in args.cases]
    payload = {
        'experiment': 'raw_riccati_ratio_pole_audit',
        'production_path_changed': False,
        'guard_policy': 'pole is terminal and reported; no fallback or restart',
        'commit': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'resource': telemetry(workers=1, threads=2, concurrent_processes=1),
        'settings': vars(args),
        'records': records,
    }
    output = ROOT / args.out
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n',
                      encoding='utf-8')
    print(json.dumps({
        'commit': payload['commit'],
        'cases': len(records),
        'modes': sum(len(record['rows']) for record in records),
        'poles': sum(record['pole_count'] for record in records),
        'tails': sum(record['tail_count'] for record in records),
        'numerical_failures': sum(record['numerical_failure_count'] for record in records),
        'out': str(output),
    }, sort_keys=True))


if __name__ == '__main__':
    main()
