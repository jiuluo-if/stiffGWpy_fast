"""Prototype an independent finite-window phase-averaged tail observable.

The production fast solver is not changed.  A mode is integrated to an initial
handoff, continued through a finite oscillatory window, and only then reduced
to a phase-averaged amplitude before the analytic tail.  The output is evidence
for an Oracle-B candidate, not a certified replacement for Oracle A.
"""

import argparse
import json
import math
import os
import sys

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

WINDOWS = (0.0, 0.5, 1.0, 2.0)


def _phase_window_observable(z_handoff, event_N, x_handoff, y_handoff,
                             n_inf, today_f_hor, freq, tensor_power):
    """Return the phase-averaged today ``Opgw`` from a handoff state."""
    coefficient_sq = 0.5 * (x_handoff * x_handoff + y_handoff * y_handoff)
    transfer_sq = math.exp(-2.0 * z_handoff + 2.0 * event_N - 2.0 * n_inf)
    x_today_sq = coefficient_sq * transfer_sq \
        * 10.0 ** (2.0 * (freq - today_f_hor))
    return {
        'observable': float(x_today_sq * tensor_power / 36.0),
        'phase_averaged': True,
        'coefficient_squared': float(coefficient_sq),
    }


def _start_state(model, freq, dn_eff, z_start):
    bg = REF._Background(model, dn_eff)
    n_start = REF._find_start_N(bg, freq, -12.0)
    f_start = REF.f_hor_abs(bg, freq, n_start)
    z0 = (freq - f_start) * REF.ln10
    y0 = math.exp(z0)
    result = solve_ivp(
        REF._tensor_orig, (n_start, bg.N_inf), (z0, 0.0, y0),
        method='DOP853', rtol=1e-10, atol=[1e-12, 1e-22, 1e-22],
        events=[REF._make_tail_event(z_start)], args=(bg,), dense_output=True)
    if not result.t_events[0].size:
        return bg, None
    event_N = float(result.t_events[0][0])
    return bg, (event_N, *map(float, result.sol(event_N)))


def solve_window(model, freq, dn_eff, z_start=5.0, window=0.0):
    """Solve one finite-window phase-averaged diagnostic mode."""
    bg, initial = _start_state(model, freq, dn_eff, z_start)
    if initial is None:
        return {'used_tail': False, 'window': float(window)}
    start_N, z0, x0, y0 = initial
    target_z = z_start + float(window)
    if window > 0.0:
        continued = solve_ivp(
            REF._tensor_orig, (start_N, bg.N_inf), (z0, x0, y0),
            method='DOP853', rtol=1e-10, atol=[1e-12, 1e-22, 1e-22],
            events=[REF._make_tail_event(target_z)], args=(bg,),
            dense_output=True)
        if not continued.t_events[0].size:
            return {'used_tail': False, 'window': float(window)}
        handoff_N = float(continued.t_events[0][0])
        zf, xf, yf = map(float, continued.sol(handoff_N))
    else:
        handoff_N, zf, xf, yf = start_N, z0, x0, y0
    p_t = model.derived_param['A_t'] \
        * (10.0 ** freq / REF.gp.f_piv) ** model.derived_param['nt']
    f_today = REF.f_hor_abs(bg, freq, bg.N_inf)
    result = _phase_window_observable(
        zf, handoff_N, xf, yf, bg.N_inf, f_today, freq, p_t)
    result.update({
        'used_tail': True,
        'window': float(window),
        'event_N': handoff_N,
        'z_handoff': zf,
        'phase_handoff': math.atan2(xf, yf),
        'amplitude_handoff': math.sqrt(0.5 * (xf * xf + yf * yf)),
    })
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--point', choices=sorted(CASES), default='default')
    parser.add_argument('--freq-count', type=int, default=4)
    parser.add_argument('--out', default='docs/oracle_b_phase_window_head.json')
    args = parser.parse_args(argv)

    FS.apply_accuracy_mode('fast')
    model = LCDM_SG(**CASES[args.point])
    FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                      frequency_quadrature='pchip')
    freqs = np.asarray(model.f, dtype=float)
    picks = np.linspace(0, freqs.size - 1, args.freq_count, dtype=int)
    freqs = freqs[picks]
    dn_eff = float(model.cosmo_param['DN_eff'])
    rows = []
    for freq in freqs:
        for window in WINDOWS:
            row = solve_window(model, float(freq), dn_eff, window=window)
            row.update({'point': args.point, 'frequency': float(freq),
                        'DN_eff': dn_eff})
            rows.append(row)
    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'point': args.point,
        'frequencies': freqs.tolist(),
        'windows': list(WINDOWS),
        'resources': telemetry(workers=1, threads=2),
        'rows': rows,
        'semantics': 'fixed DN_eff; finite DOP853 phase window; phase-averaged tail',
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
