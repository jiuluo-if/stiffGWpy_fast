"""Local feasibility screen for an algebraic characteristic-state basis."""
from __future__ import annotations

import argparse
import json
import math
import statistics
import time
from pathlib import Path

import numpy as np
from numba import njit


@njit(cache=False)
def _cartesian_step(xh, yh, z_mid, h):
    w = math.exp(z_mid)
    w2 = w * w
    if w2 >= 1.0:
        omega = math.sqrt(w2 - 1.0)
        c = math.cos(omega * h)
        si = math.sin(omega * h) / omega
    else:
        omega = math.sqrt(1.0 - w2)
        x = omega * h
        c = 1.0 + 0.5 * x * x
        si = h * (1.0 + x * x / 6.0)
    return (c - si) * xh - w * si * yh, w * si * xh + (c + si) * yh


@njit(cache=False)
def _basis_step(u, v, z_mid, h):
    w = math.exp(z_mid)
    w2 = w * w
    if w2 >= 1.0:
        omega = math.sqrt(w2 - 1.0)
        c = math.cos(omega * h)
        si = math.sin(omega * h) / omega
    else:
        omega = math.sqrt(1.0 - w2)
        x = omega * h
        c = 1.0 + 0.5 * x * x
        si = h * (1.0 + x * x / 6.0)
    a = si * (w - 1.0)
    b = si * (w + 1.0)
    return c * u - a * v, c * v + b * u


@njit(cache=False)
def _cartesian_chain(z0, z1, h, steps):
    x, y = 0.1, 1.0
    dz = (z1 - z0) / steps
    for k in range(steps):
        x, y = _cartesian_step(x, y, z0 + (k + 0.5) * dz, h)
    return x, y


@njit(cache=False)
def _basis_chain(z0, z1, h, steps):
    u, v = 1.1, 0.9
    dz = (z1 - z0) / steps
    for k in range(steps):
        u, v = _basis_step(u, v, z0 + (k + 0.5) * dz, h)
    return u, v


def probe(z0, z1, h, steps, repeats=25):
    base = _cartesian_chain(z0, z1, h, steps)
    cand = _basis_chain(z0, z1, h, steps)
    base_power = base[0] * base[0] + base[1] * base[1]
    cand_power = 0.25 * ((cand[0] - cand[1]) ** 2 +
                         (cand[0] + cand[1]) ** 2)
    bt, ct = [], []
    for _ in range(repeats):
        start = time.perf_counter()
        _cartesian_chain(z0, z1, h, steps)
        bt.append(time.perf_counter() - start)
        start = time.perf_counter()
        _basis_chain(z0, z1, h, steps)
        ct.append(time.perf_counter() - start)
    bm, cm = statistics.median(bt), statistics.median(ct)
    return {
        'z0': z0, 'z1': z1, 'h': h, 'steps': steps, 'repeats': repeats,
        'finite': bool(np.isfinite(cand).all()),
        'power_relative_error': abs(cand_power - base_power) /
        max(abs(base_power), 1e-300),
        'baseline_median_us': bm * 1e6,
        'candidate_median_us': cm * 1e6,
        'candidate_over_baseline': cm / bm,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default='docs/state_basis_transfer_round48_20260923.json')
    args = parser.parse_args()
    rows = [
        probe(0.2, 3.0, 0.01, 256),
        probe(0.2, 5.0, 0.005, 1024),
        probe(-1.2, 0.4, 0.01, 256),
    ]
    payload = {
        'experiment': 'round48_state_basis_transfer_screen',
        'production_path_changed': False,
        'commit': __import__('subprocess').check_output(
            ['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'records': rows,
    }
    Path(args.out).write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(payload, indent=2))


if __name__ == '__main__':
    main()
