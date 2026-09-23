"""Native-grid Prüfer transfer screen using the production discrete map."""
from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
import time
from pathlib import Path

import numpy as np
from numba import njit


@njit(cache=False)
def _step(xh, yh, z_mid, h):
    w = math.exp(z_mid)
    w2 = w * w
    if w2 >= 1.0:
        omega = math.sqrt(w2 - 1.0)
        c = math.cos(omega * h)
        si = math.sin(omega * h) / omega
    else:
        omega = math.sqrt(1.0 - w2)
        xx = omega * h
        c = 1.0 + 0.5 * xx * xx
        si = h * (1.0 + xx * xx / 6.0)
    return (c - si) * xh - w * si * yh, w * si * xh + (c + si) * yh


@njit(cache=False)
def _cartesian_chain(z0, z1, h, steps):
    x, y = 0.1, 1.0
    dz = (z1 - z0) / steps
    for k in range(steps):
        x, y = _step(x, y, z0 + (k + 0.5) * dz, h)
    return x, y


@njit(cache=False)
def _prufer_chain(z0, z1, h, steps):
    amplitude = math.sqrt(1.01)
    phase = math.atan2(0.1, 1.0)
    dz = (z1 - z0) / steps
    for k in range(steps):
        z = z0 + (k + 0.5) * dz
        x = amplitude * math.sin(phase)
        y = amplitude * math.cos(phase)
        x, y = _step(x, y, z, h)
        amplitude = math.sqrt(x * x + y * y)
        phase = math.atan2(x, y)
    return amplitude, phase


def probe(z0, z1, h, steps, repeats=25):
    cart = _cartesian_chain(z0, z1, h, steps)
    prufer = _prufer_chain(z0, z1, h, steps)
    cart_power = cart[0] * cart[0] + cart[1] * cart[1]
    prufer_power = prufer[0] * prufer[0]
    bt, pt = [], []
    for _ in range(repeats):
        start = time.perf_counter()
        _cartesian_chain(z0, z1, h, steps)
        bt.append(time.perf_counter() - start)
        start = time.perf_counter()
        _prufer_chain(z0, z1, h, steps)
        pt.append(time.perf_counter() - start)
    bm, pm = statistics.median(bt), statistics.median(pt)
    return {
        'z0': z0, 'z1': z1, 'h': h, 'steps': steps, 'repeats': repeats,
        'finite': bool(np.isfinite(prufer).all()),
        'power_relative_error': abs(prufer_power - cart_power) /
        max(abs(cart_power), 1e-300),
        'baseline_median_us': bm * 1e6,
        'candidate_median_us': pm * 1e6,
        'candidate_over_baseline': pm / bm,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default='docs/native_prufer_transfer_round49_20260923.json')
    args = parser.parse_args()
    rows = [
        probe(0.2, 3.0, 0.01, 256),
        probe(0.2, 5.0, 0.005, 1024),
        probe(-1.2, 0.4, 0.01, 256),
    ]
    payload = {
        'experiment': 'round49_native_prufer_transfer_screen',
        'production_path_changed': False,
        'commit': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'records': rows,
    }
    Path(args.out).write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(payload, indent=2))


if __name__ == '__main__':
    main()
