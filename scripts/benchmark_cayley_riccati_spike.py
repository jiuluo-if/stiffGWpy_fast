"""Tiny pole-free projective Riccati/Cayley feasibility screen.

The raw ``x/y`` Riccati coordinate hit physical poles in every tested mode.
This screen uses the Cayley coordinate ``q=(x+i y)/(x-i y)`` plus a log
amplitude, which has no real denominator zero.  It is intentionally a scalar
method screen, not a production-kernel replacement.
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

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

from numba import njit  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


@njit(cache=False)
def _cartesian_loop(z0, z1, h, steps):
    x_value = 0.1
    y_value = 1.0
    z_step = (z1 - z0) / steps
    for step in range(steps):
        z_mid = z0 + (step + 0.5) * z_step
        x_value, y_value = FS.scaled_step(x_value, y_value, z_mid, h)
    return x_value, y_value


@njit(cache=False)
def _cayley_loop(z0, z1, h, steps):
    x_initial = 0.1
    y_initial = 1.0
    initial_norm = math.hypot(x_initial, y_initial)
    q_value = complex(x_initial, y_initial) / complex(x_initial, -y_initial)
    log_amplitude = math.log(initial_norm)
    z_step = (z1 - z0) / steps
    for step in range(steps):
        z_mid = z0 + (step + 0.5) * z_step
        theta = 0.5 * math.atan2(q_value.imag, q_value.real)
        x_unit = math.cos(theta)
        y_unit = math.sin(theta)
        x_next, y_next = FS.scaled_step(x_unit, y_unit, z_mid, h)
        norm = math.hypot(x_next, y_next)
        log_amplitude += math.log(norm)
        q_value = complex(x_next, y_next) / complex(x_next, -y_next)
    theta = 0.5 * math.atan2(q_value.imag, q_value.real)
    amplitude = math.exp(log_amplitude)
    return amplitude * math.cos(theta), amplitude * math.sin(theta)


def _cayley_probe(z0, z1, h, steps, repeats=20):
    cartesian = _cartesian_loop(z0, z1, h, steps)
    cayley = _cayley_loop(z0, z1, h, steps)
    cartesian_power = cartesian[0] ** 2 + cartesian[1] ** 2
    cayley_power = cayley[0] ** 2 + cayley[1] ** 2
    cartesian_times = []
    cayley_times = []
    for _ in range(int(repeats)):
        start = time.perf_counter()
        _cartesian_loop(z0, z1, h, steps)
        cartesian_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _cayley_loop(z0, z1, h, steps)
        cayley_times.append(time.perf_counter() - start)
    cartesian_median = statistics.median(cartesian_times)
    cayley_median = statistics.median(cayley_times)
    return {
        "z0": float(z0),
        "z1": float(z1),
        "h": float(h),
        "steps": int(steps),
        "finite": bool(math.isfinite(cayley[0]) and math.isfinite(cayley[1])),
        "power_relative_error": abs(cayley_power - cartesian_power)
        / max(abs(cartesian_power), 1e-300),
        "cartesian_median_us": cartesian_median * 1e6,
        "cayley_median_us": cayley_median * 1e6,
        "cayley_over_cartesian": cayley_median / cartesian_median,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="docs/cayley_riccati_round45_20260923.json")
    args = parser.parse_args()
    probes = [
        _cayley_probe(0.2, 3.0, 0.01, 256),
        _cayley_probe(0.2, 5.0, 0.005, 1024),
        _cayley_probe(-0.2, 4.0, 0.005, 1024),
    ]
    payload = {
        "experiment": "round45_cayley_regularized_projective_riccati_screen",
        "production_path_changed": False,
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "records": probes,
    }
    Path(args.out).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
