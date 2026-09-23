"""Standalone two-native-interval commutator Magnus screen.

Each pair is replaced by one closed-form exponential of Omega_1 + Omega_2
for the original 2x2 Cartesian system.  This is distinct from the rejected
fourth-order two-exponential CF4 and from residual-guarded pair maps: it has
one exponential and no per-pair oracle/guard.  Production is untouched.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

import numpy as np  # noqa: E402
from numba import njit  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


@njit(cache=False, inline="always")
def _matrix_exp_traceless(a, b, c, d, h):
    """Apply exp(h[[a,b],[c,-a]]) to the state (d is ignored)."""
    _ = d
    aa = h * a
    bb = h * b
    cc = h * c
    delta2 = aa * aa + bb * cc
    if abs(delta2) < 1.0e-28:
        fac = 1.0 + delta2 / 6.0
        scalar = 1.0 + delta2 / 2.0
    elif delta2 > 0.0:
        delta = math.sqrt(delta2)
        scalar = math.cosh(delta)
        fac = math.sinh(delta) / delta
    else:
        delta = math.sqrt(-delta2)
        scalar = math.cos(delta)
        fac = math.sin(delta) / delta
    return scalar, fac, bb, cc


@njit(cache=False)
def _pair_apply(x, y, z0, z2, h):
    w0 = math.exp(z0)
    wm = math.exp(0.5 * (z0 + z2))
    w2 = math.exp(z2)
    a0, b0, c0 = -1.0, -w0, w0
    am, bm, cm = -1.0, -wm, wm
    a1, b1, c1 = -1.0, -w2, w2
    da, db, dc = a1 - a0, b1 - b0, c1 - c0
    # [A_mid, Delta] for traceless 2x2 matrices.
    ca = bm * dc - db * cm
    cb = 2.0 * (am * db - bm * da)
    cc = 2.0 * (cm * da - am * dc)
    H = 2.0 * h
    oa = H * am - (H * H / 12.0) * ca
    ob = H * bm - (H * H / 12.0) * cb
    oc = H * cm - (H * H / 12.0) * cc
    od = -oa
    scalar, fac, ob, oc = _matrix_exp_traceless(oa, ob, oc, od, 1.0)
    return scalar * x + fac * (oa * x + ob * y), scalar * y + fac * (oc * x + od * y)


@njit(cache=False)
def _midpoint_chain(starts, ends, steps):
    x, y = 0.0, 1.0
    for k in range(len(steps)):
        x, y = FS.scaled_step(x, y, 0.5 * (starts[k] + ends[k]), steps[k])
    return x, y


def _path(case, mode=20):
    model, common = _prepared(case, 2)
    _, phi, _, _, _, j0s, z0s, *_ = common
    j0 = int(j0s[mode])
    z0 = float(z0s[mode])
    starts, ends, steps = [], [], []
    for k in range(j0, len(model.Nv) - 1):
        left = z0 + float(phi[k]) - float(phi[j0])
        right = z0 + float(phi[k + 1]) - float(phi[j0])
        if left >= 5.0:
            break
        starts.append(left)
        ends.append(right)
        steps.append(float(model.Nv[k + 1] - model.Nv[k]))
    if len(starts) % 2:
        starts.pop()
        ends.pop()
        steps.pop()
    return np.asarray(starts), np.asarray(ends), np.asarray(steps)


@njit(cache=False)
def _candidate(starts, ends, steps):
    x, y = 0.0, 1.0
    for k in range(0, len(steps), 2):
        x, y = _pair_apply(x, y, starts[k], ends[k + 1], steps[k])
    return x, y


def _reference(starts, ends, steps):
    state = np.array([0.0, 1.0])
    for k in range(0, len(steps), 2):
        z_left = float(starts[k])
        z_right = float(ends[k + 1])
        width = float(steps[k] + steps[k + 1])

        def rhs(n_value, local_state):
            z = z_left + (z_right - z_left) * n_value / width
            w = math.exp(z)
            return np.array([-local_state[0] - w * local_state[1],
                             w * local_state[0] + local_state[1]])

        result = solve_ivp(rhs, (0.0, width), state, method="DOP853",
                           rtol=2e-11, atol=2e-13)
        if not result.success:
            raise RuntimeError(result.message)
        state = result.y[:, -1]
    return state


def _probe(case):
    starts, ends, steps = _path(case)
    # Use the production local constant-midpoint map as a fast baseline for
    # operation count, and an independent DOP853 chain for correctness.
    reference = _reference(starts, ends, steps)
    _midpoint_chain(starts, ends, steps)
    _candidate(starts, ends, steps)
    _midpoint_chain(starts, ends, steps)
    start = time.perf_counter()
    _midpoint_chain(starts, ends, steps)
    baseline_s = time.perf_counter() - start
    start = time.perf_counter()
    candidate = _candidate(starts, ends, steps)
    candidate_s = time.perf_counter() - start
    candidate = np.asarray(candidate)
    ref_power = float(reference @ reference)
    cand_power = float(candidate @ candidate)
    return {
        "case": case,
        "native_intervals": int(len(steps)),
        "candidate_pair_maps": int(len(steps) // 2),
        "finite": bool(np.isfinite(candidate).all()),
        "candidate_power_relative_error": abs(cand_power - ref_power) / max(abs(ref_power), 1e-300),
        "candidate_component_relative_error": float(np.max(
            np.abs(candidate - reference) / np.maximum(np.abs(reference), 1e-300))),
        "midpoint_baseline_seconds": baseline_s,
        "candidate_seconds": candidate_s,
        "candidate_over_midpoint": candidate_s / baseline_s,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = {
        "candidate": "pair_commutator_magnus_round63",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=ROOT, text=True).strip(),
        "production_unchanged": True,
        "settings": vars(args),
        "cases": [_probe(case) for case in (args.case or list(CASES))],
    }
    (ROOT / args.output).write_text(json.dumps(payload, indent=2) + "\n",
                                    encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
