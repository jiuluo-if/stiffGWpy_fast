"""Standalone phase-increment recurrence screen for the Cartesian transfer."""

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

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment
apply_environment()
from numba import njit  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.benchmark_phase_recurrence import _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


@njit(cache=False)
def _cart(z0, z1, h, steps):
    x, y = 0.1, 1.0
    dz = (z1 - z0) / steps
    for k in range(steps):
        x, y = FS.scaled_step(x, y, z0 + (k + 0.5) * dz, h)
    return x, y


@njit(cache=False)
def _recurrence(z0, z1, h, steps, reanchor):
    x, y = 0.1, 1.0
    dz = (z1 - z0) / steps
    c, s, previous = 1.0, 0.0, 0.0
    for k in range(steps):
        z = z0 + (k + 0.5) * dz
        w = math.exp(z)
        omega = math.sqrt(w * w - 1.0)
        angle = omega * h
        if k % reanchor == 0:
            c, s = math.cos(angle), math.sin(angle)
        else:
            delta = angle - previous
            cd = 1.0 - 0.5 * delta * delta
            sd = delta - delta * delta * delta / 6.0
            c, s = c * cd - s * sd, s * cd + c * sd
        previous = angle
        si = s / omega
        x, y = (c - si) * x - w * si * y, w * si * x + (c + si) * y
    return x, y


@njit(cache=False)
def _path(z_values, h_values, x, y, reanchor):
    angle_previous = 0.0
    c, s = 1.0, 0.0
    for k in range(len(z_values)):
        w = math.exp(z_values[k])
        if w * w < 1.0:
            x, y = FS.scaled_step(x, y, z_values[k], h_values[k])
            angle_previous = 0.0
            c, s = 1.0, 0.0
            continue
        omega = math.sqrt(w * w - 1.0)
        angle = omega * h_values[k]
        if k % reanchor == 0:
            c, s = math.cos(angle), math.sin(angle)
        else:
            delta = angle - angle_previous
            cd = 1.0 - 0.5 * delta * delta
            sd = delta - delta * delta * delta / 6.0
            c, s = c * cd - s * sd, s * cd + c * sd
        angle_previous = angle
        si = s / omega
        x, y = (c - si) * x - w * si * y, w * si * x + (c + si) * y
    return x, y


def _actual_path_probe(case="default", mode=20, reanchor=32):
    model, common = _prepared(case, 2)
    Nv, phi, phi_mid, _, s2inv, j0s, z0s, *_ = common
    j0 = int(j0s[mode])
    z0 = float(z0s[mode])
    phi0 = float(phi[j0])
    z_values = []
    h_values = []
    for k in range(j0, len(Nv) - 1):
        z_start = z0 + float(phi[k]) - phi0
        if z_start >= 5.0:
            break
        z_end = z0 + float(phi[k + 1]) - phi0
        zmid = z0 + float(phi_mid[k]) - phi0
        n = int(FS._phase_substeps(0.005, zmid, 0.25))
        hs = 0.005 / n
        for sub in range(n):
            z_values.append(zmid + (z_end - z_start) * (2 * sub + 1) / (2 * n))
            h_values.append(hs)
    zv = np.asarray(z_values, dtype=np.float64)
    hv = np.asarray(h_values, dtype=np.float64)
    x0 = 0.0
    y0 = math.exp(z0) * float(s2inv[j0])
    base = _path(zv, hv, x0, y0, 1)  # reanchor every step is exact sin/cos baseline
    cand = _path(zv, hv, x0, y0, reanchor)
    bp = base[0] ** 2 + base[1] ** 2
    cp = cand[0] ** 2 + cand[1] ** 2
    return {
        "case": case,
        "mode": mode,
        "substeps": len(zv),
        "reanchor": reanchor,
        "power_relative_error": abs(cp - bp) / max(abs(bp), 1e-300),
        "finite": bool(np.isfinite(cand).all()),
    }


def _probe(z0, z1, h, steps, reanchor, repeats=20):
    base = _cart(z0, z1, h, steps)
    cand = _recurrence(z0, z1, h, steps, reanchor)
    bp = base[0] * base[0] + base[1] * base[1]
    cp = cand[0] * cand[0] + cand[1] * cand[1]
    bt = []
    ct = []
    for _ in range(repeats):
        t = time.perf_counter()
        _cart(z0, z1, h, steps)
        bt.append(time.perf_counter() - t)
        t = time.perf_counter()
        _recurrence(z0, z1, h, steps, reanchor)
        ct.append(time.perf_counter() - t)
    bm = statistics.median(bt)
    cm = statistics.median(ct)
    return {
        "z0": z0,
        "z1": z1,
        "h": h,
        "steps": steps,
        "reanchor": reanchor,
        "finite": math.isfinite(cand[0]) and math.isfinite(cand[1]),
        "power_relative_error": abs(cp - bp) / max(abs(bp), 1e-300),
        "cartesian_us": bm * 1e6,
        "candidate_us": cm * 1e6,
        "candidate_over_baseline": cm / bm,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="docs/phase_increment_recurrence_round46_20260923.json")
    a = p.parse_args()
    rows = [_probe(0.2, 3, 0.01, 256, 32), _probe(0.2, 5, 0.005, 1024, 32), _probe(0.2, 5, 0.005, 1024, 64)]
    for case in ("default", "lowT", "highT", "stiff", "high_kappa"):
        for mode in (0, 20, 40, 60):
            rows.append(_actual_path_probe(case, mode, 32))
    payload = {
        "experiment": "round46_phase_increment_recurrence_screen",
        "production_path_changed": False,
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "records": rows,
    }
    Path(a.out).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
