"""Standalone combined transfer-map lookup-table feasibility screen.

The table stores the already-composed constant-z transfer coefficients for the
canonical ``h=0.005`` phase substeps.  It is not a production path: this first
screen only measures cubic interpolation error and table footprint before a
full kernel twin is justified.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

import numpy as np  # noqa: E402
from numba import njit  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402

Z_MIN = -1.0
Z_MAX = 5.1
TABLE_NODES = 4097
N_SUB_MAX = 32


def _build_tables():
    z = np.linspace(Z_MIN, Z_MAX, TABLE_NODES, dtype=np.float64)
    w = np.exp(z)
    c_table = np.empty((N_SUB_MAX + 1, TABLE_NODES), dtype=np.float64)
    si_table = np.empty_like(c_table)
    for n_sub in range(1, N_SUB_MAX + 1):
        h = 0.005 / n_sub
        w2 = w * w
        c = np.empty_like(w)
        si = np.empty_like(w)
        high = w2 > 1.0 + 1.0e-14
        low = w2 < 1.0 - 1.0e-14
        omega = np.sqrt(np.maximum(w2 - 1.0, 0.0))
        c[high] = np.cos(omega[high] * h)
        si[high] = np.sin(omega[high] * h) / omega[high]
        omega_low = np.sqrt(np.maximum(1.0 - w2, 0.0))
        x = omega_low[low] * h
        c[low] = 1.0 + 0.5 * x * x
        si[low] = h * (1.0 + x * x / 6.0)
        c[~(high | low)] = 1.0
        si[~(high | low)] = h
        c_table[n_sub] = c
        si_table[n_sub] = si
    return z, w, c_table, si_table


@njit(inline="always", cache=True)
def _cubic(y0, y1, y2, y3, t):
    t2 = t * t
    t3 = t2 * t
    return 0.5 * (2.0 * y1 + (-y0 + y2) * t
                  + (2.0 * y0 - 5.0 * y1 + 4.0 * y2 - y3) * t2
                  + (-y0 + 3.0 * y1 - 3.0 * y2 + y3) * t3)


@njit(inline="always", cache=True)
def _lookup(table, row, u):
    index = int(u)
    if index < 1:
        index = 1
    if index > TABLE_NODES - 3:
        index = TABLE_NODES - 3
    t = u - index
    return _cubic(table[row, index - 1], table[row, index],
                  table[row, index + 1], table[row, index + 2], t)


@njit(inline="always", cache=True)
def _table_step(xh, yh, z_mid, h, n_sub, z_grid, w_table, c_table, si_table):
    u = (z_mid - Z_MIN) * (TABLE_NODES - 1) / (Z_MAX - Z_MIN)
    w = _lookup(w_table.reshape(1, w_table.size), 0, u)
    c = _lookup(c_table, n_sub, u)
    si = _lookup(si_table, n_sub, u)
    return (c - si) * xh - w * si * yh, w * si * xh + (c + si) * yh


def _local_table_gate():
    z_grid, w_table, c_table, si_table = _build_tables()
    del z_grid
    worst = 0.0
    for n_sub in range(1, N_SUB_MAX + 1):
        h = 0.005 / n_sub
        for z_mid in np.linspace(Z_MIN + 1e-6, Z_MAX - 1e-6, 601):
            actual = _table_step(
                0.2, -0.1, float(z_mid), h, n_sub,
                np.array([Z_MIN, Z_MAX]), w_table, c_table, si_table)
            expected = FS.scaled_step(0.2, -0.1, float(z_mid), h)
            worst = max(worst, abs(actual[0] - expected[0]),
                        abs(actual[1] - expected[1]))
    return {
        "max_transfer_absolute_error": float(worst),
        "table_bytes": int(w_table.nbytes + c_table.nbytes + si_table.nbytes),
        "table_nodes": TABLE_NODES,
        "n_sub_max": N_SUB_MAX,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = {
        "candidate": "combined_transfer_map_cubic_table_feasibility",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "settings": vars(args),
        "metrics": _local_table_gate(),
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
