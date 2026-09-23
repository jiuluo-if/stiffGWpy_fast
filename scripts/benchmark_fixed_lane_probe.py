"""Screen an explicitly unrolled four-lane propagation shape.

This is a feasibility probe only.  It does not replace the production kernel:
the four independent scalar states are written as separate variables so LLVM
can expose any genuine SLP/vector opportunity without the grouped-array
overhead rejected in the earlier SoA experiment.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

from numba import njit, types  # noqa: E402


@njit(cache=False, fastmath=False)
def _lane4_kernel(z0, z1, z2, z3, h, steps):
    """Advance four independent phase lanes with identical scalar algebra."""
    x0 = 0.0
    y0 = 1.0
    x1 = 0.0
    y1 = 1.0
    x2 = 0.0
    y2 = 1.0
    x3 = 0.0
    y3 = 1.0
    for _ in range(steps):
        w0 = math.exp(z0)
        q0 = w0 * w0
        om0 = math.sqrt(q0 - 1.0)
        c0 = math.cos(om0 * h)
        s0 = math.sin(om0 * h) / om0
        x0, y0 = (c0 - s0) * x0 - w0 * s0 * y0, w0 * s0 * x0 + (c0 + s0) * y0

        w1 = math.exp(z1)
        q1 = w1 * w1
        om1 = math.sqrt(q1 - 1.0)
        c1 = math.cos(om1 * h)
        s1 = math.sin(om1 * h) / om1
        x1, y1 = (c1 - s1) * x1 - w1 * s1 * y1, w1 * s1 * x1 + (c1 + s1) * y1

        w2 = math.exp(z2)
        q2 = w2 * w2
        om2 = math.sqrt(q2 - 1.0)
        c2 = math.cos(om2 * h)
        s2 = math.sin(om2 * h) / om2
        x2, y2 = (c2 - s2) * x2 - w2 * s2 * y2, w2 * s2 * x2 + (c2 + s2) * y2

        w3 = math.exp(z3)
        q3 = w3 * w3
        om3 = math.sqrt(q3 - 1.0)
        c3 = math.cos(om3 * h)
        s3 = math.sin(om3 * h) / om3
        x3, y3 = (c3 - s3) * x3 - w3 * s3 * y3, w3 * s3 * x3 + (c3 + s3) * y3
    return x0 + y0, x1 + y1, x2 + y2, x3 + y3


def _lane4_output():
    return tuple(float(v) for v in _lane4_kernel(0.25, 0.35, 0.45, 0.55, 0.01, 32))


def _llvm_text():
    _lane4_output()
    signature = (types.float64, types.float64, types.float64, types.float64,
                 types.float64, types.int64)
    return _lane4_kernel.inspect_llvm(signature)


def _llvm_has_vector_lane():
    return bool(re.search(r"<\d+ x double>", _llvm_text()))


def _llvm_vector_tokens():
    return sorted(set(re.findall(r"<\d+ x double>", _llvm_text())))


def _git_head():
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True
    ).strip()


def run_probe():
    values = _lane4_output()
    llvm = _llvm_text()
    return {
        "experiment": "round41_explicit_four_lane_llvm_probe",
        "commit": _git_head(),
        "fastmath": False,
        "lane_count": 4,
        "steps": 32,
        "outputs": list(values),
        "deterministic": values == _lane4_output(),
        "llvm_vector_lane": bool(re.search(r"<\d+ x double>", llvm)),
        "llvm_vector_tokens": sorted(set(re.findall(r"<\d+ x double>", llvm))),
        "llvm_length": len(llvm),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--json",
        default="docs/fixed_lane_probe_round41_20260923.json",
    )
    args = parser.parse_args()
    row = run_probe()
    Path(args.json).write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(row, indent=2))


if __name__ == "__main__":
    main()
