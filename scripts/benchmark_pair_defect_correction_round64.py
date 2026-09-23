"""Standalone pairwise Duhamel defect-correction screen.

The base map is the production constant-midpoint trigonometric transfer over
two native intervals.  A first-order interaction-picture correction from the
linear coefficient defect is then applied algebraically.  This is distinct
from Round 63's exponential of a commutator-corrected generator: it retains
the production map and does not evaluate a second matrix exponential.
Production remains unchanged.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_pair_magnus_round63 import (  # noqa: E402
    CASES,
    _path,
    _reference,
)
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


@njit(cache=False, inline="always")
def _corrected_pair(x, y, z0, z2, h):
    H = 2.0 * h
    zm = 0.5 * (z0 + z2)
    xb, yb = FS.scaled_step(x, y, zm, H)
    w0 = math.exp(z0)
    wm = math.exp(zm)
    w2 = math.exp(z2)
    da = 0.0
    db = -(w2 - w0)
    dc = w2 - w0
    # [A_mid, A_end-A_start] for A=[[-1,-w],[w,1]].
    ca = wm * dc - db * wm
    cb = 2.0 * (-1.0 * db - (-wm) * da)
    cc = 2.0 * (wm * da - (-1.0) * dc)
    factor = -(H * H / 12.0)
    k11 = factor * ca
    k12 = factor * cb
    k21 = factor * cc
    k22 = -k11
    return xb + k11 * xb + k12 * yb, yb + k21 * xb + k22 * yb


@njit(cache=False)
def _corrected_chain(starts, ends, steps):
    x, y = 0.0, 1.0
    for k in range(0, len(steps), 2):
        x, y = _corrected_pair(x, y, starts[k], ends[k + 1], steps[k])
    return x, y


@njit(cache=False)
def _midpoint_chain(starts, ends, steps):
    x, y = 0.0, 1.0
    for k in range(len(steps)):
        x, y = FS.scaled_step(x, y, 0.5 * (starts[k] + ends[k]), steps[k])
    return x, y


def _probe(case):
    starts, ends, steps = _path(case)
    reference = _reference(starts, ends, steps)
    _midpoint_chain(starts, ends, steps)
    _corrected_chain(starts, ends, steps)
    start = time.perf_counter()
    _midpoint_chain(starts, ends, steps)
    baseline_s = time.perf_counter() - start
    start = time.perf_counter()
    candidate = np.asarray(_corrected_chain(starts, ends, steps))
    candidate_s = time.perf_counter() - start
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
        "candidate": "pair_duhamel_defect_correction_round64",
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
