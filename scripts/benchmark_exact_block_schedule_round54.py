"""Round 54 exact-transfer composition Amdahl screen.

This does not implement a candidate.  It counts the unavoidable operations of
an exact composition schedule: every phase substep still needs its scalar
coefficient functions, while composing general 2x2 maps adds matrix products.
The screen uses real prepared fast paths to prevent a purely symbolic result
from being mistaken for a production opportunity.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _prepared  # noqa: E402


def _operation_model(native_intervals, substeps, block_size):
    maps = int(native_intervals) * int(substeps)
    groups = (int(native_intervals) + int(block_size) - 1) // int(block_size)
    extra_products = 0
    for group in range(groups):
        native_in_group = min(block_size, native_intervals - group * block_size)
        maps_in_group = native_in_group * substeps
        extra_products += max(0, maps_in_group - 1)
    return {
        "native_intervals": int(native_intervals),
        "substeps_per_native": int(substeps),
        "block_size": int(block_size),
        "baseline_transcendentals": maps,
        "candidate_transcendentals": maps,
        "candidate_extra_matrix_products": int(extra_products),
        "candidate_extra_scalar_multiply_lower_bound": int(4 * extra_products),
    }


def _phase_substeps(z0, z1, h=0.005, phase_max=0.25):
    zmid = 0.5 * (z0 + z1)
    if phase_max <= 0.0 or zmid <= 0.0:
        return 1
    return max(1, int(math.ceil(h * math.exp(zmid) / phase_max)))


def _count_case(case_name, mode_limit=16, step_limit=256, block_size=2):
    _, common = _prepared(case_name, 2)
    _, phi_grid, _, _, _, j0s, z0s, *_ = common
    h = float(common[14])
    total_native = 0
    total_maps = 0
    max_substeps = 0
    modes = np.unique(np.linspace(0, len(j0s) - 1, mode_limit, dtype=int))
    for mode in modes:
        j0 = int(j0s[mode])
        phi0 = float(phi_grid[j0])
        native = 0
        k = j0
        while k + 1 < len(phi_grid) and native < step_limit:
            z0 = float(z0s[mode] + phi_grid[k] - phi0)
            z1 = float(z0s[mode] + phi_grid[k + 1] - phi0)
            if z0 >= 5.0:
                break
            count = _phase_substeps(z0, z1, h, 0.25)
            total_native += 1
            total_maps += count
            max_substeps = max(max_substeps, count)
            native += 1
            k += 1
    groups = (total_native + block_size - 1) // block_size
    extra_products = max(0, total_maps - groups)
    return {
        "case": case_name,
        "sampled_modes": int(len(modes)),
        "native_intervals": int(total_native),
        "phase_maps": int(total_maps),
        "max_phase_substeps": int(max_substeps),
        "block_size": int(block_size),
        "baseline_transcendentals": int(total_maps),
        "candidate_transcendentals": int(total_maps),
        "candidate_extra_matrix_products": int(extra_products),
        "candidate_extra_scalar_multiply_lower_bound": int(4 * extra_products),
        "strict_equivalence": True,
        "work_reduction": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = {
        "experiment": "round54_exact_block_schedule_amdahl_screen",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "operation_probe": _operation_model(8, 3, 2),
        "rows": [_count_case(case) for case in CASES],
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
