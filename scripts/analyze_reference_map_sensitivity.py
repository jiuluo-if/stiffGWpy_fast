"""Estimate the res4 fixed-map slope from two archived full Reference maps.

This is a post-processing secant estimate, not a new solver evaluation or a
rigorous enclosure of the fixed point.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reference_source_sha(commit):
    data = subprocess.check_output(
        ["git", "show", "%s:stiffgwpy_fast/reference.py" % commit],
        cwd=ROOT)
    return hashlib.sha256(data).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map-low", type=Path,
                        default=ROOT / "docs" /
                        "res4_sigma_exact_reference_map_20261104.json")
    parser.add_argument("--grid-low", type=Path,
                        default=ROOT / "docs" /
                        "res4_sigma_exact_reference_map_20261104.npz")
    parser.add_argument("--map-high", type=Path,
                        default=ROOT / "docs" /
                        "res4_reference_secant_root_tail10_20261104.json")
    parser.add_argument("--grid-high", type=Path,
                        default=ROOT / "docs" /
                        "res4_reference_secant_root_tail10_20261104.npz")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs" /
                        "res4_reference_map_sensitivity_20261104.json")
    args = parser.parse_args(argv)

    low = json.loads(args.map_low.read_text(encoding="utf-8"))
    high = json.loads(args.map_high.read_text(encoding="utf-8"))
    low_commit = str(low["source_sha"])
    high_commit = str(high["source_commit"])
    low_module_sha = _reference_source_sha(low_commit)
    high_module_sha = _reference_source_sha(high_commit)
    if low_module_sha != high_module_sha:
        raise RuntimeError("archived Reference maps came from different solver sources")

    with np.load(args.grid_low) as archive:
        low_freqs = np.asarray(archive["freqs"], dtype=float)
    with np.load(args.grid_high) as archive:
        high_freqs = np.asarray(archive["freqs"], dtype=float)
    if low_freqs.shape != high_freqs.shape or not np.array_equal(
            low_freqs, high_freqs):
        raise RuntimeError("archived Reference maps do not share the same frequency grid")
    if len(low_freqs) != int(low["frequency_count"]) \
            or len(high_freqs) != int(high["frequency_count"]):
        raise RuntimeError("frequency archive length does not match map metadata")

    x1 = float(low["dn_eff"])
    f1 = float(low["reference_dn_map"])
    x2 = float(high["input_dn_eff"])
    f2 = float(high["reference_dn_map"])
    if x1 == x2:
        raise ValueError("map inputs must be distinct for a secant estimate")
    slope = (f2 - f1) / (x2 - x1)
    residual = f2 - x2
    correction = residual / (1.0 - slope)
    if not all(math.isfinite(v) for v in (slope, residual, correction)):
        raise FloatingPointError("non-finite secant sensitivity result")

    payload = {
        "schema_version": 1,
        "analysis": "res4_cartesian_reference_fixed_map_secant_sensitivity",
        "source_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "analysis_script_sha256": _sha256(Path(__file__).resolve()),
        "reference_source_sha256_both_maps": low_module_sha,
        "reference_source_commits": [low_commit, high_commit],
        "map_inputs_dn_eff": [x1, x2],
        "map_outputs_dn_gw": [f1, f2],
        "map_residuals": [f1 - x1, residual],
        "secant_slope_dDNgw_dDNeff": slope,
        "root_shift_estimate_at_second_map": correction,
        "secant_corrected_root_estimate": x2 + correction,
        "frequency_count": int(len(low_freqs)),
        "frequency_arrays_bitwise_equal": True,
        "low_archive_npz_sha256": _sha256(args.grid_low),
        "high_archive_npz_sha256": _sha256(args.grid_high),
        "maps": [
            {"json": str(args.map_low.resolve()),
             "json_sha256": _sha256(args.map_low),
             "npz": str(args.grid_low.resolve()),
             "npz_sha256": _sha256(args.grid_low),
             "source_commit": low_commit},
            {"json": str(args.map_high.resolve()),
             "json_sha256": _sha256(args.map_high),
             "npz": str(args.grid_high.resolve()),
             "npz_sha256": _sha256(args.grid_high),
             "source_commit": high_commit},
        ],
        "protocol_note": (
            "Both records identify the res4 tail10 Reference map; the first JSON "
            "records rtol=1e-11 explicitly, while the second is the archived "
            "tail10 secant-root map. The slope is a secant average over the "
            "two stored DN_eff inputs, not a local derivative or a rigorous root bracket."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
