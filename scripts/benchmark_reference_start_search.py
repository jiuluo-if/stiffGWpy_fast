"""Paired benchmark for Reference horizon-start search reuse.

The old path calls ``f_hor_abs`` at every bisection step. The candidate caches
the normalization anchors and stops after the floating-point bracket stalls.
The script checks all archived sample2 frequencies, then compares complete
single-mode Cartesian outputs on a fixed representative subset.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import statistics
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "r": 0.19649188544593907,
    "cr": 0.0,
    "n_t": 1.0622501865308032,
    "T_re": 1761.314039059193,
    "kappa10": 0.011598483151630803,
    "DN_re": 10.0,
}
DN_EFF = 2.291411867022825
Z_TAIL = 5.0
Z_START = -12.0
RTOL = 1e-9
REFERENCE_NPZ = ROOT / "docs" / "sample2_reference_root_res8_tail5_gamma0_20261105.npz"
OUTPUT = ROOT / "docs" / "reference_start_search_benchmark_20261105.json"


def _legacy_find_start_N(bg, freq, z_start):
    lo, hi = 0.0, bg.N_inf

    def z_at(N):
        return (freq - R.f_hor_abs(bg, N)) * math.log(10.0)

    if z_at(lo) >= z_start:
        return 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if z_at(mid) < z_start:
            lo = mid
        else:
            hi = mid
    return hi


def _source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_text(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _h2_call_count(bg, freq, fn):
    original = R._H2_and_sigma
    count = 0

    def counted(*args, **kwargs):
        nonlocal count
        count += 1
        return original(*args, **kwargs)

    R._H2_and_sigma = counted
    try:
        value = fn(bg, freq, Z_START)
    finally:
        R._H2_and_sigma = original
    return value, count


def _same_mode_output(left, right):
    keys = (
        "Ogw_today", "Oj_today", "Opgw_today", "used_tail", "event_N",
        "n_steps", "z0",
    )
    return all(left[key] == right[key] for key in keys)


def main() -> None:
    archive = np.load(REFERENCE_NPZ)
    freqs = np.asarray(archive["freqs"], dtype=float)
    freq_hash = hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest()
    model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
    bg = R._Background(model, DN_EFF)
    new_search = R._find_start_N

    # Full archived grid: exact start-position equivalence and aggregate time.
    starts = {"legacy": [], "candidate": []}
    start_ms = {"legacy": [], "candidate": []}
    for pair in range(2):
        order = ("legacy", "candidate") if pair == 0 else ("candidate", "legacy")
        for method in order:
            fn = _legacy_find_start_N if method == "legacy" else new_search
            t0 = time.perf_counter()
            values = [fn(bg, float(freq), Z_START) for freq in freqs]
            start_ms[method].append((time.perf_counter() - t0) * 1000.0)
            if pair == 0:
                starts[method] = values
    if not all(a == b for a, b in zip(starts["legacy"], starts["candidate"])):
        raise AssertionError("candidate start positions differ on the archived grid")

    representative_indices = sorted({0, 200, 400, 600, len(freqs) - 1})
    h2_counts = {}
    for index in representative_indices:
        freq = float(freqs[index])
        old_value, old_calls = _h2_call_count(bg, freq, _legacy_find_start_N)
        new_value, new_calls = _h2_call_count(bg, freq, new_search)
        if old_value != new_value:
            raise AssertionError(f"start position changed at index {index}")
        h2_counts[str(index)] = {
            "freq": freq,
            "legacy_h2_evaluations": old_calls,
            "candidate_h2_evaluations": new_calls,
        }

    R._warmup_numba_reference_background(bg)
    mode_rows = []
    for index in representative_indices:
        freq = float(freqs[index])
        times = {"legacy": [], "candidate": []}
        outputs = {"legacy": None, "candidate": None}
        for pair in range(4):
            order = ("legacy", "candidate") if pair % 2 == 0 else (
                "candidate", "legacy"
            )
            for method in order:
                R._find_start_N = (
                    _legacy_find_start_N if method == "legacy" else new_search
                )
                t0 = time.perf_counter()
                outputs[method] = R.solve_reference_mode(
                    model,
                    freq,
                    DN_EFF,
                    z_tail=Z_TAIL,
                    z_start=Z_START,
                    rtol=RTOL,
                    backend="cartesian",
                )
                times[method].append(time.perf_counter() - t0)
        R._find_start_N = new_search
        if not _same_mode_output(outputs["legacy"], outputs["candidate"]):
            raise AssertionError(f"mode output changed at index {index}")
        mode_rows.append({
            "index": index,
            "freq": freq,
            "bitwise_output_match": True,
            "nfev": outputs["candidate"]["n_steps"],
            "legacy_p50_s": statistics.median(times["legacy"]),
            "legacy_p95_s": float(np.percentile(times["legacy"], 95)),
            "candidate_p50_s": statistics.median(times["candidate"]),
            "candidate_p95_s": float(np.percentile(times["candidate"], 95)),
            "paired_p50_speedup": statistics.median(times["legacy"])
            / statistics.median(times["candidate"]),
        })

    R._find_start_N = new_search
    try:
        import numba
        import scipy
        scipy_version = scipy.__version__
        numba_version = numba.__version__
    except ImportError:
        scipy_version = numba_version = None
    result = {
        "schema_version": 1,
        "experiment": "reference_start_search_anchor_reuse",
        "created_at": datetime.now().astimezone().isoformat(),
        "baseline_commit": _git_text("rev-parse", "HEAD"),
        "branch": _git_text("branch", "--show-current"),
        "worktree_dirty": bool(_git_text("status", "--porcelain")),
        "baseline_reference_sha256": hashlib.sha256(
            subprocess.check_output(
                ["git", "show", "HEAD:stiffgwpy_fast/reference.py"], cwd=ROOT
            )
        ).hexdigest(),
        "candidate_reference_sha256": _source_hash(ROOT / "stiffgwpy_fast" / "reference.py"),
        "parameters": PARAMS,
        "DN_eff": DN_EFF,
        "frequency_count": int(freqs.size),
        "frequency_sha256": freq_hash,
        "z_tail": Z_TAIL,
        "z_start": Z_START,
        "rtol": RTOL,
        "backend": "cartesian",
        "warm_jit_before_full_mode_pairs": True,
        "full_grid_start_positions_bitwise_equal": True,
        "full_grid_legacy_ms": start_ms["legacy"],
        "full_grid_candidate_ms": start_ms["candidate"],
        "full_grid_start_p50_speedup": statistics.median(start_ms["legacy"])
        / statistics.median(start_ms["candidate"]),
        "representative_h2_evaluations": h2_counts,
        "representative_mode_pairs": mode_rows,
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy_version,
            "numba": numba_version,
            "workers": 1,
        },
        "limitations": [
            "Full-grid benchmark compares start searches only; complete ODE outputs are paired on five representative modes.",
            "The frequency array is the archived sample2 res8 grid, not a global parameter-space certification.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
