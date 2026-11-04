"""Paired full-grid Reference A/B for start-search anchor reuse.

This measures the full Cartesian fixed-background map and requires bitwise
identity in the returned 800-node spectrum and ``DN_gw`` value.
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
FREQ_RES = 8.0
Z_TAIL = 5.0
Z_START = -12.0
RTOL = 1e-9
PAIR_COUNT = 3
BASELINE_NPZ = ROOT / "docs" / "sample2_reference_root_res8_tail5_gamma0_20261105.npz"
OUTPUT = ROOT / "docs" / "reference_start_search_fullgrid_ab_20261105.json"


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


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args):
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _run_map(method):
    R._find_start_N = _legacy_find_start_N if method == "legacy" else _CANDIDATE
    model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
    started = time.perf_counter()
    result = R.run_reference(
        model,
        dn_eff=DN_EFF,
        freq_res=FREQ_RES,
        z_tail=Z_TAIL,
        z_start=Z_START,
        rtol=RTOL,
        self_consistent=False,
        workers=1,
        physical_guard=False,
        backend="cartesian",
    )
    return result, time.perf_counter() - started


def main():
    global _CANDIDATE
    _CANDIDATE = R._find_start_N
    old = np.load(BASELINE_NPZ)
    old_freqs = np.asarray(old["freqs"], dtype=float)
    old_spectrum = np.asarray(old["reference_log10OmegaGW"], dtype=float)
    old_tail = np.asarray(old["reference_used_tail"], dtype=bool)
    if old_freqs.size != 800 or old_spectrum.size != 800:
        raise RuntimeError("archived sample2 Reference output is not 800 nodes")

    warm_started = time.perf_counter()
    R._warmup_numba_reference_background(
        R._Background(LCDM_SG(**PARAMS, DN_eff=DN_EFF), DN_EFF)
    )
    warmup_s = time.perf_counter() - warm_started
    runs = {"legacy": [], "candidate": []}
    outputs = {"legacy": [], "candidate": []}
    for pair in range(PAIR_COUNT):
        order = ("legacy", "candidate") if pair % 2 == 0 else (
            "candidate", "legacy"
        )
        for method in order:
            result, elapsed = _run_map(method)
            runs[method].append(elapsed)
            outputs[method].append(result)
            freqs = np.asarray(result["freqs"], dtype=float)
            spectrum = np.asarray(result["log10OmegaGW"], dtype=float)
            tail = np.asarray(result["used_tail"], dtype=bool)
            if freqs.size != 800 or spectrum.size != 800:
                raise RuntimeError(f"{method} map returned an unexpected grid")
            if not np.array_equal(freqs, old_freqs):
                raise AssertionError(f"{method} frequency grid differs from archive")
            if not np.array_equal(spectrum, old_spectrum):
                raise AssertionError(f"{method} spectrum differs from archive")
            if not np.array_equal(tail, old_tail):
                raise AssertionError(f"{method} tail mask differs from archive")

    # Pair every legacy/candidate run and assert bitwise map and spectrum equality.
    pair_checks = []
    for i in range(PAIR_COUNT):
        left, right = outputs["legacy"][i], outputs["candidate"][i]
        left_spectrum = np.asarray(left["log10OmegaGW"], dtype=float)
        right_spectrum = np.asarray(right["log10OmegaGW"], dtype=float)
        left_freqs = np.asarray(left["freqs"], dtype=float)
        right_freqs = np.asarray(right["freqs"], dtype=float)
        left_tail = np.asarray(left["used_tail"], dtype=bool)
        right_tail = np.asarray(right["used_tail"], dtype=bool)
        same_dn = float(left["DN_gw"]) == float(right["DN_gw"])
        same_spectrum = np.array_equal(left_spectrum, right_spectrum)
        same_grid = np.array_equal(left_freqs, right_freqs)
        same_tail = np.array_equal(left_tail, right_tail)
        if not (same_dn and same_spectrum and same_grid and same_tail):
            raise AssertionError(f"full-grid A/B mismatch in pair {i}")
        pair_checks.append({
            "pair": i,
            "bitwise_same_DN_gw": same_dn,
            "bitwise_same_spectrum": same_spectrum,
            "bitwise_same_grid": same_grid,
            "bitwise_same_tail_mask": same_tail,
            "DN_gw": float(right["DN_gw"]),
        })

    result = {
        "schema_version": 1,
        "experiment": "reference_start_search_fullgrid_ab",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "worktree_dirty": bool(_git("status", "--porcelain")),
        "baseline_reference_sha256": hashlib.sha256(
            subprocess.check_output(
                ["git", "show", "HEAD:stiffgwpy_fast/reference.py"], cwd=ROOT
            )
        ).hexdigest(),
        "candidate_reference_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "baseline_npz_sha256": _sha(BASELINE_NPZ),
        "parameters": PARAMS,
        "DN_eff": DN_EFF,
        "frequency_res": FREQ_RES,
        "frequency_count": 800,
        "frequency_sha256": hashlib.sha256(
            np.asarray(old_freqs, dtype="<f8").tobytes()
        ).hexdigest(),
        "z_tail": Z_TAIL,
        "z_start": Z_START,
        "rtol": RTOL,
        "backend": "cartesian",
        "workers": 1,
        "jit_warmup_s_excluded": warmup_s,
        "paired_runs": PAIR_COUNT,
        "legacy_runtime_s": runs["legacy"],
        "candidate_runtime_s": runs["candidate"],
        "legacy_p50_s": statistics.median(runs["legacy"]),
        "legacy_p95_s": float(np.percentile(runs["legacy"], 95)),
        "candidate_p50_s": statistics.median(runs["candidate"]),
        "candidate_p95_s": float(np.percentile(runs["candidate"], 95)),
        "candidate_over_legacy_p50": statistics.median(runs["candidate"])
        / statistics.median(runs["legacy"]),
        "paired_checks": pair_checks,
        "archive_max_abs_spectrum_delta_dex": float(np.max(np.abs(
            np.asarray(outputs["candidate"][0]["log10OmegaGW"], dtype=float)
            - old_spectrum
        ))),
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "workers": 1,
            "concurrent_processes": 1,
        },
        "limitations": [
            "This is one fixed cosmology, DN_eff background and frequency grid, not a global speed claim.",
            "The test holds the Cartesian Reference propagation and all solver tolerances fixed.",
        ],
    }
    R._find_start_N = _CANDIDATE
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
