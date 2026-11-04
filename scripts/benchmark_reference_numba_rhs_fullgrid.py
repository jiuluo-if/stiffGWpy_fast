"""Paired full-grid benchmark for a fused Numba Reference RHS prototype."""

from __future__ import annotations

import hashlib
import json
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
OUTPUT = ROOT / "docs" / "reference_numba_rhs_fullgrid_ab_20261105.json"


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_map(use_fused, model, python_rhs):
    R._tensor_orig = python_rhs
    R._NUMBA_TENSOR_RHS_READY = use_fused
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


def main() -> None:
    archived = np.load(BASELINE_NPZ)
    old_freqs = np.asarray(archived["freqs"], dtype=float)
    old_spectrum = np.asarray(archived["reference_log10OmegaGW"], dtype=float)
    old_tail = np.asarray(archived["reference_used_tail"], dtype=bool)
    model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
    bg = R._Background(model, DN_EFF)
    old_rhs = R._tensor_orig
    original_rhs_threshold = R._NUMBA_TENSOR_RHS_MIN_BATCH_SIZE
    # Keep the batch auto-warmup from changing the baseline route during A/B.
    R._NUMBA_TENSOR_RHS_MIN_BATCH_SIZE = 801

    compile_started = time.perf_counter()
    R._warmup_numba_reference_background(bg)
    R._warmup_numba_reference_tensor_rhs(bg)
    compile_s = time.perf_counter() - compile_started

    times = {"current_compiled_h2": [], "fused_numba_rhs": []}
    outputs = {"current_compiled_h2": [], "fused_numba_rhs": []}
    for pair in range(PAIR_COUNT):
        order = (("current_compiled_h2", "fused_numba_rhs") if pair % 2 == 0
                 else ("fused_numba_rhs", "current_compiled_h2"))
        for method in order:
            result, elapsed = _run_map(
                method == "fused_numba_rhs", model, old_rhs)
            times[method].append(elapsed)
            outputs[method].append(result)
            freqs = np.asarray(result["freqs"], dtype=float)
            spectrum = np.asarray(result["log10OmegaGW"], dtype=float)
            tail = np.asarray(result["used_tail"], dtype=bool)
            if freqs.size != 800 or not np.array_equal(freqs, old_freqs):
                raise AssertionError(f"{method} changed the archived frequency grid")
            if not np.array_equal(spectrum, old_spectrum):
                raise AssertionError(f"{method} differs from archived Cartesian Reference")
            if not np.array_equal(tail, old_tail):
                raise AssertionError(f"{method} changed the tail mask")

    pair_checks = []
    for i in range(PAIR_COUNT):
        current, fused = (outputs["current_compiled_h2"][i],
                          outputs["fused_numba_rhs"][i])
        exact = (float(current["DN_gw"]) == float(fused["DN_gw"])
                 and np.array_equal(current["log10OmegaGW"],
                                    fused["log10OmegaGW"])
                 and np.array_equal(current["used_tail"], fused["used_tail"]))
        if not exact:
            raise AssertionError(f"candidate outputs differ in pair {i}")
        pair_checks.append({
            "pair": i,
            "bitwise_same_DN_gw": True,
            "bitwise_same_spectrum": True,
            "bitwise_same_tail_mask": True,
            "DN_gw": float(fused["DN_gw"]),
        })

    R._tensor_orig = old_rhs
    R._NUMBA_TENSOR_RHS_READY = False
    R._NUMBA_TENSOR_RHS_MIN_BATCH_SIZE = original_rhs_threshold
    current_p50 = statistics.median(times["current_compiled_h2"])
    candidate_p50 = statistics.median(times["fused_numba_rhs"])
    result = {
        "schema_version": 1,
        "experiment": "reference_fused_numba_rhs_fullgrid_ab",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "reference_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "script_sha256": _sha(Path(__file__)),
        "baseline_archive_sha256": _sha(BASELINE_NPZ),
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
        "jit_compile_and_warm_s_excluded": compile_s,
        "paired_runs": PAIR_COUNT,
        "current_compiled_h2_runtime_s": times["current_compiled_h2"],
        "fused_numba_rhs_runtime_s": times["fused_numba_rhs"],
        "current_compiled_h2_p50_s": current_p50,
        "current_compiled_h2_p95_s": float(
            np.percentile(times["current_compiled_h2"], 95)),
        "fused_numba_rhs_p50_s": candidate_p50,
        "fused_numba_rhs_p95_s": float(
            np.percentile(times["fused_numba_rhs"], 95)),
        "candidate_over_current_p50": candidate_p50 / current_p50,
        "pair_checks": pair_checks,
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "workers": 1,
            "concurrent_processes": 1,
        },
        "limitations": [
            "One cosmology and one fixed 800-frequency background only.",
            "The fused RHS is an experiment-only dispatch in this run; no production default is enabled.",
            "The separate JIT warmup is reported, and this run does not measure cold Reference startup or multi-worker throughput.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
