"""Screen warmup crossover for the fused Reference RHS on stratified batches."""

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
Z_TAIL = 5.0
Z_START = -12.0
RTOL = 1e-9
BATCH_SIZES = (32, 64, 128, 256)
PAIR_COUNT = 2
BASELINE_NPZ = ROOT / "docs" / "sample2_reference_root_res8_tail5_gamma0_20261105.npz"
OUTPUT = ROOT / "docs" / "reference_numba_rhs_batch_crossover_20261105.json"


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run(use_fused, python_rhs, model, freqs):
    R._tensor_orig = python_rhs
    R._NUMBA_TENSOR_RHS_READY = use_fused
    started = time.perf_counter()
    values = R.spectrum_reference(
        model, freqs, DN_EFF, z_tail=Z_TAIL, z_start=Z_START,
        rtol=RTOL, workers=1, backend="cartesian",
    )
    return values, time.perf_counter() - started


def main() -> None:
    archive = np.load(BASELINE_NPZ)
    all_freqs = np.asarray(archive["freqs"], dtype=float)
    archive_log = np.asarray(archive["reference_log10OmegaGW"], dtype=float)
    archive_tail = np.asarray(archive["reference_used_tail"], dtype=bool)
    model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
    bg = R._Background(model, DN_EFF)
    old_rhs = R._tensor_orig
    original_threshold = R._NUMBA_TENSOR_RHS_MIN_BATCH_SIZE
    R._NUMBA_TENSOR_RHS_MIN_BATCH_SIZE = max(BATCH_SIZES) + 1

    R._warmup_numba_reference_background(bg)
    compile_started = time.perf_counter()
    R._warmup_numba_reference_tensor_rhs(bg)
    compile_s = time.perf_counter() - compile_started

    results = []
    for size in BATCH_SIZES:
        idx = np.unique(np.linspace(0, all_freqs.size - 1, size).round().astype(int))
        freqs = all_freqs[idx]
        if freqs.size != size:
            raise AssertionError("stratified batch produced duplicate frequency indices")
        times = {"current": [], "fused": []}
        paired = []
        for pair in range(PAIR_COUNT):
            order = ("current", "fused") if pair % 2 == 0 else ("fused", "current")
            run_values = {}
            for method in order:
                run_values[method], elapsed = _run(
                    method == "fused", old_rhs, model, freqs)
                times[method].append(elapsed)
            for left, right in zip(run_values["current"], run_values["fused"]):
                if not np.array_equal(left, right):
                    raise AssertionError(f"paired spectrum mismatch at size={size}")
            paired.append({
                "pair": pair,
                "bitwise_same_Ogw": np.array_equal(
                    run_values["current"][0], run_values["fused"][0]),
                "bitwise_same_Oj": np.array_equal(
                    run_values["current"][1], run_values["fused"][1]),
                "bitwise_same_Opgw": np.array_equal(
                    run_values["current"][2], run_values["fused"][2]),
                "bitwise_same_tail_mask": np.array_equal(
                    run_values["current"][3], run_values["fused"][3]),
            })

        # Independent archived Cartesian values check that this is the same grid.
        archived_log_subset = archive_log[idx]
        current_omega = run_values["current"][0] - run_values["current"][1]
        if not np.array_equal(np.log10(current_omega), archived_log_subset):
            raise AssertionError(f"archived Reference mismatch at size={size}")
        if not np.array_equal(run_values["current"][3], archive_tail[idx]):
            raise AssertionError(f"archived tail mismatch at size={size}")

        current_p50 = statistics.median(times["current"])
        fused_p50 = statistics.median(times["fused"])
        results.append({
            "frequency_count": size,
            "frequency_index_sha256": hashlib.sha256(
                np.asarray(idx, dtype="<i8").tobytes()
            ).hexdigest(),
            "current_runtime_s": times["current"],
            "fused_runtime_s": times["fused"],
            "current_p50_s": current_p50,
            "fused_p50_s": fused_p50,
            "fused_over_current_p50": fused_p50 / current_p50,
            "cold_fused_first_batch_estimate_s": fused_p50 + compile_s,
            "cold_fused_over_current_p50": (fused_p50 + compile_s) / current_p50,
            "paired_checks": paired,
        })

    R._tensor_orig = old_rhs
    R._NUMBA_TENSOR_RHS_READY = False
    R._NUMBA_TENSOR_RHS_MIN_BATCH_SIZE = original_threshold
    result = {
        "schema_version": 1,
        "experiment": "reference_fused_numba_rhs_batch_crossover",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "reference_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "script_sha256": _sha(Path(__file__)),
        "baseline_archive_sha256": _sha(BASELINE_NPZ),
        "parameters": PARAMS,
        "DN_eff": DN_EFF,
        "frequency_grid_sha256": hashlib.sha256(
            np.asarray(all_freqs, dtype="<f8").tobytes()
        ).hexdigest(),
        "batch_sizes": BATCH_SIZES,
        "z_tail": Z_TAIL,
        "z_start": Z_START,
        "rtol": RTOL,
        "atol": [1e-12, 1e-22, 1e-22],
        "workers": 1,
        "jit_compile_and_warm_s": compile_s,
        "paired_runs_per_size": PAIR_COUNT,
        "results": results,
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "workers": 1,
            "concurrent_processes": 1,
        },
        "limitations": [
            "One cosmology and stratified subsets of one archived 800-node grid.",
            "Cold estimate adds the measured one-time JIT warmup to warm P50; it is not a separate process-start measurement.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
