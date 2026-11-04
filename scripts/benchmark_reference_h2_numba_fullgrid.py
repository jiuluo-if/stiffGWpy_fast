"""Paired full-grid sample2 benchmark for the integrated H2/sigma Numba path."""

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
OUTPUT = ROOT / "docs" / "reference_h2_numba_fullgrid_ab_20261105.json"
_PYTHON_EVAL = R._H2_and_sigma_python
_CANDIDATE_EVAL = R._H2_and_sigma


def _git(*args):
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_map(method):
    R._H2_and_sigma = _PYTHON_EVAL if method == "python" else _CANDIDATE_EVAL
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
    old_archive = np.load(BASELINE_NPZ)
    old_freqs = np.asarray(old_archive["freqs"], dtype=float)
    old_spectrum = np.asarray(old_archive["reference_log10OmegaGW"], dtype=float)
    old_tail = np.asarray(old_archive["reference_used_tail"], dtype=bool)
    model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
    bg = R._Background(model, DN_EFF)
    compile_started = time.perf_counter()
    R._warmup_numba_reference_background(bg)
    compile_s = time.perf_counter() - compile_started

    times = {"python": [], "compiled": []}
    outputs = {"python": [], "compiled": []}
    for pair in range(PAIR_COUNT):
        order = ("python", "compiled") if pair % 2 == 0 else (
            "compiled", "python"
        )
        for method in order:
            result, elapsed = _run_map(method)
            times[method].append(elapsed)
            outputs[method].append(result)
            freqs = np.asarray(result["freqs"], dtype=float)
            spectrum = np.asarray(result["log10OmegaGW"], dtype=float)
            tail = np.asarray(result["used_tail"], dtype=bool)
            if freqs.size != 800 or not np.array_equal(freqs, old_freqs):
                raise AssertionError(f"{method} changed the archived frequency grid")
            if not np.array_equal(spectrum, old_spectrum):
                raise AssertionError(f"{method} spectrum differs from archived Reference")
            if not np.array_equal(tail, old_tail):
                raise AssertionError(f"{method} tail mask differs from archived Reference")

    pair_checks = []
    for i in range(PAIR_COUNT):
        py, nb = outputs["python"][i], outputs["compiled"][i]
        exact_spectrum = np.array_equal(py["log10OmegaGW"], nb["log10OmegaGW"])
        exact_dn = float(py["DN_gw"]) == float(nb["DN_gw"])
        exact_tail = np.array_equal(py["used_tail"], nb["used_tail"])
        if not (exact_spectrum and exact_dn and exact_tail):
            raise AssertionError(f"H2/sigma implementations differ in pair {i}")
        pair_checks.append({
            "pair": i,
            "bitwise_same_DN_gw": exact_dn,
            "bitwise_same_spectrum": exact_spectrum,
            "bitwise_same_tail_mask": exact_tail,
            "DN_gw": float(nb["DN_gw"]),
        })

    R._H2_and_sigma = _PYTHON_EVAL
    result = {
        "schema_version": 1,
        "experiment": "reference_h2_sigma_numba_fullgrid_ab",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "reference_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "candidate_kernel": "Reference._h2_sigma_numba_core via _H2_and_sigma dispatch",
        "candidate_reference_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
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
        "jit_compile_s_excluded": compile_s,
        "paired_runs": PAIR_COUNT,
        "python_runtime_s": times["python"],
        "compiled_runtime_s": times["compiled"],
        "python_p50_s": statistics.median(times["python"]),
        "python_p95_s": float(np.percentile(times["python"], 95)),
        "compiled_p50_s": statistics.median(times["compiled"]),
        "compiled_p95_s": float(np.percentile(times["compiled"], 95)),
        "compiled_over_python_p50": statistics.median(times["compiled"])
        / statistics.median(times["python"]),
        "archive_max_abs_spectrum_delta_dex": float(np.max(np.abs(
            np.asarray(outputs["compiled"][0]["log10OmegaGW"], dtype=float)
            - old_spectrum
        ))),
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
            "This is one fixed cosmology, DN_eff and 800-node grid; it does not certify the complete parameter domain.",
            "It changes only the scalar H2/sigma evaluator; the ODE solver, equations, tolerances and frequency nodes are unchanged.",
            "This run validates the integrated candidate only on the recorded fixed background; multi-background and cold-start checks remain required.",
        ],
    }
    R._H2_and_sigma = _PYTHON_EVAL
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
