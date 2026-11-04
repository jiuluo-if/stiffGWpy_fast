"""Compare phase transfer with Cartesian Reference on sample2's 800-bin root grid.

中文：在 sample2 的高分辨率 Reference 根背景上逐频率比较相位传播；只有完整频谱通过时才计算候选积分映射。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import benchmark_nonosc_phase_res4_fullgrid as phase_grid  # noqa: E402
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
DN_ROOT = 2.291411867022825
FREQ_RES = 8.0
Z_TAIL = 5.0
RTOL = 1e-9
Z_START = -12.0
REF_JSON = ROOT / "docs" / "sample2_reference_root_res8_tail5_20261104.json"
REF_NPZ = ROOT / "docs" / "sample2_reference_root_res8_tail5_20261104.npz"
OUTPUT = ROOT / "docs" / "sample2_phase_reference_root_res8_tail5_20261104.json"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-gamma", type=float, default=10.0)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    if not math.isfinite(args.minimum_gamma) or args.minimum_gamma < 0.0:
        raise ValueError("minimum-gamma must be finite and nonnegative")

    phase = phase_grid.phase
    phase.PARAMS.update(PARAMS)
    phase.ZTAIL = Z_TAIL
    phase.PHASE_GAMMA_MIN = args.minimum_gamma
    model = LCDM_SG(**PARAMS, DN_eff=DN_ROOT)
    reference_source = ROOT / "stiffgwpy_fast" / "reference.py"
    reference_module_sha = hashlib.sha256(reference_source.read_bytes()).hexdigest()

    spline_warmup_start = time.perf_counter()
    R._warmup_numba_spline()
    spline_warmup_s = time.perf_counter() - spline_warmup_start
    reference_started = time.perf_counter()
    reference = R.run_reference(
        model,
        dn_eff=DN_ROOT,
        freq_res=FREQ_RES,
        z_tail=Z_TAIL,
        rtol=RTOL,
        z_start=Z_START,
        self_consistent=False,
        workers=1,
        physical_guard=False,
        backend="cartesian",
    )
    reference_seconds = time.perf_counter() - reference_started
    freqs = np.asarray(reference["freqs"], dtype=float)
    reference_logs = np.asarray(reference["log10OmegaGW"], dtype=float)
    reference_used_tail = np.asarray(reference["used_tail"], dtype=bool)
    if not (freqs.ndim == reference_logs.ndim == reference_used_tail.ndim == 1):
        raise RuntimeError("Reference returned arrays with unexpected dimensions")
    if not (freqs.size == reference_logs.size == reference_used_tail.size == 800):
        raise RuntimeError("Reference did not return the archived 800-bin res8 grid")
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(reference_logs)):
        raise FloatingPointError("Reference grid or log spectrum contains non-finite values")

    np.savez_compressed(
        REF_NPZ,
        freqs=freqs,
        reference_log10OmegaGW=reference_logs,
        reference_used_tail=reference_used_tail,
        reference_Ogw=np.asarray(reference["Ogw"], dtype=float),
        reference_Oj=np.asarray(reference["Oj"], dtype=float),
    )
    reference_payload = {
        "schema_version": 1,
        "experiment": "sample2_fixed_background_cartesian_reference_res8_tail5",
        "source_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "reference_module_sha256": reference_module_sha,
        "params": {**PARAMS, "DN_eff": DN_ROOT},
        "frequency_res": FREQ_RES,
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.tobytes()).hexdigest(),
        "z_tail": Z_TAIL,
        "z_start": Z_START,
        "rtol": RTOL,
        "backend": "cartesian",
        "spline_jit_warmup_s": spline_warmup_s,
        "runtime_s": reference_seconds,
        "reference_DN_gw": float(reference["DN_gw"]),
        "fixed_point_relative_residual": float((DN_ROOT - reference["DN_gw"]) / DN_ROOT),
        "quadrature_error": reference["quadrature_error"],
        "interpolation_error": reference["interpolation_error"],
        "tail_count": int(reference_used_tail.sum()),
        "log10OmegaGW_min": float(reference_logs.min()),
        "log10OmegaGW_max": float(reference_logs.max()),
        "npz": REF_NPZ.name,
        "npz_sha256": hashlib.sha256(REF_NPZ.read_bytes()).hexdigest(),
        "resources": telemetry(workers=1, threads=1),
        "limitations": [
            "DN_eff is the previously converged res8 tail5 Reference root anchor; this is a fixed-background map, not a fresh outer root solve.",
            "The 800-bin sampled spectrum does not certify interpolation between native nodes.",
        ],
    }
    REF_JSON.write_text(json.dumps(reference_payload, indent=2) + "\n", encoding="utf-8")

    rows = []
    phase_started = time.perf_counter()
    for index, (freq, log_omega, used_tail) in enumerate(zip(freqs, reference_logs, reference_used_tail)):
        reference_power = float(reference["Ogw"][index] - reference["Oj"][index])
        ref_row = {
            "index": index,
            "freq_log10_hz": float(freq),
            "cartesian_log10_omega": float(log_omega),
            "used_tail": bool(used_tail),
        }
        if not math.isfinite(reference_power) or reference_power <= 0.0 or log_omega <= -40.0:
            result = {
                **ref_row,
                "status": "reference_nonpositive_or_floored_spectrum",
                "reference_OmegaGW": reference_power,
            }
        else:
            result = phase_grid.run_mode(ref_row, {**PARAMS, "DN_eff": DN_ROOT})
        rows.append(result)
        if (index + 1) % 25 == 0:
            checkpoint = ROOT / ".planning" / "2026-10-30-expand-fast-domain" / "sample2_root_phase_checkpoint.jsonl"
            with checkpoint.open("w" if index == 24 else "a", encoding="utf-8") as stream:
                stream.write(json.dumps(result, allow_nan=False) + "\n")
                stream.flush()
        if (index + 1) % 50 == 0:
            print(json.dumps({"completed_modes": index + 1, "status": result["status"]}), flush=True)

    completed = [row for row in rows if row["status"] == "completed"]
    eligible = [row for row in rows if row["status"] != "reference_nonpositive_or_floored_spectrum"]
    order32_errors = [row["orders"][2]["abs_dex_error_vs_cartesian"] for row in completed]
    phase_map = None
    phase_map_diagnostics = None
    if len(completed) == len(rows):
        logs32 = np.array([row["orders"][2]["log10_omega"] for row in rows], dtype=float)
        omega32 = np.power(10.0, logs32)
        g2, quadrature_error, interpolation_error = R.integrate_spectrum(
            freqs, omega32, np.zeros_like(omega32)
        )
        omega_nu = R.gp.Omega_nh2 / model.derived_param["h"] ** 2
        phase_map = float(R.gp.Neff0 * g2 / omega_nu)
        phase_map_diagnostics = {
            "quadrature_error": quadrature_error,
            "interpolation_error": interpolation_error,
            "fixed_point_relative_residual": float((DN_ROOT - phase_map) / DN_ROOT),
        }

    payload = {
        "schema_version": 1,
        "experiment": "sample2_phase_transfer_vs_cartesian_root_res8_tail5",
        "source_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "phase_script_sha256": hashlib.sha256(Path(phase.__file__).read_bytes()).hexdigest(),
        "benchmark_driver_sha256": hashlib.sha256(Path(phase_grid.__file__).read_bytes()).hexdigest(),
        "reference_module_sha256": reference_module_sha,
        "reference_json_sha256": hashlib.sha256(REF_JSON.read_bytes()).hexdigest(),
        "reference_npz_sha256": hashlib.sha256(REF_NPZ.read_bytes()).hexdigest(),
        "params": {**PARAMS, "DN_eff": DN_ROOT},
        "frequency_res": FREQ_RES,
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.tobytes()).hexdigest(),
        "z_tail": Z_TAIL,
        "z_start": Z_START,
        "reference_rtol": RTOL,
        "minimum_gamma": args.minimum_gamma,
        "reference_runtime_s": reference_seconds,
        "phase_runtime_s": time.perf_counter() - phase_started,
        "mode_count_completed": len(completed),
        "mode_count_unsupported_or_failed": len(rows) - len(completed),
        "mode_count_reference_nonpositive_or_floored": len(rows) - len(eligible),
        "max_order32_abs_dex_error": max(order32_errors, default=None),
        "p95_order32_abs_dex_error": (
            float(np.percentile(order32_errors, 95)) if order32_errors else None
        ),
        "order32_candidate_seconds_p50": (
            float(np.percentile([row["orders"][2]["total_seconds"] for row in completed], 50))
            if completed
            else None
        ),
        "order32_candidate_seconds_p95": (
            float(np.percentile([row["orders"][2]["total_seconds"] for row in completed], 95))
            if completed
            else None
        ),
        "reference_DN_gw": float(reference["DN_gw"]),
        "phase_DN_gw": phase_map,
        "phase_map_diagnostics": phase_map_diagnostics,
        "reference_root_anchor_residual": reference_payload["fixed_point_relative_residual"],
        "resources": telemetry(workers=1, threads=1),
        "limitations": [
            "This is a fixed-background full native spectrum, not a fresh self-consistency root solve.",
            "Accuracy statistics include only modes with all three phase orders completed.",
            "Unsupported phase modes remain explicit and are never replaced by Cartesian values in candidate map statistics.",
        ],
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items() if key != "rows"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
