"""Estimate and check the Reference fixed point around Fast's N_inf>120 result."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import subprocess
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numba
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARENT = ROOT / "docs" / "fast_ninf122_tail8_target1e5_validation_20261105.json"
OUTPUT = ROOT / "docs" / "fast_ninf122_reference_root_secant_20261105.json"
RTOL = 1e-9
Z_TAIL = 8.0
FREQ_RES = 2.0
X1_RELATIVE_OFFSET = 1e-3


def _git(*args: str) -> str | None:
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(record: dict) -> None:
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": record["status"], "maps": list(record["maps"]), "output": str(OUTPUT)}), flush=True)


def _run(params: dict, x: float, *, freqs=None, freq_res=1.0):
    model = LCDM_SG(**params)
    start = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_reference(
            model,
            dn_eff=x,
            freq_res=freq_res,
            z_tail=Z_TAIL,
            rtol=RTOL,
            freq_subset=freqs,
            self_consistent=False,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
    if caught:
        raise RuntimeError(f"Reference emitted warnings at DN_eff={x}: {[str(w.message) for w in caught]}")
    values = np.asarray(result["freqs"], dtype=float)
    raw = np.asarray(result["Ogw"], dtype=float) - np.asarray(result["Oj"], dtype=float)
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(raw)):
        raise FloatingPointError(f"non-finite Reference map at DN_eff={x}")
    return result, values, raw, time.perf_counter() - start


def _summary(result: dict, freqs: np.ndarray, raw: np.ndarray, elapsed: float) -> dict:
    return {
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "frequencies_log10_hz": freqs.tolist(),
        "omega_gw": raw.tolist(),
        "DN_gw": float(result["DN_gw"]),
        "quadrature_error": float(result["quadrature_error"]),
        "interpolation_error": float(result["interpolation_error"]),
        "runtime_s": elapsed,
        "used_tail_count": int(np.count_nonzero(result["used_tail"])),
    }


def main() -> None:
    parent_sha = _sha(PARENT)
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    params = parent["parameters"]
    x0 = float(parent["fast"]["DN_eff_final"])
    x1 = x0 * (1.0 + X1_RELATIVE_OFFSET)
    base_map = parent["reference_maps"]["native_res2"]
    base_adaptive_map = parent["reference_maps"]["adaptive_fast_grid"]
    base_freqs = np.asarray(base_map["frequencies_log10_hz"], dtype=float)
    base_hash = str(base_map["frequency_sha256"])
    adaptive_freqs = np.asarray(parent["fast"]["frequencies_log10_hz"], dtype=float)
    record = {
        "schema_version": 1,
        "experiment": "fast_ninf122_reference_fixed_point_secant_and_check",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": params,
        "equation": "x = DN_eff_input + DN_gw(x)",
        "DN_eff_input": float(params["DN_eff"]),
        "parent_artifact": PARENT.name,
        "parent_artifact_sha256": parent_sha,
        "fast_x0": x0,
        "x1_relative_offset": X1_RELATIVE_OFFSET,
        "x1": x1,
        "reference": {"backend": "prufer", "rtol": RTOL, "z_tail": Z_TAIL, "workers": 1, "native_freq_res": FREQ_RES},
        "x0_native_res2": {
            "frequency_count": int(base_freqs.size),
            "frequency_sha256": base_hash,
            "DN_gw": float(base_map["DN_gw"]),
            "runtime_s": float(base_map["runtime_s"]),
        },
        "g_x0_adaptive_grid": x0 - float(base_adaptive_map["DN_gw"]),
        "maps": {},
        "status": "x1_native_res2_running",
    }
    _save(record)

    print(f"Native res2 Reference map at x1={x1:.16g}", flush=True)
    r1, f1, s1, t1 = _run(params, x1, freq_res=FREQ_RES)
    hash1 = hashlib.sha256(f1.astype("<f8").tobytes()).hexdigest()
    if hash1 != base_hash:
        raise RuntimeError(f"native grid changed with DN_eff: {hash1} != {base_hash}")
    record["maps"]["x1_native_res2"] = _summary(r1, f1, s1, t1)
    f0 = float(base_map["DN_gw"])
    f_x1 = float(r1["DN_gw"])
    slope = (f_x1 - f0) / (x1 - x0)
    g0 = x0 - f0
    g1 = x1 - f_x1
    denominator = 1.0 - slope
    if not math.isfinite(denominator) or denominator == 0.0:
        raise FloatingPointError(f"degenerate fixed-point secant slope {slope}")
    x_sec = x0 - g0 / denominator
    if not math.isfinite(x_sec) or x_sec <= 0.0:
        raise FloatingPointError(f"invalid secant root estimate {x_sec}")
    record.update({
        "secant_slope": slope,
        "g_x0": g0,
        "g_x1": g1,
        "x_secant": x_sec,
        "secant_relative_shift_from_fast_x0": abs(x_sec - x0) / x0,
        "x0_x1_sign_change": bool(g0 * g1 <= 0.0),
        "status": "secant_estimate_ready",
    })
    _save(record)

    print(f"Native res2 Reference map at secant candidate x={x_sec:.16g}", flush=True)
    r_sec, f_sec, s_sec, t_sec = _run(params, x_sec, freq_res=FREQ_RES)
    hash_sec = hashlib.sha256(f_sec.astype("<f8").tobytes()).hexdigest()
    if hash_sec != base_hash:
        raise RuntimeError(f"native grid changed at secant candidate: {hash_sec} != {base_hash}")
    record["maps"]["x_secant_native_res2"] = _summary(r_sec, f_sec, s_sec, t_sec)
    g_sec = x_sec - float(r_sec["DN_gw"])
    record["g_x_secant_native_res2"] = g_sec
    record["secant_corrected_relative_root_residual"] = abs(g_sec / denominator) / x_sec
    record["status"] = "secant_native_res2_complete"
    _save(record)

    print(f"Adaptive-grid Reference map at secant candidate x={x_sec:.16g}", flush=True)
    r_adapt, f_adapt, s_adapt, t_adapt = _run(params, x_sec, freqs=adaptive_freqs)
    record["maps"]["x_secant_adaptive_grid"] = _summary(r_adapt, f_adapt, s_adapt, t_adapt)
    record["adaptive_vs_native_res2_relative_DN_gw_at_secant"] = abs(
        float(r_adapt["DN_gw"]) - float(r_sec["DN_gw"])
    ) / max(abs(float(r_adapt["DN_gw"])), abs(float(r_sec["DN_gw"])), np.finfo(float).tiny)
    g_adapt = x_sec - float(r_adapt["DN_gw"])
    record["g_x_secant_adaptive_grid"] = g_adapt
    g0_adapt = float(record["g_x0_adaptive_grid"])
    left_signs_agree = g0 * g0_adapt > 0.0
    right_signs_agree = g_sec * g_adapt > 0.0
    f0_adapt = float(base_adaptive_map["DN_gw"])
    f_sec_adapt = float(r_adapt["DN_gw"])
    adaptive_slope = (f_sec_adapt - f0_adapt) / (x_sec - x0)
    adaptive_denominator = 1.0 - adaptive_slope
    if not math.isfinite(adaptive_denominator) or adaptive_denominator == 0.0:
        raise FloatingPointError(f"degenerate adaptive secant slope {adaptive_slope}")
    x_adaptive_secant = x0 - g0_adapt / adaptive_denominator
    two_grid_root_delta = abs(x_sec - x_adaptive_secant) / max(
        abs(x_sec), abs(x_adaptive_secant), np.finfo(float).tiny
    )
    record["adaptive_secant_slope"] = adaptive_slope
    record["x_adaptive_secant"] = x_adaptive_secant
    record["adaptive_secant_corrected_relative_root_residual"] = abs(
        g_adapt / adaptive_denominator
    ) / x_adaptive_secant
    record["two_grid_secant_root_relative_delta"] = two_grid_root_delta
    record["secant_candidate_bracket"] = {
        "left": min(x0, x_sec),
        "right": max(x0, x_sec),
        "g_left_native_res2": g0 if x0 <= x_sec else g_sec,
        "g_right_native_res2": g_sec if x_sec >= x0 else g0,
        "relative_width": abs(x_sec - x0) / x0,
        "sign_change": bool(g0 * g_sec <= 0.0),
        "left_sign_agrees_across_grids": left_signs_agree,
        "right_sign_agrees_across_grids": right_signs_agree,
        "two_grid_sign_change": bool(left_signs_agree and right_signs_agree and g0 * g_sec <= 0.0),
        "two_grid_width_le_2e-4": abs(x_sec - x0) / x0 <= 2e-4,
        "adaptive_grid_sign_change_x0_to_x_secant": bool(g0_adapt * g_adapt <= 0.0),
        "adaptive_grid_width_le_2e-4": abs(x_sec - x0) / x0 <= 2e-4,
        "two_grid_secant_roots_agree_le_2e-4": two_grid_root_delta <= 2e-4,
    }
    record["hashes"] = {
        "fast_sgwb_py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "reference_py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "script_py": _sha(Path(__file__)),
    }
    record["resources"] = {
        **telemetry(threads=1),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "numba": numba.__version__,
    }
    record["limitations"] = [
        "The derivative is an interval secant, not a local derivative or rigorous root bound.",
        "Two-grid convergence is checked at the secant candidate, but Reference frequency error is not a formal interval enclosure.",
        "No production Fast code or Reference performance path is modified by this diagnostic.",
    ]
    record["status"] = "complete"
    _save(record)
    print(json.dumps({
        "output": str(OUTPUT),
        "x0": x0,
        "x1": x1,
        "F_x0": f0,
        "F_x1": f_x1,
        "secant_slope": slope,
        "x_secant": x_sec,
        "bracket": record["secant_candidate_bracket"],
        "g_x_secant_native_res2": g_sec,
        "g_x_secant_adaptive": g_adapt,
        "relative_grid_delta_at_secant": record["adaptive_vs_native_res2_relative_DN_gw_at_secant"],
        "x_adaptive_secant": x_adaptive_secant,
        "two_grid_secant_root_relative_delta": two_grid_root_delta,
        "adaptive_grid_x0_xsecant_sign_change": record["secant_candidate_bracket"]["adaptive_grid_sign_change_x0_to_x_secant"],
    }, indent=2))


if __name__ == "__main__":
    main()
