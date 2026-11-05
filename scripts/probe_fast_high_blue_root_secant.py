"""Two-grid secant screen for the high-blue Fast fixed point."""

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
BASELINE = ROOT / "docs" / "fast_high_blue_target1e7_prufer_20261105.json"
OUTPUT = ROOT / "docs" / "fast_high_blue_root_secant_20261105.json"
OFFSET = 1e-3
RTOL = 1e-9
Z_TAIL = 8.0
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _git(*args: str) -> str | None:
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(record: dict) -> None:
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": record["status"], "maps": list(record["maps"]), "output": str(OUTPUT)}), flush=True)


def _map(model: LCDM_SG, dn_eff: float, *, freqs=None, freq_res=1.0) -> dict:
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_reference(
            model,
            dn_eff=dn_eff,
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
        raise RuntimeError(f"Reference emitted warnings at DN_eff={dn_eff}: {[str(item.message) for item in caught]}")
    f = np.asarray(result["freqs"], dtype=float)
    y = np.asarray(result["Ogw"], dtype=float) - np.asarray(result["Oj"], dtype=float)
    if not np.all(np.isfinite(f)) or not np.all(np.isfinite(y)) or not np.all(y > 0.0):
        raise FloatingPointError(f"invalid Reference map at DN_eff={dn_eff}")
    return {
        "frequency_count": int(f.size),
        "frequency_sha256": hashlib.sha256(f.astype("<f8").tobytes()).hexdigest(),
        "frequencies_log10_hz": f.tolist(),
        "omega_gw": y.tolist(),
        "DN_gw": float(result["DN_gw"]),
        "quadrature_error": float(result["quadrature_error"]),
        "interpolation_error": float(result["interpolation_error"]),
        "runtime_s": time.perf_counter() - started,
    }


def _secant_root(x0: float, x1: float, f0: float, f1: float) -> tuple[float, float, float, float]:
    g0, g1 = x0 - f0, x1 - f1
    slope = (f1 - f0) / (x1 - x0)
    denominator = 1.0 - slope
    if not math.isfinite(denominator) or denominator == 0.0:
        raise FloatingPointError(f"degenerate fixed-point secant denominator {denominator}")
    root = x0 - g0 / denominator
    if not math.isfinite(root) or root <= 0.0:
        raise FloatingPointError(f"invalid secant root estimate {root}")
    return root, slope, g0, g1


def main() -> None:
    baseline_sha = _sha(BASELINE)
    base = json.loads(BASELINE.read_text(encoding="utf-8"))
    x0 = float(base["fast"]["DN_eff_final"])
    x1 = x0 * (1.0 + OFFSET)
    adaptive_f0 = float(base["reference_maps"]["adaptive_exact_fast_grid"]["DN_gw"])
    res4_f0 = float(base["reference_maps"]["native_res4"]["DN_gw"])
    freqs = np.asarray(base["fast"]["frequencies_log10_hz"], dtype=float)
    model = LCDM_SG(**base["parameters"])
    record = {
        "schema_version": 1,
        "experiment": "fast_high_blue_reference_fixed_point_two_grid_secant_screen",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": base["parameters"],
        "derived": base["derived"],
        "equation": "x = DN_eff_input + DN_gw(x)",
        "DN_eff_input": float(base["parameters"]["DN_eff"]),
        "fast_x0": x0,
        "x1_relative_offset": OFFSET,
        "x1": x1,
        "parent_artifact": BASELINE.name,
        "parent_artifact_sha256": baseline_sha,
        "reference": {"backend": "prufer", "rtol": RTOL, "z_tail": Z_TAIL, "workers": 1, "native_res4": 4.0},
        "x0_maps": {
            "adaptive_fast_grid": {"DN_gw": adaptive_f0, "runtime_s": base["reference_maps"]["adaptive_exact_fast_grid"]["runtime_s"]},
            "native_res2": {"DN_gw": base["reference_maps"]["native_res2"]["DN_gw"]},
            "native_res4": {"DN_gw": res4_f0, "runtime_s": base["reference_maps"]["native_res4"]["runtime_s"]},
        },
        "maps": {},
        "status": "adaptive_x1_running",
    }
    _save(record)

    print(f"Adaptive Prüfer map at x1={x1:.16g}", flush=True)
    adaptive_x1 = _map(model, x1, freqs=freqs)
    if adaptive_x1["frequency_sha256"] != base["fast"]["frequency_sha256"]:
        raise AssertionError("adaptive frequency grid changed with DN_eff")
    record["maps"]["x1_adaptive_fast_grid"] = adaptive_x1
    record["status"] = "adaptive_x1_complete_res4_running"
    _save(record)

    print(f"Native res4 Prüfer map at x1={x1:.16g}", flush=True)
    res4_x1 = _map(model, x1, freq_res=4.0)
    record["maps"]["x1_native_res4"] = res4_x1
    adaptive_root, adaptive_slope, adaptive_g0, adaptive_g1 = _secant_root(
        x0, x1, adaptive_f0, float(adaptive_x1["DN_gw"])
    )
    res4_root, res4_slope, res4_g0, res4_g1 = _secant_root(
        x0, x1, res4_f0, float(res4_x1["DN_gw"])
    )
    root_delta = abs(adaptive_root - res4_root) / max(abs(adaptive_root), abs(res4_root), np.finfo(float).tiny)
    record["root_estimates"] = {
        "adaptive_secant_root": adaptive_root,
        "native_res4_secant_root": res4_root,
        "adaptive_secant_slope": adaptive_slope,
        "native_res4_secant_slope": res4_slope,
        "adaptive_g_x0": adaptive_g0,
        "adaptive_g_x1": adaptive_g1,
        "native_res4_g_x0": res4_g0,
        "native_res4_g_x1": res4_g1,
        "two_grid_root_relative_delta": root_delta,
        "adaptive_sign_change": adaptive_g0 * adaptive_g1 <= 0.0,
        "native_res4_sign_change": res4_g0 * res4_g1 <= 0.0,
        "native_res4_frequency_grid_hash_stable": (
            res4_x1["frequency_sha256"] == base["reference_maps"]["native_res4"]["frequency_sha256"]
        ),
        "native_res4_frequency_count_x0": int(base["reference_maps"]["native_res4"]["frequency_count"]),
        "native_res4_frequency_count_x1": int(res4_x1["frequency_count"]),
        "fast_to_res4_root_relative_shift": abs(x0 - res4_root) / x0,
    }
    record["status"] = "x1_maps_complete_candidate_residual_running"
    _save(record)

    print(f"Adaptive Prüfer map at res4 secant candidate x={res4_root:.16g}", flush=True)
    candidate = _map(model, res4_root, freqs=freqs)
    candidate_g = res4_root - float(candidate["DN_gw"])
    corrected = abs(candidate_g / (1.0 - adaptive_slope)) / res4_root
    record["maps"]["native_res4_secant_candidate_on_adaptive_grid"] = candidate
    record["root_estimates"]["adaptive_residual_at_native_res4_candidate"] = candidate_g
    record["root_estimates"]["adaptive_corrected_relative_root_residual"] = corrected
    record["acceptance"] = {
        "two_grid_secant_roots_agree_le_2e-4": root_delta <= 2e-4,
        "adaptive_candidate_corrected_residual_le_2e-4": corrected <= 2e-4,
    }
    record["status"] = "two_grid_secant_screen_complete"
    record["resources"] = {
        **telemetry(threads=1),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "numba": numba.__version__,
        "reference_workers": 1,
    }
    record["hashes"] = {"reference_py": _sha(ROOT / "stiffgwpy_fast/reference.py"), "script_py": _sha(Path(__file__))}
    record["limitations"] = [
        "This is a two-grid interval-secant screen, not a rigorous root enclosure or full iterative Reference solution.",
        "Res4 frequency convergence is checked at x0; root sensitivity compares res4 scalar maps, whose frequency nodes may shift with DN_eff.",
    ]
    _save(record)
    print(json.dumps({"status": record["status"], "root_estimates": record["root_estimates"], "acceptance": record["acceptance"], "output": str(OUTPUT)}), flush=True)


if __name__ == "__main__":
    main()
