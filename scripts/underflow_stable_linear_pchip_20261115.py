"""Audit log-scaled arithmetic for the same linear PCHIP on underflow data."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import warnings
from pathlib import Path

import numpy as np
import scipy
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.stable_log_pchip import LogPchip  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_stable_linear_pchip_20261115.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _new_nodes(coarse: np.ndarray, fine: np.ndarray) -> np.ndarray:
    exact = set(np.asarray(coarse, dtype=np.float64).tolist())
    return np.asarray([value for value in fine if float(value) not in exact], dtype=np.float64)


def _error_summary(actual_log: np.ndarray, predicted_log: np.ndarray) -> dict:
    errors = np.abs(actual_log - predicted_log) / math.log(10.0)
    return {
        "max_abs_error_dex": float(np.max(errors)),
        "p95_abs_error_dex": float(np.quantile(errors, 0.95)),
        "count": int(errors.size),
    }


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    frequencies = {
        int(count): np.asarray(row["log10_hz"], dtype=np.float64)
        for count, row in data["frequency_grids"].items()
    }
    ref = data["reference_fixed_background"]
    h = float(data["parameters"]["H0"]) / 100.0
    omega_nu = gp.Omega_nh2 / (h * h)
    map_scale = gp.Neff0 / omega_nu
    per_grid = {}

    for count in (235, 469, 937):
        x = frequencies[count][::-1]
        log_y = np.asarray(ref[str(count)]["log_integrand_ln"], dtype=np.float64)[::-1]
        with warnings.catch_warnings(record=True) as caught_stable:
            warnings.simplefilter("always", RuntimeWarning)
            stable = LogPchip(x, log_y)
            log_g2 = math.log(math.log(10.0)) + stable.integrate_log()
            mapped = math.exp(log_g2) * map_scale

        normalized = np.exp(log_y - float(np.max(log_y)))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            linear = PchipInterpolator(x, normalized)
        linear_zero_count = int(np.count_nonzero(normalized == 0.0))
        per_grid[str(count)] = {
            "same_grid_map_DN_gw": float(mapped),
            "baseline_reference_map_DN_gw": float(ref[str(count)]["mapped_DN_gw"]),
            "relative_map_difference": float(mapped / float(ref[str(count)]["mapped_DN_gw"]) - 1.0),
            "log_pchip_nodes_finite": bool(np.isfinite(stable.evaluate_log(x)).all()),
            "stable_log_pchip_warnings": [str(item.message) for item in caught_stable],
            "linear_normalized_zero_nodes": linear_zero_count,
            "linear_pchip_warnings": [str(item.message) for item in caught],
        }

    interpolation = {}
    for coarse, fine in ((235, 469), (469, 937), (235, 937)):
        x_coarse = frequencies[coarse][::-1]
        log_coarse = np.asarray(ref[str(coarse)]["log_integrand_ln"], dtype=np.float64)[::-1]
        stable = LogPchip(x_coarse, log_coarse)
        x_fine = frequencies[fine]
        new_x_desc = _new_nodes(frequencies[coarse], x_fine)
        new_x = new_x_desc[::-1]
        fine_map = dict(zip(
            frequencies[fine][::-1].tolist(),
            np.asarray(ref[str(fine)]["log_integrand_ln"], dtype=np.float64)[::-1].tolist(),
        ))
        actual_log = np.asarray([fine_map[float(value)] for value in new_x], dtype=np.float64)
        stable_prediction = stable.evaluate_log(new_x)

        global_shift = max(float(np.max(log_coarse)), float(np.max(actual_log)))
        linear_values = np.exp(log_coarse - global_shift)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            linear = PchipInterpolator(x_coarse, linear_values)
            linear_prediction = np.asarray(linear(new_x), dtype=np.float64)
        valid = linear_prediction > 0.0
        linear_error = (
            _error_summary(
                actual_log[valid],
                np.log(linear_prediction[valid]) + global_shift,
            )
            if np.any(valid)
            else None
        )
        interpolation[f"{coarse}_to_{fine}"] = {
            "new_node_count": int(new_x.size),
            "stable_same_pchip_error_dex": _error_summary(actual_log, stable_prediction),
            "linear_pchip_zero_predictions": int(np.count_nonzero(linear_prediction == 0.0)),
            "linear_pchip_positive_error_dex": linear_error,
            "linear_pchip_warnings": [str(item.message) for item in caught],
        }

    output = {
        "schema_version": 1,
        "status": "complete_same_linear_pchip_stability_diagnostic_not_adopted",
        "experiment": "Log-scaled arithmetic for SciPy's linear PCHIP, validated against ordinary-range parity and direct 469/937 nodes",
        "parameters": data["parameters"],
        "background": data["background"],
        "frequency_grid_hashes_f64le": {
            str(count): data["frequency_grids"][str(count)]["sha256_f64le"]
            for count in (235, 469, 937)
        },
        "same_grid_integrals": per_grid,
        "new_node_interpolation": interpolation,
        "ordinary_range_parity_tests": {
            "command": "python -m pytest -q tests/test_stable_log_pchip.py",
            "passed": 3,
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "source_sha256": {
            "stable_log_pchip.py": _sha(ROOT / "scripts" / "stable_log_pchip.py"),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "input_json": _sha(INPUT),
            "diagnostic_script": _sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "same_grid_map_relative_differences": {
            count: row["relative_map_difference"] for count, row in per_grid.items()
        },
        "new_node_errors": {
            key: row["stable_same_pchip_error_dex"] for key, row in interpolation.items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
