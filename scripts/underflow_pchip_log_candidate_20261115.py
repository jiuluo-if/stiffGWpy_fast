"""Measure log-PCHIP interpolation bias against direct 469/937 Reference nodes."""

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
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_pchip_log_candidate_20261115.json"

from stiffgwpy_fast import global_param as gp  # noqa: E402


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _new_nodes(coarse: np.ndarray, fine: np.ndarray) -> np.ndarray:
    values = set(np.asarray(coarse, dtype=np.float64).tolist())
    return np.asarray([x for x in fine if float(x) not in values], dtype=np.float64)


def _error_summary(actual_log: np.ndarray, predicted_log: np.ndarray) -> dict:
    errors = np.abs(actual_log - predicted_log) / math.log(10.0)
    return {
        "max_abs_error_dex": float(np.max(errors)),
        "p95_abs_error_dex": float(np.quantile(errors, 0.95)),
        "count": int(errors.size),
    }


def _log_pchip_integral(freqs: np.ndarray, log_integrand_ln: np.ndarray) -> dict:
    x = np.asarray(freqs, dtype=np.float64)[::-1]
    y = np.asarray(log_integrand_ln, dtype=np.float64)[::-1]
    interpolator = PchipInterpolator(x, y)
    shift = float(np.max(y)) + 1.0
    value, error = quad(
        lambda point: math.exp(float(interpolator(point)) - shift),
        float(x[0]),
        float(x[-1]),
        epsabs=0.0,
        epsrel=1e-11,
        limit=max(800, x.size + 50),
        points=x,
    )
    log_g2 = math.log(math.log(10.0)) + shift + math.log(value)
    return {
        "log_g2": log_g2,
        "quadrature_relative_error": float(error / value),
    }


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    grids = {
        int(n): np.asarray(row["log10_hz"], dtype=np.float64)
        for n, row in data["frequency_grids"].items()
    }
    reference = data["reference_fixed_background"]
    common_shift = max(
        float(np.max(reference[n]["log_integrand_ln"])) for n in ("235", "469", "937")
    )
    results = {}
    for coarse_n, fine_n in ((235, 469), (469, 937), (235, 937)):
        coarse_x = grids[coarse_n]
        fine_x = grids[fine_n]
        new_x = _new_nodes(coarse_x, fine_x)
        coarse_logs = np.asarray(reference[str(coarse_n)]["log_integrand_ln"], dtype=np.float64)
        fine_logs_by_x = dict(zip(
            np.asarray(fine_x, dtype=np.float64).tolist(),
            np.asarray(reference[str(fine_n)]["log_integrand_ln"], dtype=np.float64).tolist(),
        ))
        actual_log = np.asarray([fine_logs_by_x[float(x)] for x in new_x], dtype=np.float64)
        log_interpolator = PchipInterpolator(coarse_x[::-1], coarse_logs[::-1])
        log_prediction = np.asarray(log_interpolator(new_x), dtype=np.float64)

        coarse_normalized = np.exp(coarse_logs - common_shift)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            linear_interpolator = PchipInterpolator(coarse_x[::-1], coarse_normalized[::-1])
            linear_prediction = np.asarray(linear_interpolator(new_x), dtype=np.float64)
        finite_mask = linear_prediction > 0.0
        linear_error = (
            _error_summary(
                actual_log[finite_mask],
                np.log(linear_prediction[finite_mask]) + common_shift,
            )
            if np.any(finite_mask)
            else None
        )
        log_result = _log_pchip_integral(coarse_x, coarse_logs)
        coarse_map = float(reference[str(coarse_n)]["mapped_DN_gw"])
        fine_map = float(reference[str(fine_n)]["mapped_DN_gw"])
        h = float(data["parameters"]["H0"]) / 100.0
        omega_nu = gp.Omega_nh2 / (h * h)
        log_candidate_map = log_result["log_g2"] + math.log(gp.Neff0 / omega_nu)
        log_map_error = log_candidate_map - math.log(coarse_map)
        results[f"{coarse_n}_to_{fine_n}"] = {
            "new_node_count": int(new_x.size),
            "log_pchip_new_node_error_dex": _error_summary(actual_log, log_prediction),
            "linear_pchip_zero_prediction_count": int(np.count_nonzero(linear_prediction == 0.0)),
            "linear_pchip_positive_new_node_error_dex": linear_error,
            "linear_pchip_runtime_warnings": [str(item.message) for item in caught],
            "log_pchip_integral_relative_to_direct_coarse_map": float(math.expm1(log_map_error)),
            "direct_fine_map_relative_to_direct_coarse_map": float(fine_map / coarse_map - 1.0),
            "log_pchip_quadrature_relative_error": log_result["quadrature_relative_error"],
        }

    output = {
        "schema_version": 1,
        "status": "complete_log_pchip_bias_diagnostic_not_adopted",
        "experiment": "Log-PCHIP versus direct nested Reference nodes; linear PCHIP retained as baseline",
        "parameters": data["parameters"],
        "background": data["background"],
        "frequency_grid_hashes_f64le": {
            str(n): data["frequency_grids"][str(n)]["sha256_f64le"] for n in grids
        },
        "comparisons": results,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "source_sha256": {
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "input_json": _sha(INPUT),
            "diagnostic_script": _sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), "comparisons": results}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
