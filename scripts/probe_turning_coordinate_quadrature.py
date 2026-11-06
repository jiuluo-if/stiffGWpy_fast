"""Screen a horizon-turning coordinate for the low-frequency integral panel."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "lowT": ROOT / "docs" / "cr1_lowT_reference_frequency_panels_20261106.json",
    "holdout01": ROOT / "docs" / "cr1_holdout01_reference_frequency_panels_20261106.json",
}
OUTPUT = ROOT / "docs" / "turning_coordinate_quadrature_20261106.json"
PANEL_COUNTS = (1, 2, 4, 8, 16, 32, 75)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _turning_coordinate(log10_f: np.ndarray, low_edge: float) -> np.ndarray:
    """Use s=log(1+sqrt((f/f_min)^2-1)) at the present-horizon turn."""
    u = np.log(10.0) * (np.asarray(log10_f, dtype=np.float64) - float(low_edge))
    if np.any(u < -1e-14):
        raise ValueError("frequency lies below the present-horizon support edge")
    u = np.maximum(u, 0.0)
    omega = np.sqrt(np.expm1(2.0 * u))
    return np.log1p(omega)


def _jac_du_ds(s: float) -> float:
    omega = np.expm1(float(s))
    return omega * (1.0 + omega) / (1.0 + omega * omega)


def _load_case(path: Path) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    record = json.loads(path.read_text(encoding="utf-8"))
    raw_path = ROOT / record["raw_arrays"]["path"]
    if _sha(raw_path) != record["raw_arrays"]["sha256"]:
        raise RuntimeError(f"raw Reference checksum mismatch: {raw_path}")
    raw = np.load(raw_path)
    freqs = np.asarray(raw["freqs"], dtype=np.float64)
    integrand = np.asarray(raw["Ogw"] - raw["Oj"], dtype=np.float64)
    if hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest() != record["grids"]["nested_sha256"]:
        raise RuntimeError("raw nested frequency checksum differs from the reference panel record")
    if hashlib.sha256(freqs[::2].astype("<f8").tobytes()).hexdigest() != record["grids"]["native_sha256"]:
        raise RuntimeError("raw coarse frequency checksum differs from the reference panel record")
    if freqs[0] > freqs[-1]:
        freqs = freqs[::-1]
        integrand = integrand[::-1]
    if freqs.size != 151 or not np.all(np.diff(freqs) > 0.0):
        raise RuntimeError("expected the exact ascending nested 151-node support")
    return record, freqs[::2], integrand[::2], freqs[1::2], integrand[1::2]


def _screen_case(path: Path) -> dict:
    record, base_f, base_g, mid_f, mid_g = _load_case(path)
    s = _turning_coordinate(base_f, base_f[0])
    smid = _turning_coordinate(mid_f, base_f[0])
    pchip_base = PchipInterpolator(base_f, base_g)
    nested_f = np.sort(np.concatenate((base_f, mid_f)))
    nested_g = np.empty(nested_f.size, dtype=np.float64)
    nested_g[::2] = base_g
    nested_g[1::2] = mid_g
    pchip_nested = PchipInterpolator(nested_f, nested_g)
    actual_integral = float(pchip_nested.integrate(float(nested_f[0]), float(nested_f[-1])) * np.log(10.0))
    base_integral = float(pchip_base.integrate(float(base_f[0]), float(base_f[-1])) * np.log(10.0))
    nested_map = float(record["reference_151_fixed_map"]["nested_dn_gw"])
    reference_root = float(record["reference_301_root"])
    map_scale = nested_map / actual_integral

    rows = []
    for count in PANEL_COUNTS:
        count = min(count, base_f.size - 1)
        pchip_s = PchipInterpolator(s[:count + 1], base_g[:count + 1])
        transformed_integral = 0.0
        for panel in range(count):
            transformed_integral += quad(
                lambda value: float(pchip_s(value)) * _jac_du_ds(value),
                float(s[panel]), float(s[panel + 1]), epsabs=1e-30,
                epsrel=2e-12,
            )[0]
        rest_integral = 0.0
        if count < base_f.size - 1:
            rest_integral = float(
                pchip_base.integrate(float(base_f[count]), float(base_f[-1])) * np.log(10.0)
            )
        predicted = np.empty(mid_f.size, dtype=np.float64)
        for panel in range(mid_f.size):
            if panel < count:
                predicted[panel] = float(pchip_s(smid[panel]))
            else:
                predicted[panel] = float(pchip_base(mid_f[panel]))
        positive = (predicted > 0.0) & (mid_g > 0.0)
        midpoint_error = np.full(mid_f.size, np.nan)
        midpoint_error[positive] = np.abs(
            np.log10(predicted[positive]) - np.log10(mid_g[positive])
        )
        integral = transformed_integral + rest_integral
        mapped = map_scale * integral
        rows.append({
            "turning_chart_panel_count": int(count),
            "positive_midpoint_predictions": int(np.count_nonzero(positive)),
            "midpoint_count": int(mid_f.size),
            "max_midpoint_error_dex_positive_only": float(np.nanmax(midpoint_error)),
            "p95_midpoint_error_dex_positive_only": float(np.nanpercentile(midpoint_error, 95)),
            "integral_relative_difference_vs_nested_151_reference": float(integral / actual_integral - 1.0),
            "dn_gw_map_relative_difference_vs_nested_151_reference": float(mapped / nested_map - 1.0),
            "dn_gw_map_relative_difference_vs_reference_301_root": float(abs(mapped - reference_root) / abs(reference_root)),
            "transformed_boundary_integral": float(transformed_integral),
            "ordinary_pchip_remainder_integral": float(rest_integral),
        })

    return {
        "case": record["case"],
        "parameters": record["parameters"],
        "reference_tail_convention": record["reference_151_fixed_map"]["z_tail"],
        "base_grid_count": int(base_f.size),
        "nested_grid_count": int(nested_f.size),
        "base_grid_sha256": record["grids"]["native_sha256"],
        "nested_grid_sha256": record["grids"]["nested_sha256"],
        "base_pchip_map_relative_difference_vs_nested_reference": float(base_integral / actual_integral - 1.0),
        "base_pchip_map_relative_difference_vs_reference_301_root": float(abs(map_scale * base_integral - reference_root) / abs(reference_root)),
        "sweep": rows,
        "limitations": [
            "This is a quadrature-only screen on stored z_tail=10 Reference samples, not a Fast production solve.",
            "The chart changes only the endpoint interpolation variable; it does not add mode evaluations.",
            "Panel counts are diagnostic sensitivity values, not a proposed input-dependent selector.",
        ],
        "source_hashes": {
            "panel_artifact": _sha(path),
            "raw_reference_array": _sha(ROOT / record["raw_arrays"]["path"]),
        },
    }


def main() -> int:
    cases = {name: _screen_case(path) for name, path in CASES.items()}
    panel_accuracy_passes = {
        name: any(
            row["dn_gw_map_relative_difference_vs_nested_151_reference"] <= 2e-4
            and row["max_midpoint_error_dex_positive_only"] <= 1e-3
            for row in case["sweep"]
        )
        for name, case in cases.items()
    }
    output = {
        "schema_version": 1,
        "status": "complete",
        "experiment": "present-horizon turning coordinate for low-frequency spectrum integration",
        "coordinate": "u=ln(f/f_min), omega=sqrt(exp(2u)-1), s=ln(1+omega), du/ds=omega*(1+omega)/(1+omega^2)",
        "script_sha256": _sha(Path(__file__)),
        "cases": cases,
        "panel_accuracy_passes_by_case": panel_accuracy_passes,
        "decision": "PASS_SCREEN" if all(panel_accuracy_passes.values()) else "REJECTED_AT_QUADRATURE_FEASIBILITY_SCREEN",
    }
    OUTPUT.write_text(json.dumps(output, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        name: {
            "base_map_rel": case["base_pchip_map_relative_difference_vs_nested_reference"],
            "sweep": [
                {
                    "panels": row["turning_chart_panel_count"],
                    "map_rel": row["dn_gw_map_relative_difference_vs_nested_151_reference"],
                    "max_dex": row["max_midpoint_error_dex_positive_only"],
                }
                for row in case["sweep"]
            ],
        }
        for name, case in output["cases"].items()
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
