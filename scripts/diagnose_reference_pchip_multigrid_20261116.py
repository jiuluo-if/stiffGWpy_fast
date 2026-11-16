"""Compare stable linear-PCHIP arithmetic on direct Reference grids."""

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

from scripts.stable_log_pchip import LogPchip, _signed_log_sum  # noqa: E402
from stiffgwpy_fast.reference import integrate_spectrum  # noqa: E402

RESULTS = ROOT / "docs" / "research" / "results"
GRID_235_469 = RESULTS / "underflow_min_subnormal_reference_grid_20261115.json"
REFERENCE_937 = RESULTS / "underflow_min_subnormal_reference_z10_20261115.json"
LOW_T_GRID = RESULTS / "adaptive_handoff_background_holdout_20261115.json"
LOW_T_REFERENCE = RESULTS / "adaptive_handoff_background_z10_20261115.json"
OUTPUT = RESULTS / "reference_pchip_multigrid_diagnostic_20261116.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def diagnose_case(
    label: str,
    x_raw: np.ndarray,
    log10_y_raw: np.ndarray,
    *,
    parameters: dict[str, object],
    frozen_dn: float,
    z_tail: float,
    rtol: float,
    source_sha256: dict[str, str],
) -> dict[str, object]:
    x_raw = np.asarray(x_raw, dtype=np.float64)
    log10_y_raw = np.asarray(log10_y_raw, dtype=np.float64)
    if x_raw.ndim != 1 or x_raw.size < 4 or x_raw.shape != log10_y_raw.shape:
        raise ValueError(f"{label}: frequency and spectrum arrays are incompatible")
    if not np.isfinite(x_raw).all() or not np.isfinite(log10_y_raw).all():
        raise ValueError(f"{label}: source nodes contain non-finite values")
    order = np.argsort(x_raw)
    x = x_raw[order]
    log10_y = log10_y_raw[order]
    if np.any(np.diff(x) <= 0.0):
        raise ValueError(f"{label}: frequency coordinates are not unique")
    log_y = log10_y * math.log(10.0)
    shift = float(np.max(log_y))
    with np.errstate(under="ignore"):
        y = np.exp(log_y - shift)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        g2, qerr, ierr = integrate_spectrum(x, y, np.zeros_like(y))
    original_log_integral = shift + math.log(g2)

    stable = LogPchip(x, log_y)
    stable_log_integral = math.log(math.log(10.0)) + stable.integrate_log()
    panel_logs = []
    for i, width in enumerate(np.diff(x)):
        log_width = math.log(float(width))
        terms = [
            (1, log_width - math.log(2.0) + float(log_y[i])),
            (1, log_width - math.log(2.0) + float(log_y[i + 1])),
        ]
        if stable.derivative_sign[i] != 0:
            terms.append((int(stable.derivative_sign[i]), 2.0 * log_width - math.log(12.0) + float(stable.derivative_log[i])))
        if stable.derivative_sign[i + 1] != 0:
            terms.append((-int(stable.derivative_sign[i + 1]), 2.0 * log_width - math.log(12.0) + float(stable.derivative_log[i + 1])))
        sign, log_area = _signed_log_sum(terms)
        if sign <= 0 or not math.isfinite(log_area):
            raise FloatingPointError(f"{label}: stable PCHIP panel {i} has non-positive area")
        panel_logs.append(log_area)
    panel_logs = np.asarray(panel_logs, dtype=np.float64)
    zero_adjacent = (y[:-1] == 0.0) | (y[1:] == 0.0)
    zero_panel_log_fraction = float(
        np.logaddexp.reduce(panel_logs[zero_adjacent]) - np.logaddexp.reduce(panel_logs)
    ) if zero_adjacent.any() else None

    midpoint = (x[:-1] + x[1:]) / 2.0
    with warnings.catch_warnings(record=True) as midpoint_warnings:
        warnings.simplefilter("always", RuntimeWarning)
        linear_midpoint = PchipInterpolator(x, y)(midpoint)
    stable_log_midpoint = stable.evaluate_log(midpoint)
    representable = linear_midpoint > 0.0
    error_dex = np.full(midpoint.shape, np.nan, dtype=np.float64)
    error_dex[representable] = (
        np.log10(linear_midpoint[representable])
        - (stable_log_midpoint[representable] - shift) / math.log(10.0)
    )
    worst = int(np.nanargmax(np.abs(error_dex))) if representable.any() else None
    return {
        "label": label,
        "frequency_count": int(x.size),
        "frequency_sha256_f64le": digest(x),
        "parameters": parameters,
        "frozen_DN_eff": frozen_dn,
        "z_tail": z_tail,
        "rtol": rtol,
        "log10_signal_range_dex": float(np.max(log10_y) - np.min(log10_y)),
        "normalized_zero_nodes": int(np.count_nonzero(y == 0.0)),
        "normalized_positive_subnormal_nodes": int(np.count_nonzero((y > 0.0) & (y < np.finfo(np.float64).tiny))),
        "original_linear_pchip": {
            "warnings": [str(item.message) for item in caught],
            "log_integral": original_log_integral,
            "quadrature_error_relative": float(qerr / g2),
            "half_density_interpolation_error_relative": float(ierr / g2),
        },
        "stable_log_arithmetic_same_linear_pchip": {
            "warnings": [],
            "log_integral": stable_log_integral,
            "relative_integral_difference": math.expm1(stable_log_integral - original_log_integral),
            "zero_adjacent_panel_count": int(np.count_nonzero(zero_adjacent)),
            "log10_fraction_total_integral_in_zero_adjacent_panels": (
                None if zero_panel_log_fraction is None else zero_panel_log_fraction / math.log(10.0)
            ),
        },
        "midpoint_interpolation": {
            "count": int(midpoint.size),
            "linear_zero_count": int(np.count_nonzero(~representable)),
            "warnings": [str(item.message) for item in midpoint_warnings],
            "finite_linear_vs_stable_max_abs_dex": float(np.nanmax(np.abs(error_dex))),
            "finite_linear_vs_stable_p95_abs_dex": float(np.nanpercentile(np.abs(error_dex[representable]), 95.0)),
            "worst_finite": None if worst is None else {
                "log10_frequency": float(midpoint[worst]),
                "difference_dex": float(error_dex[worst]),
                "left_node_log10_signal": float(log10_y[worst]),
                "right_node_log10_signal": float(log10_y[worst + 1]),
            },
        },
        "input_sha256": source_sha256,
    }


def main() -> int:
    ref_sha = sha(ROOT / "stiffgwpy_fast" / "reference.py")
    cases: list[dict[str, object]] = []

    grid_data = json.loads(GRID_235_469.read_text(encoding="utf-8"))
    if grid_data["source_sha256"]["reference.py"] != ref_sha:
        raise RuntimeError("469-node Reference source SHA differs from current Reference")
    case_469 = next(case for case in grid_data["cases"] if case["frequency_count"] == 469)
    cases.append(diagnose_case(
        "cr0_minimum_subnormal_469_direct_z8",
        np.asarray(case_469["frequency_log10_hz"], dtype=np.float64),
        np.asarray(case_469["reference"]["log10_OmegaGW"], dtype=np.float64),
        parameters=grid_data["parameters"],
        frozen_dn=float(case_469["frozen_DN_eff"]),
        z_tail=float(case_469["reference"]["z_tail"]),
        rtol=float(case_469["reference"]["rtol"]),
        source_sha256={"grid_and_reference_input": sha(GRID_235_469)},
    ))

    ref_937 = json.loads(REFERENCE_937.read_text(encoding="utf-8"))
    if ref_937["source_sha256"]["reference.py"] != ref_sha:
        raise RuntimeError("937-node Reference source SHA differs from current Reference")
    case_937 = next(case for case in grid_data["cases"] if case["frequency_count"] == 937)
    grid_937 = np.asarray(case_937["frequency_log10_hz"], dtype=np.float64)
    if digest(grid_937) != grid_data["frequency_grids"]["937"]["sha256_f64le"]:
        raise RuntimeError("canonical 937-node grid data fails the recorded grid hash")
    cases.append(diagnose_case(
        "cr0_minimum_subnormal_937_direct_z8",
        grid_937,
        np.asarray(case_937["reference"]["log10_OmegaGW"], dtype=np.float64),
        parameters=grid_data["parameters"],
        frozen_dn=float(case_937["frozen_DN_eff"]),
        z_tail=float(case_937["reference"]["z_tail"]),
        rtol=float(case_937["reference"]["rtol"]),
        source_sha256={"grid_and_reference_input": sha(GRID_235_469)},
    ))
    if digest(grid_937) != ref_937["frequency_sha256_f64le"]:
        raise RuntimeError("937-node Reference spectrum grid does not match its direct grid")
    cases.append(diagnose_case(
        "cr0_minimum_subnormal_937_direct_z10",
        grid_937,
        np.asarray(ref_937["reference"]["log10_OmegaGW"], dtype=np.float64),
        parameters=ref_937["parameters"],
        frozen_dn=float(ref_937["frozen_DN_eff"]),
        z_tail=float(ref_937["reference"]["z_tail"]),
        rtol=float(ref_937["reference"]["rtol"]),
        source_sha256={"grid_input": sha(GRID_235_469), "reference_input": sha(REFERENCE_937)},
    ))

    low_t_grid = json.loads(LOW_T_GRID.read_text(encoding="utf-8"))
    low_t_ref = json.loads(LOW_T_REFERENCE.read_text(encoding="utf-8"))
    if low_t_ref["source_sha256"]["reference.py"] != ref_sha:
        raise RuntimeError("low-temperature Reference source SHA differs from current Reference")
    freq_221 = np.asarray(low_t_grid["frequency_log10_hz"], dtype=np.float64)
    if digest(freq_221) != low_t_ref["frequency_sha256_f64le"]:
        raise RuntimeError("low-temperature direct Reference grid hash mismatch")
    cases.append(diagnose_case(
        "cr0_minimum_subnormal_lowT_221_direct_z10",
        freq_221,
        np.asarray(low_t_ref["reference"]["10.0"]["log10_OmegaGW"], dtype=np.float64),
        parameters=low_t_ref["parameters"],
        frozen_dn=float(low_t_ref["frozen_DN_eff"]),
        z_tail=float(low_t_ref["reference"]["10.0"]["z_tail"]),
        rtol=float(low_t_ref["reference"]["10.0"]["rtol"]),
        source_sha256={"grid_input": sha(LOW_T_GRID), "reference_input": sha(LOW_T_REFERENCE)},
    ))

    output = {
        "schema_version": 1,
        "status": "complete_direct_reference_pchip_multigrid_diagnostic",
        "scope": "Independent direct Reference spectra: same-background 469/937 z8 grids and 937 z10, plus a distinct low-temperature 221 z10 background. Compare the existing Reference linear-PCHIP quadrature with log arithmetic for the same positive linear-PCHIP polynomial; do not modify Reference or interpolate spectra between grids.",
        "cases": cases,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "source_sha256": {
            "reference.py": ref_sha,
            "stable_log_pchip.py": sha(ROOT / "scripts" / "stable_log_pchip.py"),
            "diagnostic_script": sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps([{
        "label": case["label"],
        "n": case["frequency_count"],
        "zeros": case["normalized_zero_nodes"],
        "warnings": len(case["original_linear_pchip"]["warnings"]),
        "integral_delta": case["stable_log_arithmetic_same_linear_pchip"]["relative_integral_difference"],
        "zero_panel_log10_fraction": case["stable_log_arithmetic_same_linear_pchip"]["log10_fraction_total_integral_in_zero_adjacent_panels"],
        "midpoint_zeros": case["midpoint_interpolation"]["linear_zero_count"],
        "midpoint_max_dex": case["midpoint_interpolation"]["finite_linear_vs_stable_max_abs_dex"],
    } for case in cases], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
