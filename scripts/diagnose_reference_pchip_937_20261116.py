"""Diagnose Reference PCHIP underflow arithmetic on the direct 937-node data."""

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

GRID_INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
REFERENCE_INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_z10_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "reference_pchip_937_diagnostic_20261116.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    grid_data = json.loads(GRID_INPUT.read_text(encoding="utf-8"))
    reference = json.loads(REFERENCE_INPUT.read_text(encoding="utf-8"))
    freq = np.asarray(grid_data["frequency_grids"]["937"]["log10_hz"], dtype=np.float64)
    log10_y = np.asarray(reference["reference"]["log10_OmegaGW"], dtype=np.float64)
    grid_hash = hashlib.sha256(np.asarray(freq, dtype="<f8").tobytes()).hexdigest()
    if freq.size != 937 or grid_hash != grid_data["frequency_grids"]["937"]["sha256_f64le"]:
        raise RuntimeError("canonical direct 937-point frequency grid hash mismatch")
    if grid_hash != reference["frequency_sha256_f64le"]:
        raise RuntimeError("direct Reference grid does not match the canonical 937-point grid")
    if log10_y.shape != freq.shape or not np.isfinite(log10_y).all():
        raise RuntimeError("direct Reference log spectrum is incomplete or non-finite")

    order = np.argsort(freq)
    x = freq[order]
    log_y = log10_y[order] * math.log(10.0)
    shift = float(np.max(log_y))
    with np.errstate(under="ignore"):
        y = np.exp(log_y - shift)

    with warnings.catch_warnings(record=True) as original_warnings:
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
            raise FloatingPointError(f"stable PCHIP panel {i} has non-positive area")
        panel_logs.append(log_area)
    panel_logs = np.asarray(panel_logs, dtype=np.float64)
    zero_adjacent_panels = (y[:-1] == 0.0) | (y[1:] == 0.0)
    zero_panel_log_fraction = float(
        np.logaddexp.reduce(panel_logs[zero_adjacent_panels])
        - np.logaddexp.reduce(panel_logs)
    )

    midpoint = (x[:-1] + x[1:]) / 2.0
    with warnings.catch_warnings(record=True) as midpoint_warnings:
        warnings.simplefilter("always", RuntimeWarning)
        linear_midpoint = PchipInterpolator(x, y)(midpoint)
    stable_log_midpoint = stable.evaluate_log(midpoint)
    representable = linear_midpoint > 0.0
    midpoint_error = np.full(midpoint.shape, np.nan, dtype=np.float64)
    midpoint_error[representable] = (
        np.log10(linear_midpoint[representable])
        - (stable_log_midpoint[representable] - shift) / math.log(10.0)
    )
    worst_midpoint = int(np.nanargmax(np.abs(midpoint_error))) if representable.any() else None

    output = {
        "schema_version": 1,
        "status": "complete_reference_937_linear_pchip_warning_diagnostic",
        "scope": "Reuse direct Reference z10 log nodes on the exact canonical 937-point grid; compare original normalized linear PCHIP quadrature with log-arithmetic evaluation/integration of the positive linear-PCHIP Hermite polynomial. Reference production source is unchanged.",
        "frequency_count": int(freq.size),
        "frequency_sha256_f64le": grid_hash,
        "frozen_DN_eff": reference["frozen_DN_eff"],
        "source_values": {
            "minimum_log10_OmegaGW": float(np.min(log10_y)),
            "maximum_log10_OmegaGW": float(np.max(log10_y)),
            "finite_count": int(np.isfinite(log10_y).sum()),
        },
        "normalization": {
            "shift_natural_log": shift,
            "normalized_zero_nodes": int(np.count_nonzero(y == 0.0)),
            "normalized_positive_subnormal_nodes": int(np.count_nonzero(
                (y > 0.0) & (y < np.finfo(np.float64).tiny)
            )),
        },
        "original_linear_pchip": {
            "runtime_warnings": [str(item.message) for item in original_warnings],
            "log_integral": original_log_integral,
            "quadrature_error_relative": float(qerr / g2),
            "half_density_interpolation_error_relative": float(ierr / g2),
        },
        "stable_log_arithmetic_same_pchip": {
            "runtime_warnings": [],
            "log_integral": stable_log_integral,
            "relative_integral_difference": math.expm1(stable_log_integral - original_log_integral),
            "interpolant_definition": "Same positive linear PCHIP Hermite polynomial from log-domain secants and derivatives; not PCHIP applied to log(y).",
            "panels_touching_a_normalized_zero_node": int(np.count_nonzero(zero_adjacent_panels)),
            "log10_fraction_of_total_integral_in_those_panels": zero_panel_log_fraction / math.log(10.0),
        },
        "midpoint_interpolation": {
            "point_count": int(midpoint.size),
            "linear_zero_count": int(np.count_nonzero(~representable)),
            "runtime_warnings": [str(item.message) for item in midpoint_warnings],
            "finite_linear_vs_stable_log_pchip_max_abs_dex": float(np.nanmax(np.abs(midpoint_error))),
            "finite_linear_vs_stable_log_pchip_p95_abs_dex": float(np.nanpercentile(np.abs(midpoint_error[representable]), 95.0)),
            "worst_finite_midpoint": None if worst_midpoint is None else {
                "index": worst_midpoint,
                "log10_frequency": float(midpoint[worst_midpoint]),
                "difference_dex": float(midpoint_error[worst_midpoint]),
                "left_node_log10_signal": float(log10_y[order][worst_midpoint]),
                "right_node_log10_signal": float(log10_y[order][worst_midpoint + 1]),
            },
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "source_sha256": {
            "reference.py": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "stable_log_pchip.py": sha(ROOT / "scripts" / "stable_log_pchip.py"),
            "diagnostic_script": sha(Path(__file__)),
            "grid_input": sha(GRID_INPUT),
            "reference_input": sha(REFERENCE_INPUT),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "zero_nodes": output["normalization"]["normalized_zero_nodes"],
        "original_warning_count": len(output["original_linear_pchip"]["runtime_warnings"]),
        "integral_relative_difference": output["stable_log_arithmetic_same_pchip"]["relative_integral_difference"],
        "midpoint_zero_count": output["midpoint_interpolation"]["linear_zero_count"],
        "midpoint_max_abs_dex": output["midpoint_interpolation"]["finite_linear_vs_stable_log_pchip_max_abs_dex"],
        "midpoint_p95_abs_dex": output["midpoint_interpolation"]["finite_linear_vs_stable_log_pchip_p95_abs_dex"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
