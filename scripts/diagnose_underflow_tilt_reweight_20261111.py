"""Reweight frozen z8 Reference moments over the saved underflow tilt line.

This is an offline fixed-background diagnostic. It performs no Fast or
Reference propagation and makes no continuous-domain accuracy claim.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import warnings
from pathlib import Path

import numpy as np
import scipy
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
INPUT_NPZ = ROOT / "docs" / "underflow_phase_midpoint_reference_arrays_20261107.npz"
CONVERGENCE_JSON = ROOT / "docs" / "underflow_phase_map_z8_235_469_convergence_20261107.json"
OUTPUT_JSON = ROOT / "docs" / "underflow_tilt_reweight_235_469_20261111.json"
TILTS = (17.30, 17.38, 17.45)
LN10 = math.log(10.0)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_f64(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def main() -> None:
    from stiffgwpy_fast.reference import integrate_spectrum

    convergence = json.loads(CONVERGENCE_JSON.read_text(encoding="utf-8"))
    with np.load(INPUT_NPZ) as data:
        source_frequency = np.asarray(data["frequency_log10_hz"], dtype=np.float64)
        order = np.argsort(source_frequency)
        frequency = source_frequency[order]
        reference_g_unit = np.asarray(
            data["reference_Ogw_unit"] - data["reference_Oj_unit"],
            dtype=np.float64,
        )[order]

    if frequency.size != 469 or reference_g_unit.shape != frequency.shape:
        raise RuntimeError("expected 469 frozen Reference samples")
    if not np.all(np.diff(frequency) > 0.0) or not np.all(reference_g_unit > 0.0):
        raise RuntimeError("frequency order or positive Reference-integrand contract failed")

    coarse_indices = np.arange(0, frequency.size, 2)
    if coarse_indices.size != 235 or coarse_indices[-1] != frequency.size - 1:
        raise RuntimeError("expected the exact nested 235-node support")
    parent_grids = convergence["frequency_grids"]
    if sha256_f64(source_frequency) != parent_grids["frequency_sha256_469"]:
        raise RuntimeError("469 input frequency hash does not match the parent artifact")
    if sha256_f64(source_frequency[::2]) != parent_grids["frequency_sha256_235"]:
        raise RuntimeError("235 input frequency hash does not match the parent artifact")

    results: list[dict[str, float | int | list[str]]] = []
    for tilt in TILTS:
        # r*A_s/f_piv**tilt is constant over frequency and cancels in all
        # relative coarse/fine and panel-cancellation metrics.
        log_weighted_g = np.log(reference_g_unit) + tilt * LN10 * frequency
        weighted = np.exp(log_weighted_g - float(np.max(log_weighted_g)))
        underflowed_nodes = int(np.count_nonzero(weighted == 0.0))

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            coarse_integral, coarse_quad_error, coarse_interp_error = integrate_spectrum(
                frequency[coarse_indices],
                weighted[coarse_indices],
                np.zeros(coarse_indices.size, dtype=np.float64),
            )
            fine_integral, fine_quad_error, fine_interp_error = integrate_spectrum(
                frequency,
                weighted,
                np.zeros(frequency.size, dtype=np.float64),
            )

            coarse_spline = PchipInterpolator(
                frequency[coarse_indices], weighted[coarse_indices], extrapolate=False
            )
            fine_spline = PchipInterpolator(frequency, weighted, extrapolate=False)
            panel_delta = np.empty(234, dtype=np.float64)
            for panel in range(panel_delta.size):
                left, mid, right = frequency[2 * panel : 2 * panel + 3]
                coarse_panel = coarse_spline.integrate(left, right)
                fine_panel = fine_spline.integrate(left, mid) + fine_spline.integrate(mid, right)
                panel_delta[panel] = LN10 * (fine_panel - coarse_panel)

        signed_delta = float(fine_integral - coarse_integral)
        absolute_panel_delta = float(np.sum(np.abs(panel_delta)))
        results.append(
            {
                "n_t": tilt,
                "normalized_integral_235": float(coarse_integral),
                "normalized_integral_469": float(fine_integral),
                "fine_minus_coarse_relative": float(signed_delta / fine_integral),
                "sum_absolute_panel_delta_relative": float(absolute_panel_delta / fine_integral),
                "panel_cancellation_ratio": float(
                    absolute_panel_delta / abs(signed_delta) if signed_delta else math.inf
                ),
                "signed_panel_delta_relative": float(np.sum(panel_delta) / fine_integral),
                "quadrature_error_235_normalized": (
                    None if coarse_quad_error is None else float(coarse_quad_error)
                ),
                "quadrature_error_469_normalized": (
                    None if fine_quad_error is None else float(fine_quad_error)
                ),
                "interpolation_error_235_normalized": (
                    None if coarse_interp_error is None else float(coarse_interp_error)
                ),
                "interpolation_error_469_normalized": (
                    None if fine_interp_error is None else float(fine_interp_error)
                ),
                "weighted_nodes_underflowed_to_zero": underflowed_nodes,
                "scipy_warning_messages": sorted({str(w.message) for w in caught}),
            }
        )

    baseline = next(row for row in results if row["n_t"] == 17.38)
    expected_baseline = float(
        convergence["reference_frequency_convergence"]["converged_roots"][
            "relative_469_minus_235"
        ]
    )
    output = {
        "schema_version": 1,
        "status": "complete_offline_fixed_background_reweighting",
        "experiment": "z8 Reference 235/469 PCHIP tilt reweighting",
        "parameters_inherited_from_convergence": convergence["parameters"],
        "background_scope": {
            "fixed_DN_eff": convergence["reference_configuration"]["frozen_DN_eff"],
            "tail_convention": "z_tail=8",
            "only_parameter_varied": "n_t",
            "tilts": list(TILTS),
        },
        "frequency_grid_hashes": {
            "469_source_order_sha256": sha256_f64(source_frequency),
            "235_source_order_sha256": sha256_f64(source_frequency[::2]),
        },
        "method": {
            "input": "saved unit Reference Ogw-Oj arrays; no new propagation",
            "coarse_grid": "every other point of the nested 469 grid (235 nodes)",
            "fine_grid": "all 469 saved points",
            "integral": "Reference integrate_spectrum (PCHIP in log10 f, dln f)",
            "panel_delta": "sum fine two-subpanel PCHIP integrals minus coarse PCHIP integral per parent panel",
            "normalization": "subtract max log integrand; positive constants r*A_s/f_piv**n_t cancel relative metrics",
        },
        "reproduction_check": {
            "n_t": 17.38,
            "saved_nested_root_relative_469_minus_235": expected_baseline,
            "reweighted_fixed_background_relative_469_minus_235": float(
                baseline["fine_minus_coarse_relative"]
            ),
            "absolute_difference": abs(
                float(baseline["fine_minus_coarse_relative"]) - expected_baseline
            ),
        },
        "results": results,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "interpretation": (
            "The fixed-background 235/469 signed integral delta stays near 1.61e-6 relative "
            "over this saved one-dimensional tilt line, while the sum of absolute panel deltas "
            "is about 10.9 times larger. Global normalization underflows 24-25 negligible-weight "
            "nodes to zero and SciPy reports an overflow warning while forming PCHIP slopes. "
            "This is a narrow reweighting robustness screen only; "
            "it is not a Fast error bound, root certificate, or multidimensional region result."
        ),
        "source_hashes": {
            "script_sha256": sha256(Path(__file__)),
            "input_npz_sha256": sha256(INPUT_NPZ),
            "convergence_json_sha256": sha256(CONVERGENCE_JSON),
            "reference_py_sha256": sha256(ROOT / "stiffgwpy_fast" / "reference.py"),
        },
    }
    OUTPUT_JSON.write_text(
        json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
