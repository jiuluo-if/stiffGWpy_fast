"""Reconstruct paired full-grid Reference metrics from per-mode records.

中文：从逐模配对日志重建当前 Cartesian 与 Prüfer 的频谱差、积分映射及根残差。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _stats(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    return {
        "median": float(np.percentile(values, 50)),
        "p95": float(np.percentile(values, 95)),
        "max": float(np.max(values)),
    }


def _dn_gw(freqs: np.ndarray, log_omega: np.ndarray, model: LCDM_SG):
    omega = np.power(10.0, log_omega)
    g2, quad_error, interp_error = REF.integrate_spectrum(freqs, omega, np.zeros_like(omega))
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    value = gp.Neff0 * g2 / omega_nu
    return float(value), quad_error, interp_error


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--paired", type=Path, default=ROOT / "docs" / "reference_prufer_cartesian_paired_fullgrid_20261104.json"
    )
    parser.add_argument(
        "--archive", type=Path, default=ROOT / "docs" / "res4_reference_secant_root_tail10_20261104.npz"
    )
    parser.add_argument(
        "--reference-map", type=Path, default=ROOT / "docs" / "res4_reference_secant_root_tail10_20261104.json"
    )
    parser.add_argument("--secant", type=Path, default=ROOT / "docs" / "res4_reference_map_sensitivity_20261104.json")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "docs" / "reference_prufer_cartesian_paired_fullgrid_analysis_20261104.json",
    )
    args = parser.parse_args(argv)

    paired = json.loads(args.paired.read_text(encoding="utf-8"))
    reference_map = json.loads(args.reference_map.read_text(encoding="utf-8"))
    secant = json.loads(args.secant.read_text(encoding="utf-8"))
    with np.load(args.archive) as archive:
        freqs = np.asarray(archive["freqs"], dtype=float)
        archived_log = np.asarray(archive["reference_log"], dtype=float)

    rows = paired["rows"]
    current_reference_hash = _sha256(ROOT / "stiffgwpy_fast" / "reference.py")
    if current_reference_hash != paired["reference_module_sha256"]:
        raise ValueError(
            "current Reference module hash differs from the paired run; "
            "refusing to reconstruct with mismatched integration code"
        )
    if len(rows) != freqs.size or len(rows) != 260:
        raise ValueError("paired rows must cover the complete 260-mode archive")
    if not np.array_equal(np.asarray([row["freq_log10_hz"] for row in rows]), freqs):
        raise ValueError("paired frequency order differs from the archived grid")
    if not np.array_equal(np.asarray([row["index"] for row in rows]), np.arange(freqs.size)):
        raise ValueError("paired mode indices are not complete and ordered")

    prufer_log = np.asarray([row["candidate_log10_omega"] for row in rows], dtype=float)
    cartesian_log = np.asarray([row["cartesian_log10_omega"] for row in rows], dtype=float)
    candidate_times = np.asarray([row["candidate_runtime_s"] for row in rows], dtype=float)
    cartesian_times = np.asarray([row["cartesian_runtime_s"] for row in rows], dtype=float)
    paired_dex = np.abs(prufer_log - cartesian_log)
    archive_cart_dex = np.abs(cartesian_log - archived_log)
    archive_prufer_dex = np.abs(prufer_log - archived_log)
    speedups = cartesian_times / candidate_times
    if not all(
        np.all(np.isfinite(values))
        for values in (prufer_log, cartesian_log, archived_log, candidate_times, cartesian_times)
    ):
        raise FloatingPointError("paired full-grid rows contain non-finite values")
    if any(np.any(values <= -40.0) for values in (archived_log, prufer_log, cartesian_log)):
        raise ValueError("the benchmark's -40 log floor would prevent exact reconstruction")
    if np.any(candidate_times <= 0.0) or np.any(cartesian_times <= 0.0):
        raise ValueError("paired per-mode runtimes must be positive")

    params = paired["parameters"]
    model = LCDM_SG(
        **{key: params[key] for key in ("r", "cr", "n_t", "T_re", "kappa10", "DN_re")}, DN_eff=params["DN_eff"]
    )
    candidate_dn, candidate_quad, candidate_interp = _dn_gw(freqs, prufer_log, model)
    cartesian_dn, cartesian_quad, cartesian_interp = _dn_gw(freqs, cartesian_log, model)
    candidate_recorded = float(paired["candidate_DN_gw"])
    map_reconstruction_relative = abs(candidate_dn - candidate_recorded) / abs(candidate_recorded)
    if map_reconstruction_relative > 5e-10:
        raise ArithmeticError(
            "log-spectrum map reconstruction differs from benchmark map by %.3g" % map_reconstruction_relative
        )

    root_input = float(params["DN_eff"])
    cartesian_residual = cartesian_dn - root_input
    candidate_residual = candidate_dn - root_input
    slope = float(secant["secant_slope_dDNgw_dDNeff"])
    root_shift_estimate = (candidate_dn - cartesian_dn) / (1.0 - slope)
    payload = {
        "schema_version": 1,
        "analysis": "paired_current_cartesian_vs_prufer_full_260_mode_res4",
        "source_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "analysis_script_sha256": _sha256(Path(__file__).resolve()),
        "reference_module_sha256": paired["reference_module_sha256"],
        "parameters": params,
        "z_tail": paired["z_tail"],
        "candidate_rtol": paired["candidate_rtol"],
        "cartesian_rtol": paired["cartesian_rtol"],
        "frequency_count": int(freqs.size),
        "successful_paired_modes": int(len(rows)),
        "candidate_used_tail_modes": int(sum(bool(row["used_tail"]) for row in rows)),
        "paired_same_input_spectrum_abs_dex": _stats(paired_dex),
        "prufer_vs_archived_cartesian_abs_dex": _stats(archive_prufer_dex),
        "current_cartesian_vs_archived_cartesian_abs_dex": _stats(archive_cart_dex),
        "per_mode_speedup_cartesian_over_prufer": _stats(speedups),
        "timing": {
            "prufer_sum_s": float(np.sum(candidate_times)),
            "cartesian_sum_s": float(np.sum(cartesian_times)),
            "prufer_mode_runtime_s": _stats(candidate_times),
            "cartesian_mode_runtime_s": _stats(cartesian_times),
            "prufer_rhs_evaluations_total": int(sum(row["candidate_nfev"] for row in rows)),
            "cartesian_rhs_evaluations_total": int(sum(row["cartesian_nfev"] for row in rows)),
            "aggregate_speedup": float(np.sum(cartesian_times) / np.sum(candidate_times)),
            "prufer_modes_per_second_serial": float(freqs.size / np.sum(candidate_times)),
            "cartesian_modes_per_second_serial": float(freqs.size / np.sum(cartesian_times)),
            "jit_warmup_s": float(paired["spline_jit_warmup_seconds"]),
            "aggregate_speedup_including_jit": float(
                np.sum(cartesian_times) / (np.sum(candidate_times) + paired["spline_jit_warmup_seconds"])
            ),
        },
        "integral_map": {
            "root_input_dn_eff": root_input,
            "candidate_dn_gw": candidate_dn,
            "candidate_recorded_dn_gw": candidate_recorded,
            "candidate_map_reconstruction_relative_error": (map_reconstruction_relative),
            "current_cartesian_dn_gw": cartesian_dn,
            "archived_cartesian_dn_gw": float(reference_map["reference_dn_map"]),
            "candidate_vs_current_cartesian_map_relative_difference": (
                abs(candidate_dn - cartesian_dn) / abs(cartesian_dn)
            ),
            "current_cartesian_vs_archived_map_relative_difference": (
                abs(cartesian_dn - reference_map["reference_dn_map"]) / abs(reference_map["reference_dn_map"])
            ),
            "candidate_fixed_point_residual_relative": (candidate_residual / root_input),
            "current_cartesian_fixed_point_residual_relative": (cartesian_residual / root_input),
            "candidate_map_quadrature_error": candidate_quad,
            "candidate_map_interpolation_error": candidate_interp,
            "current_cartesian_map_quadrature_error": cartesian_quad,
            "current_cartesian_map_interpolation_error": cartesian_interp,
            "archived_interval_secant_slope": slope,
            "candidate_cartesian_root_shift_first_order_estimate": (root_shift_estimate),
            "root_shift_method_limit": (
                "uses an archived interval secant average, not a local derivative or a rigorous root enclosure"
            ),
        },
        "input_provenance": {
            "paired_json": str(args.paired.resolve()),
            "paired_json_sha256": _sha256(args.paired),
            "frequency_spectrum_npz": str(args.archive.resolve()),
            "frequency_spectrum_npz_sha256": _sha256(args.archive),
            "reference_map_json": str(args.reference_map.resolve()),
            "reference_map_json_sha256": _sha256(args.reference_map),
            "secant_analysis_json": str(args.secant.resolve()),
            "secant_analysis_json_sha256": _sha256(args.secant),
        },
        "semantics": (
            "post-processing of 260 paired per-mode DOP853 Reference runs; "
            "log-spectrum inversion is exact for this archive because all nodes "
            "are above the benchmark's -40 floor"
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
