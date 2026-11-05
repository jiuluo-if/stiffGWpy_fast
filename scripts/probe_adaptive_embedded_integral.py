"""Probe an embedded global PCHIP integral delta on a fixed Fast background.

This diagnostic takes the archived adaptive Fast frequency grid, evaluates
its actual log-frequency midpoints with the same Fast propagation kernel, and
compares the coarse/refined integral delta with archived independent-grid
changes. It does not change the production adaptive rule or call Reference.
"""

from __future__ import annotations

import hashlib
import json
import os
import statistics
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, telemetry  # noqa: E402

apply_environment()

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

SOURCE_NAME = os.environ.get(
    "ADAPTIVE_EMBEDDED_SOURCE",
    "fast_high_blue_target1e7_prufer_20261105.json",
)
SOURCE = ROOT / "docs" / SOURCE_NAME
OUTPUT_NAME = (
    "adaptive_embedded_integral_fast_midpoints_20261105.json"
    if SOURCE_NAME == "fast_high_blue_target1e7_prufer_20261105.json"
    else "adaptive_cr1_prufer_embedded_midpoints_20261105.json"
)
OUTPUT = ROOT / "docs" / OUTPUT_NAME
REPEATS = 7
_SOURCE_DATA = json.loads(SOURCE.read_text(encoding="utf-8"))
SOURCE_Z_TAIL = float(
    _SOURCE_DATA["fast"].get(
        "z_tail", _SOURCE_DATA.get("reference", {}).get("z_tail", 5.0)
    )
)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True,
                            text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def _prepared(model: LCDM_SG, freqs: np.ndarray, dn_eff: float):
    config = FS.resolve_config(
        "fast", z_tail=SOURCE_Z_TAIL, threads=2
    )
    h = float(config.h)
    if not np.isclose(float(model.cosmo_param["DN_eff"]), dn_eff, rtol=0.0, atol=1e-14):
        raise RuntimeError("model background does not match the requested fixed DN_eff")
    Nv = np.asarray(model.Nv, dtype=np.float64)
    Sv, _, j0s, z0s, fp_minus = FS.prep_frequency_only(model, Nv, freqs)
    from stiffgwpy_fast.exact_background import exact_phi_s2_split

    Phi_grid, Phi_mid, S2, S2inv, kink_index, kink_fraction, phi_re = (
        exact_phi_s2_split(model, Nv, float(dn_eff), sigma_nodes=model.sigma)
    )
    col_step = int(config.col_step)
    idx_out = np.unique(np.append(np.arange(0, len(Nv), col_step), len(Nv) - 1))
    P_t = model.derived_param["A_t"] * np.power(
        (10.0 ** freqs) / gp.f_piv, model.derived_param["nt"]
    )
    ev_minus = np.exp(-Nv)
    fp_freq = np.power(10.0, freqs)
    shape = (len(freqs), len(idx_out))
    outputs = [np.zeros(shape, dtype=np.float64) for _ in range(3)]
    args = (
        Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
        fp_minus, fp_freq, 1, len(idx_out), col_step, h,
        float(config.z_tail), *outputs, None, model.sigma,
        float(config.phase_max), np.full(len(freqs), -1.0),
        kink_index, kink_fraction, phi_re,
        FS._uniform_step_terminal_width(Nv, h),
    )
    return args, outputs[0][:, -1] - outputs[1][:, -1]


def _evaluate(args, outputs):
    for output in outputs:
        output.fill(0.0)
    FS.solve_kernel(*args)
    return outputs[0][:, -1] - outputs[1][:, -1]


def _dn_gw(freqs: np.ndarray, spectrum: np.ndarray, model: LCDM_SG) -> float:
    g2 = FS.integrate_frequency_pchip(freqs, spectrum) * FS.ln10
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    return float(gp.Neff0 * g2 / omega_nu)


def main() -> None:
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    source = _SOURCE_DATA
    params = source["parameters"]
    dn_eff = float(source.get("fixed_DN_eff", source["fast"]["DN_eff_final"]))
    base_freqs = np.asarray(source["fast"]["frequencies_log10_hz"], dtype=np.float64)
    x = np.sort(base_freqs)
    mids = 0.5 * (x[:-1] + x[1:])
    refined_freqs = np.sort(np.concatenate((x, mids)))[::-1].copy()
    base_freqs = np.sort(base_freqs)[::-1].copy()
    model = LCDM_SG(**params)
    started = time.perf_counter()
    fast_result = model.SGWB_iter(
        engine="fast",
        fallback=False,
        tol=1e-6,
        z_tail=SOURCE_Z_TAIL,
        freq_res=1.0,
        threads=2,
        accuracy_mode="fast",
        physical_guard=False,
        freq_grid="adaptive",
        freq_grid_target=float(source["fast"].get("target_dex", 1e-4)),
        freq_grid_max_points=int(source["fast"].get("max_points", 1500)),
        sigma_exact=bool(source["fast"].get("sigma_exact", True)),
    )
    outer_replay_seconds = time.perf_counter() - started
    if fast_result is not model:
        raise RuntimeError(f"Fast baseline replay failed: {model.fast_failure_reason}")
    if not np.array_equal(np.asarray(model.f), base_freqs):
        raise RuntimeError("Fast baseline replay did not reproduce the archived frequency grid")
    if not np.isclose(float(model.cosmo_param["DN_eff"]), dn_eff, rtol=0.0, atol=2e-12):
        raise RuntimeError("Fast baseline replay did not reproduce the archived fixed-point background")
    base_args, _ = _prepared(model, base_freqs, dn_eff)
    refined_args, _ = _prepared(model, refined_freqs, dn_eff)
    base_outputs = [base_args[16], base_args[17], base_args[18]]
    refined_outputs = [refined_args[16], refined_args[17], refined_args[18]]

    # Warm both compiled signatures before collecting paired kernel timings.
    base_spectrum = _evaluate(base_args, base_outputs)
    refined_spectrum = _evaluate(refined_args, refined_outputs)
    base_dn = _dn_gw(base_freqs, base_spectrum, model)
    refined_dn = _dn_gw(refined_freqs, refined_spectrum, model)
    archived_spectrum = np.asarray(source["fast"]["omega_gw"], dtype=np.float64)
    if archived_spectrum.size != base_spectrum.size:
        raise RuntimeError("archived Fast spectrum length does not match archived frequency grid")
    spectrum_abs_dex = np.abs(np.log10(base_spectrum / archived_spectrum))
    ref_maps = source["reference_maps"]
    adaptive_ref_key = (
        "adaptive_exact_fast_grid"
        if "adaptive_exact_fast_grid" in ref_maps
        else "adaptive_fast_grid"
    )
    dense_ref_key = (
        "native_res4"
        if "native_res4" in ref_maps
        else "native_res2_reference"
    )
    archived_ref = ref_maps[adaptive_ref_key]["DN_gw"]
    archived_dense = (
        ref_maps[dense_ref_key]["DN_gw"]
        if dense_ref_key in ref_maps
        else source[dense_ref_key]["DN_gw"]
    )
    _, coarse_panel_integrals = FS.pchip_integral_breakdown(base_freqs, base_spectrum)
    _, refined_panel_integrals = FS.pchip_integral_breakdown(refined_freqs, refined_spectrum)
    if refined_panel_integrals.size != 2 * coarse_panel_integrals.size:
        raise RuntimeError("midpoint grid does not split every base panel exactly once")
    panel_delta = refined_panel_integrals.reshape(-1, 2).sum(axis=1) - coarse_panel_integrals
    if not np.isclose(
        float(np.sum(panel_delta)),
        float(np.sum(refined_panel_integrals) - np.sum(coarse_panel_integrals)),
        rtol=1e-10,
        atol=1e-18,
    ):
        raise RuntimeError("panel integral deltas do not sum to the global integral delta")
    panel_order = np.argsort(np.abs(panel_delta))[::-1]
    panel_selection = []
    for fraction in (0.25, 0.5, 0.75):
        count = max(1, int(np.ceil(fraction * mids.size)))
        selected_mids = mids[panel_order[:count]]
        candidate_freqs = np.sort(np.concatenate((x, selected_mids)))[::-1].copy()
        candidate_args, _ = _prepared(model, candidate_freqs, dn_eff)
        candidate_outputs = [candidate_args[16], candidate_args[17], candidate_args[18]]
        candidate_spectrum = _evaluate(candidate_args, candidate_outputs)
        candidate_dn = _dn_gw(candidate_freqs, candidate_spectrum, model)
        candidate_times = []
        for _ in range(3):
            start = time.perf_counter()
            _evaluate(candidate_args, candidate_outputs)
            candidate_times.append(time.perf_counter() - start)
        panel_selection.append({
            "selected_fraction": fraction,
            "selected_midpoint_count": int(count),
            "frequency_count": int(candidate_freqs.size),
            "DN_gw": candidate_dn,
            "relative_delta_vs_base": abs(candidate_dn - base_dn) / max(abs(candidate_dn), abs(base_dn)),
            "relative_delta_vs_all_midpoints": abs(candidate_dn - refined_dn) / max(abs(candidate_dn), abs(refined_dn)),
            "relative_delta_vs_archived_dense_reference": abs(candidate_dn - float(archived_dense)) / max(abs(candidate_dn), abs(float(archived_dense))),
            "kernel_median_s": statistics.median(candidate_times),
            "kernel_p95_s": float(np.percentile(candidate_times, 95)),
        })
    base_order = np.argsort(base_freqs)
    x_base = base_freqs[base_order]
    y_base = np.log10(base_spectrum[base_order])
    curvature = np.abs(np.asarray(PchipInterpolator(x_base, y_base).derivative(2)(x_base)))
    widths = np.diff(x_base)
    omega_base = base_spectrum[base_order]
    panel_weights = 0.5 * (omega_base[:-1] + omega_base[1:]) * widths
    curvature_score = curvature[:-1] * widths * widths / 8.0 * panel_weights
    curvature_order = np.argsort(curvature_score)[::-1]
    curvature_selection = []
    for fraction in (0.25, 0.5, 0.75):
        count = max(1, int(np.ceil(fraction * mids.size)))
        selected_indices = curvature_order[:count]
        selected_mids = mids[selected_indices]
        candidate_freqs = np.sort(np.concatenate((x, selected_mids)))[::-1].copy()
        candidate_args, _ = _prepared(model, candidate_freqs, dn_eff)
        candidate_outputs = [candidate_args[16], candidate_args[17], candidate_args[18]]
        candidate_spectrum = _evaluate(candidate_args, candidate_outputs)
        candidate_dn = _dn_gw(candidate_freqs, candidate_spectrum, model)
        candidate_times = []
        for _ in range(3):
            start = time.perf_counter()
            _evaluate(candidate_args, candidate_outputs)
            candidate_times.append(time.perf_counter() - start)
        true_abs_delta_share = float(
            np.sum(np.abs(panel_delta[selected_indices])) / np.sum(np.abs(panel_delta))
        )
        curvature_selection.append({
            "selected_fraction": fraction,
            "selected_midpoint_count": int(count),
            "frequency_count": int(candidate_freqs.size),
            "absolute_panel_delta_share_captured": true_abs_delta_share,
            "DN_gw": candidate_dn,
            "relative_delta_vs_base": abs(candidate_dn - base_dn) / max(abs(candidate_dn), abs(base_dn)),
            "relative_delta_vs_all_midpoints": abs(candidate_dn - refined_dn) / max(abs(candidate_dn), abs(refined_dn)),
            "relative_delta_vs_archived_dense_reference": abs(candidate_dn - float(archived_dense)) / max(abs(candidate_dn), abs(float(archived_dense))),
            "kernel_median_s": statistics.median(candidate_times),
            "kernel_p95_s": float(np.percentile(candidate_times, 95)),
        })
    base_times, refined_times = [], []
    for _ in range(REPEATS):
        start = time.perf_counter()
        _evaluate(base_args, base_outputs)
        base_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _evaluate(refined_args, refined_outputs)
        refined_times.append(time.perf_counter() - start)

    record = {
        "schema_version": 1,
        "experiment": "fixed_background_fast_actual_midpoint_embedded_integral_probe",
        "source_commit": _git("rev-parse", "HEAD"),
        "created_at": datetime.now().astimezone().isoformat(),
        "parameters": params,
        "fixed_DN_eff": dn_eff,
        "fast_config": {
            "h": float(FS.resolve_config("fast", z_tail=SOURCE_Z_TAIL, threads=2).h),
            "z_tail": SOURCE_Z_TAIL,
            "col_step": int(FS.resolve_config("fast", z_tail=SOURCE_Z_TAIL, threads=2).col_step),
            "phase_max": float(FS.resolve_config("fast", z_tail=SOURCE_Z_TAIL, threads=2).phase_max),
            "sigma_exact": bool(source["fast"].get("sigma_exact", True)),
        },
        "fast_outer_replay_seconds": outer_replay_seconds,
        "source_artifact": SOURCE.name,
        "adaptive_reference_key": adaptive_ref_key,
        "dense_reference_key": dense_ref_key,
        "source_artifact_sha256": _sha(SOURCE.read_bytes()),
        "base_frequency_count": int(base_freqs.size),
        "actual_midpoint_count": int(mids.size),
        "refined_frequency_count": int(refined_freqs.size),
        "base_DN_gw": base_dn,
        "refined_DN_gw": refined_dn,
        "archived_fast_DN_gw": float(source["fast"]["DN_gw"]),
        "base_relative_delta_vs_archived_fast": abs(base_dn - float(source["fast"]["DN_gw"])) / max(abs(base_dn), abs(float(source["fast"]["DN_gw"]))),
        "base_spectrum_max_abs_dex_vs_archived_fast": float(np.max(spectrum_abs_dex)),
        "base_spectrum_p95_abs_dex_vs_archived_fast": float(np.quantile(spectrum_abs_dex, 0.95)),
        "embedded_global_relative_delta": abs(refined_dn - base_dn) / max(abs(refined_dn), abs(base_dn)),
        "archived_reference_adaptive_DN_gw": float(archived_ref),
        "archived_dense_reference_DN_gw": float(archived_dense),
        "archived_reference_grid_relative_delta": abs(archived_dense - archived_ref) / max(abs(archived_dense), abs(archived_ref)),
        "embedded_to_archived_grid_delta_ratio": (
            abs(refined_dn - base_dn) / max(abs(refined_dn), abs(base_dn))
            / (abs(archived_dense - archived_ref) / max(abs(archived_dense), abs(archived_ref)))
        ),
        "panel_integral_delta": {
            "sum_signed_relative": abs(float(np.sum(panel_delta)) * FS.ln10 * gp.Neff0 / (gp.Omega_nh2 / model.derived_param["h"] ** 2)) / max(abs(base_dn), abs(refined_dn)),
            "sum_absolute_relative": float(np.sum(np.abs(panel_delta)) * FS.ln10 * gp.Neff0 / (gp.Omega_nh2 / model.derived_param["h"] ** 2)) / max(abs(base_dn), abs(refined_dn)),
            "largest_panel_fraction_of_absolute_delta": float(np.max(np.abs(panel_delta)) / np.sum(np.abs(panel_delta))),
            "selection_results": panel_selection,
            "curvature_ranked_selection_results": curvature_selection,
            "limitation": "Panels are ranked using all midpoint values; this is an upper-bound feasibility screen, not the work of an online selector.",
        },
        "kernel_seconds": {
            "base_median": statistics.median(base_times),
            "base_p95": float(np.percentile(base_times, 95)),
            "refined_median": statistics.median(refined_times),
            "refined_p95": float(np.percentile(refined_times, 95)),
            "median_ratio_refined_over_base": statistics.median(refined_times) / statistics.median(base_times),
            "paired_repeats": REPEATS,
            "warmup_excluded": True,
        },
        "resources": telemetry(threads=2),
        "hashes": {
            "script_py": _sha(Path(__file__).read_bytes()),
            "fast_sgwb_py": _sha((ROOT / "stiffgwpy_fast" / "fast_sgwb.py").read_bytes()),
            "exact_background_py": _sha((ROOT / "stiffgwpy_fast" / "exact_background.py").read_bytes()),
        },
        "limitations": [
            "Fixed-background, one-parameter-point probe only; not a root-level accuracy result.",
            "The midpoint values are actual Fast propagation outputs; no interpolated midpoint values are used.",
            "The archived independent reference grids are nonnested and provide a convergence comparison, not a rigorous bound.",
            "Kernel-only timing excludes background preparation and full adaptive-loop orchestration.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
