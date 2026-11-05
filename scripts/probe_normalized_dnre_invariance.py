"""Probe reheating-aligned background invariance under ``DN_re`` changes.

This diagnostic separates a background-coordinate identity from its
consequences for the complete Fast spectrum and self-consistent root. It does
not compare either model to Reference and is not an accuracy certification.

中文：比较不同 ``DN_re`` 下的归一化背景段，并单独报告完整 Fast 频谱/自洽根的变化。
该探针不与 Reference 比较，也不构成精度认证。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import time
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

import numba  # noqa: E402
import numpy as np  # noqa: E402

from stiffgwpy_fast import LCDM_SG  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.exact_background import (  # noqa: E402
    H2_vec,
    build_transition_grid,
    exact_phi_s2_grid,
    sigma_vec,
)

ROOT = Path(__file__).resolve().parents[1]
DN_RE_VALUES = (8.0, 12.0)
S_SAMPLES = np.array([0.05, 0.25, 0.5, 0.75, 0.95], dtype=np.float64)
F_OFFSETS = np.array([-0.75, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75], dtype=np.float64)
FIXED_MAP_OFFSETS = np.array([-8.0, -6.0, -4.0, -3.0, -2.0, 0.75], dtype=np.float64)
PARAMS = {
    "r": 1e-2,
    "cr": 0.0,
    "n_t": 0.0,
    "T_re": 2e3,
    "kappa10": 1e-2,
}
CONFIG = FS.FastSolverConfig(
    h=0.005,
    col_step=8,
    z_tail=5.0,
    phase_max=0.25,
    freq_grid="goal",
    threads=1,
    kink_split=True,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json_safe(value):
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def _run_fixed_map(dn_re: float, probe_freqs: np.ndarray) -> dict:
    """Evaluate the complete Fast mode kernel at the same fixed DN_eff input."""
    model = LCDM_SG(**PARAMS, DN_re=dn_re)
    Nv, sigma, f_hor, _ = build_transition_grid(model, CONFIG.h)
    from stiffgwpy_fast.exact_background import exact_phi_s2_split

    Phi, Phi_mid, S2, S2inv, kink_index, kink_fraction, phi_re = exact_phi_s2_split(
        model, Nv, model.cosmo_param["DN_eff"], sigma_nodes=sigma
    )
    _, _, j0s, z0s, fp_minus = FS.prep_frequency_only(model, Nv, probe_freqs)
    derived = model.derived_param
    P_t = derived["A_t"] * np.power(
        (10.0 ** probe_freqs) / gp.f_piv, derived["nt"]
    )
    ev_minus = np.exp(-Nv)
    fp_freq = np.power(10.0, probe_freqs)
    col_step = CONFIG.col_step
    idx_out = np.unique(np.append(np.arange(0, len(Nv), col_step), len(Nv) - 1))
    outputs = [np.zeros((probe_freqs.size, idx_out.size), dtype=np.float64) for _ in range(3)]
    handoff_eps = np.full(probe_freqs.size, -1.0, dtype=np.float64)
    threads_before = FS.get_num_threads()
    if threads_before != CONFIG.threads:
        numba.set_num_threads(CONFIG.threads)
    try:
        FS.solve_kernel(
            Nv, Phi, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus, fp_minus,
            fp_freq, 1, len(idx_out), col_step, CONFIG.h, CONFIG.z_tail,
            *outputs, np.diff(Nv), sigma, CONFIG.phase_max, handoff_eps,
            kink_index, kink_fraction, phi_re,
        )
    finally:
        if threads_before != CONFIG.threads:
            numba.set_num_threads(threads_before)
    omega = outputs[0][:, -1] - outputs[1][:, -1]
    n_re_abs = float(model.derived_param["N_inf"] - model.derived_param["N_re"])
    n_re = float(model.derived_param["N_re"])
    s_start = (Nv[j0s] - n_re_abs) / n_re
    return {
        "DN_eff_fixed": float(model.cosmo_param["DN_eff"]),
        "frequency_log10_hz": probe_freqs.tolist(),
        "j0_indices": j0s.astype(int).tolist(),
        "j0_normalized_s": s_start.tolist(),
        "starts_after_reheating": (s_start >= 0.0).tolist(),
        "log10_OmegaGW": np.log10(omega).tolist(),
    }


def _run_case(dn_re: float, probe_freqs: np.ndarray) -> dict:
    model = LCDM_SG(**PARAMS, DN_re=dn_re)
    dn_eff_initial = float(model.cosmo_param["DN_eff"])
    Nv, sigma, f_hor, _ = build_transition_grid(model, CONFIG.h)
    derived = model.derived_param
    n_re_abs = float(derived["N_inf"] - derived["N_re"])
    n_re = float(derived["N_re"])
    sample_N = n_re_abs + S_SAMPLES * n_re

    sigma_sample = sigma_vec(sample_N, model, model.cosmo_param["DN_eff"])
    h2_sample = H2_vec(sample_N, model, model.cosmo_param["DN_eff"])
    phi, _, s2, _, _ = exact_phi_s2_grid(
        model, Nv, model.cosmo_param["DN_eff"]
    )
    s_grid = (Nv - n_re_abs) / n_re
    phi_sample = np.interp(S_SAMPLES, s_grid, phi)
    log_s2_sample = np.interp(S_SAMPLES, s_grid, np.log(s2))
    f_hor_sample = np.interp(S_SAMPLES, s_grid, f_hor)

    start = time.perf_counter()
    result = FS.SGWB_iter_fast(
        model,
        tol=1e-6,
        freq_res=0.25,
        config=CONFIG,
        eval_freqs=probe_freqs,
        physical_guard=True,
    )
    solve_seconds = time.perf_counter() - start

    record = {
        "input": {**PARAMS, "DN_re": dn_re, "DN_eff_initial": dn_eff_initial},
        "derived": {
            "N_re": n_re,
            "N_inf": float(derived["N_inf"]),
            "N_re_abs": n_re_abs,
        },
        "background": {
            "sigma_s": sigma_sample.tolist(),
            "H2_s": h2_sample.tolist(),
            "Phi_s": phi_sample.tolist(),
            "log_S2_s": log_s2_sample.tolist(),
            "f_hor_s": f_hor_sample.tolist(),
            "f_re": float(model.f_re),
            "fmax": float(f_hor[0]),
        },
        "fast": {
            "success": result is model,
            "runtime_s": solve_seconds,
            "failure_reason": getattr(model, "fast_failure_reason", None),
            "n_frequency": int(model.f.size),
            "frequency_range": [float(model.f.min()), float(model.f.max())],
        },
    }
    if result is model:
        indices = [int(np.argmin(np.abs(model.f - f))) for f in probe_freqs]
        record["fast"].update(
            DN_gw_final=float(model.DN_gw[-1]),
            DN_eff_final=float(model.cosmo_param["DN_eff"]),
            matched_probe_frequencies=model.f[indices].astype(float).tolist(),
            probe_log10_OmegaGW=model.log10OmegaGW[indices].astype(float).tolist(),
        )
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "docs" / "normalized_background_dnre_probe_20261105.json",
    )
    args = parser.parse_args()

    probe_models = [LCDM_SG(**PARAMS, DN_re=dn_re) for dn_re in DN_RE_VALUES]
    for model in probe_models:
        build_transition_grid(model, CONFIG.h)
    f_re = float(probe_models[0].f_re)
    probe_freqs = f_re + F_OFFSETS
    fixed_map_freqs = f_re + FIXED_MAP_OFFSETS

    numba_threads_before_solve = FS.get_num_threads()
    cases = [_run_case(dn_re, probe_freqs) for dn_re in DN_RE_VALUES]
    fixed_maps = [_run_fixed_map(dn_re, fixed_map_freqs) for dn_re in DN_RE_VALUES]
    numba_threads_after_solve = FS.get_num_threads()
    lhs, rhs = (case["background"] for case in cases)

    def max_abs(a: list[float], b: list[float]) -> float:
        return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))

    phi_offset = [
        np.asarray(case["background"]["Phi_s"]) - 0.5 * DN_RE_VALUES[i]
        for i, case in enumerate(cases)
    ]
    log_s2_offset = [
        np.asarray(case["background"]["log_S2_s"]) + DN_RE_VALUES[i]
        for i, case in enumerate(cases)
    ]
    probe_delta = None
    if all(case["fast"]["success"] for case in cases):
        probe_delta = max_abs(
            cases[0]["fast"]["probe_log10_OmegaGW"],
            cases[1]["fast"]["probe_log10_OmegaGW"],
        )
    fixed_map_delta = np.abs(
        np.asarray(fixed_maps[0]["log10_OmegaGW"])
        - np.asarray(fixed_maps[1]["log10_OmegaGW"])
    )
    fixed_map_after_re = np.logical_and(
        fixed_maps[0]["starts_after_reheating"], fixed_maps[1]["starts_after_reheating"]
    )
    fixed_map_interior = np.logical_and(
        np.asarray(fixed_maps[0]["j0_indices"]) > 0,
        np.asarray(fixed_maps[1]["j0_indices"]) > 0,
    )

    try:
        source_head = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        source_head = None

    payload = {
        "experiment": "reheating-aligned normalized-coordinate and full-Fast probe",
        "source_head_before_probe": source_head,
        "source_sha256": {
            "LCDM_stiff_Neff.py": _sha256(ROOT / "stiffgwpy_fast" / "LCDM_stiff_Neff.py"),
            "exact_background.py": _sha256(ROOT / "stiffgwpy_fast" / "exact_background.py"),
            "fast_sgwb.py": _sha256(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "probe_script.py": _sha256(Path(__file__).resolve()),
        },
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "numba": numba.__version__,
            "numba_thread_limit": numba.config.NUMBA_NUM_THREADS,
            "numba_threads_before_solve": numba_threads_before_solve,
            "configured_threads": CONFIG.threads,
            "numba_threads_after_solve": numba_threads_after_solve,
            "thread_environment": {
                name: os.environ.get(name)
                for name in ("NUMBA_NUM_THREADS", "NUMBA_THREADING_LAYER", "FAST_THREADS")
            },
        },
        "configuration": {
            "h": CONFIG.h,
            "col_step": CONFIG.col_step,
            "z_tail": CONFIG.z_tail,
            "phase_max": CONFIG.phase_max,
            "freq_grid": CONFIG.freq_grid,
            "freq_res": 0.25,
            "tol": 1e-6,
            "physical_guard": True,
            "frequency_quadrature": "pchip",
        },
        "normalized_coordinate": "s=(N-N_re_abs)/N_re on N>=N_re_abs; N_inf-N=(1-s)N_re",
        "s_samples": S_SAMPLES.tolist(),
        "probe_frequencies_log10_hz": probe_freqs.tolist(),
        "cases": cases,
        "fixed_DN_eff_complete_kernel": fixed_maps,
        "comparison": {
            "N_re_abs_shift": float(cases[1]["derived"]["N_re_abs"] - cases[0]["derived"]["N_re_abs"]),
            "max_abs_sigma_difference": max_abs(lhs["sigma_s"], rhs["sigma_s"]),
            "max_relative_H2_difference": float(
                np.max(np.abs(np.asarray(lhs["H2_s"]) - np.asarray(rhs["H2_s"]))
                       / np.maximum(np.abs(lhs["H2_s"]), 1e-300))
            ),
            "max_abs_f_hor_difference": max_abs(lhs["f_hor_s"], rhs["f_hor_s"]),
            "max_abs_Phi_after_0p5DNre_shift": float(np.max(np.abs(phi_offset[0] - phi_offset[1]))),
            "max_abs_logS2_after_DNre_shift": float(np.max(np.abs(log_s2_offset[0] - log_s2_offset[1]))),
            "max_abs_log10_spectrum_difference_at_shared_probe_frequencies": probe_delta,
            "fixed_map_max_abs_log10_spectrum_difference_all": float(np.max(fixed_map_delta)),
            "fixed_map_max_abs_log10_spectrum_difference_when_both_start_after_reheating": (
                float(np.max(fixed_map_delta[fixed_map_after_re]))
                if np.any(fixed_map_after_re) else None
            ),
            "fixed_map_after_reheating_mode_count": int(np.count_nonzero(fixed_map_after_re)),
            "fixed_map_interior_start_mode_count": int(np.count_nonzero(fixed_map_interior)),
            "fixed_map_max_abs_log10_spectrum_difference_for_interior_starts": (
                float(np.max(fixed_map_delta[fixed_map_interior]))
                if np.any(fixed_map_interior) else None
            ),
            "interpretation_limit": "Cross-parameter output differences are physical sensitivity, not numerical error; this is not an independent accuracy or parameter-domain certification.",
            "timing_limit": "Sequential single-call wall times are recorded per case; the first call can include JIT cost and must not be compared as a performance result.",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(_json_safe(payload), indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"output": str(args.output), "comparison": payload["comparison"],
                      "case_success": [case["fast"]["success"] for case in cases]}, indent=2))
    return 0 if all(case["fast"]["success"] for case in cases) else 1


if __name__ == "__main__":
    raise SystemExit(main())
