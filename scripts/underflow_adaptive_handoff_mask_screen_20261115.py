"""Screen a z_today-gated per-mode WKB handoff using paired Fast outputs."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.reference import _Background, f_hor_abs  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_grid_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_adaptive_handoff_mask_screen_20261115.json"
BASE_GRID_INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
BASE_TAIL_INPUT = ROOT / "docs" / "research" / "results" / "underflow_tail_z5_z8_compare_20261115.json"
BASE_REF_INPUT = ROOT / "docs" / "research" / "results" / "underflow_reference_tail_z8_z10_20261115.json"
THRESHOLDS = (5.0, 6.0, 7.0, 8.0, 9.0, 10.0)


def run_fast(params: dict[str, float], freqs: np.ndarray, root: float, z_tail: float) -> tuple[dict, float]:
    model = LCDM_SG(**params)
    model.cosmo_param["DN_eff"] = root
    model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = root
        return original_gen(*args, **kwargs)

    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    model.gen_expansion = fixed_background_gen
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        try:
            result = model.SGWB_iter(
                engine="fast", fallback=False, tol=0.99, z_tail=z_tail,
                threads=1, accuracy_mode="fast", freq_grid="goal",
                physical_guard=False,
            )
        finally:
            FA.goal_oriented_freqs = original_goal
            FA.construct_fast_frequency_grid = original_construct
    elapsed = time.perf_counter() - started
    if result is None:
        raise RuntimeError(f"Fast z_tail={z_tail} failed: {model.fast_failure_reason}")
    return {
        "log10_OmegaGW": np.asarray(model.log10OmegaGW, dtype=np.float64),
        "Ogw_today": np.asarray(model.Ogw_today, dtype=np.float64),
        "Oj_today": np.asarray(model.Oj_today, dtype=np.float64),
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "runtime_s": elapsed,
        "runtime_warnings": [str(item.message) for item in caught],
        "frequency_count": int(model.f.size),
    }, float(model.derived_param["h"])


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    row = next(item for item in data["cases"] if item["frequency_count"] == 937)
    params = data["parameters"]
    freqs = np.asarray(row["frequency_log10_hz"], dtype=np.float64)
    root = float(row["frozen_DN_eff"])
    fast5, h = run_fast(params, freqs, root, 5.0)
    fast8, _ = run_fast(params, freqs, root, 8.0)

    bg_model = LCDM_SG(**params)
    bg_model.cosmo_param["DN_eff"] = root
    bg_model.gen_expansion()
    bg = _Background(bg_model, root)
    f_hor_today = f_hor_abs(bg, bg.N_inf)
    z_today = (freqs - f_hor_today) * gp.ln10
    ref_log = np.asarray(row["reference"]["log10_OmegaGW"], dtype=np.float64)
    map_scale = gp.Neff0 * gp.ln10 / (gp.Omega_nh2 / (h * h))
    candidates = []
    for threshold in THRESHOLDS:
        use_z8 = z_today < threshold
        spectrum = np.where(use_z8, fast8["log10_OmegaGW"], fast5["log10_OmegaGW"])
        errors = spectrum - ref_log
        integrand = np.where(use_z8, fast8["Ogw_today"] - fast8["Oj_today"], fast5["Ogw_today"] - fast5["Oj_today"])
        with warnings.catch_warnings(record=True) as caught_map:
            warnings.simplefilter("always", RuntimeWarning)
            map_value = map_scale * float(np.sum(FS._pchip_integrals_vectorized(freqs, integrand)))
        candidates.append({
            "z_today_threshold": threshold,
            "z8_mode_count": int(np.count_nonzero(use_z8)),
            "z8_mode_fraction": float(np.mean(use_z8)),
            "max_abs_spectrum_error_vs_reference_z8_dex": float(np.max(np.abs(errors))),
            "p95_abs_spectrum_error_vs_reference_z8_dex": float(np.quantile(np.abs(errors), 0.95)),
            "max_error_frequency_log10_hz": float(freqs[np.argmax(np.abs(errors))]),
            "fixed_background_mapped_DN_gw": float(map_value),
            "relative_map_error_vs_reference_z8": float(map_value / float(row["reference"]["mapped_DN_gw"]) - 1.0),
            "relative_map_residual_vs_frozen_root": float(map_value / root - 1.0),
            "mapping_warnings": [str(item.message) for item in caught_map],
        })
    base_grid_data = json.loads(BASE_GRID_INPUT.read_text(encoding="utf-8"))
    base_tail_data = json.loads(BASE_TAIL_INPUT.read_text(encoding="utf-8"))
    base_ref_data = json.loads(BASE_REF_INPUT.read_text(encoding="utf-8"))
    base_freqs = np.asarray(base_grid_data["frequency_grids"]["937"]["log10_hz"], dtype=np.float64)
    base_z5 = np.asarray(base_grid_data["fast_fixed_background"]["937"]["log10_OmegaGW"], dtype=np.float64)
    base_z8 = np.asarray(base_tail_data["grids"]["937"]["z_tail_8"]["log10_OmegaGW"], dtype=np.float64)
    base_ref10 = np.asarray(base_ref_data["grids"]["937"]["z_tail_10"]["log10_OmegaGW"], dtype=np.float64)
    base_root = float(base_grid_data["background"]["common_DN_eff"])
    base_params = base_grid_data["parameters"]
    base_model = LCDM_SG(**base_params)
    base_model.gen_expansion()
    base_bg = _Background(base_model, base_root)
    base_f_hor_today = f_hor_abs(base_bg, base_bg.N_inf)
    base_z_today = (base_freqs - base_f_hor_today) * gp.ln10
    base_h = float(base_params["H0"]) / 100.0
    base_map_scale = gp.Neff0 * gp.ln10 / (gp.Omega_nh2 / (base_h * base_h))
    baseline_candidates = []
    for threshold in THRESHOLDS:
        use_z8 = base_z_today < threshold
        spectrum = np.where(use_z8, base_z8, base_z5)
        errors = spectrum - base_ref10
        ln_spectrum = spectrum * gp.ln10
        shift = float(np.max(ln_spectrum))
        scaled = np.exp(ln_spectrum - shift)
        with warnings.catch_warnings(record=True) as caught_map:
            warnings.simplefilter("always", RuntimeWarning)
            scaled_integral = float(np.sum(FS._pchip_integrals_vectorized(base_freqs, scaled)))
        mapped = base_map_scale * math.exp(shift) * scaled_integral
        baseline_candidates.append({
            "z_today_threshold": threshold,
            "z8_mode_count": int(np.count_nonzero(use_z8)),
            "z8_mode_fraction": float(np.mean(use_z8)),
            "max_abs_spectrum_error_vs_reference_z10_dex": float(np.max(np.abs(errors))),
            "p95_abs_spectrum_error_vs_reference_z10_dex": float(np.quantile(np.abs(errors), 0.95)),
            "max_error_frequency_log10_hz": float(base_freqs[np.argmax(np.abs(errors))]),
            "fixed_background_mapped_DN_gw_log_reconstructed": float(mapped),
            "relative_map_error_vs_reference_z10": float(
                mapped / float(base_ref_data["grids"]["937"]["z_tail_10"]["mapped_DN_gw"]) - 1.0
            ),
            "relative_map_residual_vs_frozen_root": float(mapped / base_root - 1.0),
            "mapping_warnings": [str(item.message) for item in caught_map],
        })
    output = {
        "schema_version": 1,
        "status": "complete_z_today_gated_handoff_hypothesis_screen_not_production",
        "scope": "Compose per-frequency z8 rows only when today's z is below a threshold, otherwise z5 rows; both inputs are direct same-root Fast runs and scores use direct Reference z8 nodes.",
        "parameters": params,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": row["frequency_sha256_f64le"],
        "frozen_DN_eff": root,
        "f_hor_today_log10_hz": f_hor_today,
        "z_today_range": [float(np.min(z_today)), float(np.max(z_today))],
        "fast_direct_runs": {
            "z5_runtime_s": fast5["runtime_s"],
            "z5_runtime_warnings": fast5["runtime_warnings"],
            "z5_mapped_DN_gw": fast5["mapped_DN_gw"],
            "z8_runtime_s": fast8["runtime_s"],
            "z8_runtime_warnings": fast8["runtime_warnings"],
            "z8_mapped_DN_gw": fast8["mapped_DN_gw"],
        },
        "reference_z8_map_DN_gw": row["reference"]["mapped_DN_gw"],
        "threshold_candidates": candidates,
        "amplitude_holdout": {
            "label": "r=1e-316 baseline amplitude at its own frozen 235-root",
            "root_DN_eff": base_root,
            "frequency_count": int(base_freqs.size),
            "frequency_sha256_f64le": base_grid_data["frequency_grids"]["937"]["sha256_f64le"],
            "f_hor_today_log10_hz": base_f_hor_today,
            "z_today_range": [float(np.min(base_z_today)), float(np.max(base_z_today))],
            "reference_z10_map_DN_gw": base_ref_data["grids"]["937"]["z_tail_10"]["mapped_DN_gw"],
            "threshold_candidates": baseline_candidates,
        },
        "runtime_note": "These z5 then z8 single calls are not alternated or cold/warm controlled and are not comparable timings. Composition is an accuracy screen, not the runtime or memory cost of a per-mode production kernel.",
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            name: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            for name, path in {
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "input_json": "docs/research/results/underflow_min_subnormal_reference_grid_20261115.json",
                "base_grid_input": "docs/research/results/underflow_grid_235_469_937_20261115.json",
                "base_tail_input": "docs/research/results/underflow_tail_z5_z8_compare_20261115.json",
                "base_reference_input": "docs/research/results/underflow_reference_tail_z8_z10_20261115.json",
                "diagnostic_script": "scripts/underflow_adaptive_handoff_mask_screen_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "direct_call_seconds_not_comparable": [fast5["runtime_s"], fast8["runtime_s"]],
        "direct_call_warning_counts": [len(fast5["runtime_warnings"]), len(fast8["runtime_warnings"])],
        "candidates": [{
            "threshold": row["z_today_threshold"],
            "z8_fraction": row["z8_mode_fraction"],
            "max_spectrum_error": row["max_abs_spectrum_error_vs_reference_z8_dex"],
            "p95_spectrum_error": row["p95_abs_spectrum_error_vs_reference_z8_dex"],
            "map_error": row["relative_map_error_vs_reference_z8"],
            "root_residual": row["relative_map_residual_vs_frozen_root"],
        } for row in candidates],
        "amplitude_holdout_candidates": [{
            "threshold": row["z_today_threshold"],
            "z8_fraction": row["z8_mode_fraction"],
            "max_spectrum_error": row["max_abs_spectrum_error_vs_reference_z10_dex"],
            "p95_spectrum_error": row["p95_abs_spectrum_error_vs_reference_z10_dex"],
            "map_error": row["relative_map_error_vs_reference_z10"],
        } for row in baseline_candidates],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
