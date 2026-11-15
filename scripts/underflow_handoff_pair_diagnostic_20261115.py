"""Record Fast/Reference hand-off diagnostics at the underflow spectrum worst node."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
REFERENCE_TAIL_INPUT = ROOT / "docs" / "research" / "results" / "underflow_reference_tail_z8_z10_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_handoff_pair_diagnostic_20261115.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fast_fixed(data: dict, freqs: np.ndarray, frozen_dn: float, z_tail: float) -> dict:
    params = dict(data["parameters"])
    model = LCDM_SG(**params)
    model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal, original_construct = FA.goal_oriented_freqs, FA.construct_fast_frequency_grid

    def fixed_grid(*args, **kwargs):
        return freqs.copy()

    def locked_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = float(frozen_dn)
        return original_gen(*args, **kwargs)

    FA.goal_oriented_freqs = fixed_grid
    FA.construct_fast_frequency_grid = fixed_grid
    model.gen_expansion = locked_gen
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=0.99,
            z_tail=z_tail,
            threads=1,
            accuracy_mode="fast",
            freq_grid="goal",
            physical_guard=False,
        )
    finally:
        FA.goal_oriented_freqs, FA.construct_fast_frequency_grid = original_goal, original_construct
    elapsed = time.perf_counter() - started
    if result is None or not np.array_equal(model.f, freqs):
        raise RuntimeError(f"Fast handoff run failed: {model.fast_failure_reason}")
    return {
        "z_tail": z_tail,
        "runtime_s": elapsed,
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "fixed_map_residual": float(model.DN_gw[-1] - frozen_dn),
        "log10_OmegaGW": np.asarray(model.log10OmegaGW, dtype=np.float64).tolist(),
        "Ogw_today": np.asarray(model.Ogw_today, dtype=np.float64).tolist(),
        "Oj_today": np.asarray(model.Oj_today, dtype=np.float64).tolist(),
        "Opgw_today": np.asarray(model.Opgw_today, dtype=np.float64).tolist(),
        "linear_OmegaGW": np.asarray(model.Ogw_today - model.Oj_today, dtype=np.float64).tolist(),
        "handoff_eps": np.asarray(model.handoff_eps, dtype=np.float64).tolist(),
    }


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    reference_tail = json.loads(REFERENCE_TAIL_INPUT.read_text(encoding="utf-8"))
    freqs = np.asarray(data["frequency_grids"]["235"]["log10_hz"], dtype=np.float64)
    common_dn = float(data["background"]["common_DN_eff"])
    reference_z10_log = np.asarray(reference_tail["z_tail_10"]["log10_OmegaGW"], dtype=np.float64)
    fast5 = _fast_fixed(data, freqs, common_dn, 5.0)
    fast8 = _fast_fixed(data, freqs, common_dn, 8.0)
    fast12 = _fast_fixed(data, freqs, common_dn, 12.0)
    error5 = np.abs(np.asarray(fast5["log10_OmegaGW"]) - reference_z10_log)
    error8 = np.abs(np.asarray(fast8["log10_OmegaGW"]) - reference_z10_log)
    error12 = np.abs(np.asarray(fast12["log10_OmegaGW"]) - reference_z10_log)
    max_index = int(np.argmax(error12))
    frequency = float(freqs[max_index])
    old_worst_index = 209

    unit_params = dict(data["parameters"])
    unit_params.update({"A_s": 1.0, "r": 1.0, "n_t": 0.0})
    unit_model = LCDM_SG(**unit_params)
    unit_model.gen_expansion()
    def reference_mode_log(freq: float, z_tail: float) -> tuple[dict, float]:
        mode = REF.solve_reference_mode(
            unit_model,
            freq,
            common_dn,
            z_tail=z_tail,
            rtol=1e-9,
            backend="cartesian",
        )
        log_amp = _log_tensor_amplitude(data["parameters"], {"A_t": 0.0})
        log_power = log_amp + float(data["parameters"]["n_t"]) * (
            math.log(10.0) * freq - math.log(REF.gp.f_piv)
        )
        base = 3.0 * mode["Opgw_today"] if mode["used_tail"] else mode["Ogw_today"] - mode["Oj_today"]
        return mode, (math.log(base) + log_power) / math.log(10.0)

    ref_max_z8, ref_max_log_z8 = reference_mode_log(frequency, 8.0)
    ref_max_z10, ref_max_log_z10 = reference_mode_log(frequency, 10.0)
    old_frequency = float(freqs[old_worst_index])
    ref_old_z10, ref_old_log_z10 = reference_mode_log(old_frequency, 10.0)

    def reference_unit_moments(mode: dict, freq: float) -> dict:
        power_log = _log_tensor_amplitude(data["parameters"], {"A_t": 0.0}) + float(
            data["parameters"]["n_t"]
        ) * (math.log(10.0) * freq - math.log(REF.gp.f_piv))
        scale = math.exp(power_log)
        return {
            "unit_Ogw": float(mode["Ogw_today"]),
            "unit_Oj": float(mode["Oj_today"]),
            "unit_Opgw": float(mode["Opgw_today"]),
            "scaled_Ogw": float(mode["Ogw_today"] * scale),
            "scaled_Oj": float(mode["Oj_today"] * scale),
            "scaled_Opgw": float(mode["Opgw_today"] * scale),
            "log_P_t": float(power_log),
        }

    output = {
        "schema_version": 1,
        "status": "complete_single_mode_handoff_diagnostic",
        "parameters": data["parameters"],
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": data["frequency_grids"]["235"]["sha256_f64le"],
        "frozen_DN_eff": common_dn,
        "worst_node": {
            "index": max_index,
            "log10_f_Hz": frequency,
            "Reference_z8_used_tail": bool(ref_max_z8["used_tail"]),
            "Reference_z8_event_N": ref_max_z8["event_N"],
            "Reference_z8_handoff_eps": ref_max_z8["eps_handoff"],
            "Reference_z8_n_steps": ref_max_z8["n_steps"],
            "Reference_z10_used_tail": bool(ref_max_z10["used_tail"]),
            "Reference_z10_event_N": ref_max_z10["event_N"],
            "Reference_z10_handoff_eps": ref_max_z10["eps_handoff"],
            "Reference_z10_n_steps": ref_max_z10["n_steps"],
            "Reference_z10_log10_OmegaGW": ref_max_log_z10,
            "Fast_z5_log10_OmegaGW": fast5["log10_OmegaGW"][max_index],
            "Fast_z8_log10_OmegaGW": fast8["log10_OmegaGW"][max_index],
            "Fast_z12_log10_OmegaGW": fast12["log10_OmegaGW"][max_index],
            "Fast_z5_handoff_eps": fast5["handoff_eps"][max_index],
            "Fast_z8_handoff_eps": fast8["handoff_eps"][max_index],
            "Fast_z12_handoff_eps": fast12["handoff_eps"][max_index],
            "Fast_z5_terminal_moments": {
                name: fast5[name][max_index] for name in ("Ogw_today", "Oj_today", "Opgw_today")
            },
            "Fast_z8_terminal_moments": {
                name: fast8[name][max_index] for name in ("Ogw_today", "Oj_today", "Opgw_today")
            },
            "Fast_z12_terminal_moments": {
                name: fast12[name][max_index] for name in ("Ogw_today", "Oj_today", "Opgw_today")
            },
            "Reference_z8_unit_terminal_moments": reference_unit_moments(ref_max_z8, frequency),
            "Reference_z10_unit_terminal_moments": reference_unit_moments(ref_max_z10, frequency),
            "Fast_z5_abs_error_vs_Reference_z10_dex": float(error5[max_index]),
            "Fast_z8_abs_error_vs_Reference_z10_dex": float(error8[max_index]),
            "Fast_z12_abs_error_vs_Reference_z10_dex": float(error12[max_index]),
        },
        "fast_spectrum_summary": {
            "Reference_tail": 10.0,
            "z5_max_p95_error_dex": [float(np.max(error5)), float(np.quantile(error5, 0.95))],
            "z8_max_p95_error_dex": [float(np.max(error8)), float(np.quantile(error8, 0.95))],
            "z12_max_p95_error_dex": [float(np.max(error12)), float(np.quantile(error12, 0.95))],
            "z8_minus_z5_max_abs_dex": float(np.max(np.abs(np.asarray(fast8["log10_OmegaGW"]) - np.asarray(fast5["log10_OmegaGW"])) ))
        },
        "previous_worst_node_z8_reference": {
            "index": old_worst_index,
            "log10_f_Hz": old_frequency,
            "Reference_z10_log10_OmegaGW": ref_old_log_z10,
            "Reference_z10_handoff_eps": ref_old_z10["eps_handoff"],
            "Fast_z8_log10_OmegaGW": fast8["log10_OmegaGW"][old_worst_index],
            "Fast_z12_log10_OmegaGW": fast12["log10_OmegaGW"][old_worst_index],
            "Fast_z8_handoff_eps": fast8["handoff_eps"][old_worst_index],
            "Fast_z12_handoff_eps": fast12["handoff_eps"][old_worst_index],
            "Fast_z8_terminal_moments": {
                name: fast8[name][old_worst_index] for name in ("Ogw_today", "Oj_today", "Opgw_today")
            },
            "Fast_z12_terminal_moments": {
                name: fast12[name][old_worst_index] for name in ("Ogw_today", "Oj_today", "Opgw_today")
            },
            "Reference_z10_unit_terminal_moments": reference_unit_moments(ref_old_z10, old_frequency),
        },
        "fast_z5": fast5,
        "fast_z8": fast8,
        "fast_z12": fast12,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive.py": _sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "input_json": _sha(INPUT),
            "reference_tail_input_json": _sha(REFERENCE_TAIL_INPUT),
            "diagnostic_script": _sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(output["worst_node"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
