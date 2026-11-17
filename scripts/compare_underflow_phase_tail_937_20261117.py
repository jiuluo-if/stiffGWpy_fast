"""English description: Isolate the phase-resolved terminal correction on one frozen 937-node underflow case. 中文说明：在固定背景、频率网格和线程下单独比较相位尾部修正。"""

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

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude, _requires_log_tensor_power  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_current_20261117.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_phase_tail_937_ab_20261117.json"
FREQUENCY_COUNT = 937
Z_TAIL = 5.0
THREADS = 1


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes(order="C")).hexdigest()


def error_stats(a: np.ndarray, b: np.ndarray) -> dict[str, float | int]:
    delta = np.abs(np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64))
    if delta.size == 0 or np.any(~np.isfinite(delta)):
        raise FloatingPointError("Spectrum comparison contains non-finite values")
    return {
        "count": int(delta.size),
        "max_abs_dex": float(np.max(delta)),
        "p95_abs_dex": float(np.quantile(delta, 0.95)),
    }


def run_case(
    params: dict[str, object], freqs: np.ndarray, frozen_dn: float, phase_enabled: bool
) -> dict[str, object]:
    model = LCDM_SG(**params)
    model.cosmo_param["DN_eff"] = frozen_dn
    model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    original_phase_weights = FS._phase_resolved_tail_weights
    phase_calls: list[dict[str, object]] = []

    def fixed_grid(*args, **kwargs):
        return freqs.copy()

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = frozen_dn
        return original_gen(*args, **kwargs)

    if phase_enabled:
        def capture_phase_weights(z_today, z_tail, enabled):
            weights = original_phase_weights(z_today, z_tail, enabled)
            phase_calls.append({
                "mode_count": int(weights.size),
                "active_mode_count": int(np.count_nonzero(weights > 0.0)),
                "max_weight": float(np.max(weights)) if weights.size else 0.0,
            })
            return weights
        phase_callable = capture_phase_weights
    else:
        def phase_callable(z_today, z_tail, enabled):
            weights = np.zeros_like(np.asarray(z_today, dtype=np.float64))
            phase_calls.append({"mode_count": int(weights.size), "active_mode_count": 0, "max_weight": 0.0})
            return weights

    FA.goal_oriented_freqs = fixed_grid
    FA.construct_fast_frequency_grid = fixed_grid
    model.gen_expansion = fixed_background_gen
    FS._phase_resolved_tail_weights = phase_callable
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=0.99,
            z_tail=Z_TAIL,
            threads=THREADS,
            accuracy_mode="fast",
            freq_grid="goal",
            physical_guard=False,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
        FS._phase_resolved_tail_weights = original_phase_weights
    runtime = time.perf_counter() - started
    if result is None:
        raise RuntimeError(f"Fast phase-tail run failed: {model.fast_failure_reason}")
    if not np.array_equal(np.asarray(model.f, dtype=np.float64), freqs):
        raise AssertionError("Fast did not use the exact requested frequency grid")
    logs = np.asarray(model.log10OmegaGW, dtype=np.float64)
    if logs.size != freqs.size or not np.isfinite(logs).all():
        raise FloatingPointError("Fast log spectrum is invalid")
    mapped = float(model.DN_gw[-1])
    if not math.isfinite(mapped) or mapped <= 0.0:
        raise FloatingPointError("Fast fixed-background map is invalid")
    derived = model.derived_param
    log_amplitude = _log_tensor_amplitude(model.cosmo_param, derived)
    expected_log_mode = (
        (float(derived["A_t"]) == 0.0 or not math.isfinite(float(derived["A_t"])))
        or _requires_log_tensor_power(log_amplitude, derived["nt"], freqs, FS.gp.f_piv)
    )
    return {
        "phase_enabled": phase_enabled,
        "runtime_s": runtime,
        "frozen_DN_eff": frozen_dn,
        "mapped_DN_gw": mapped,
        "map_residual_relative": mapped / frozen_dn - 1.0,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": digest(freqs),
        "log_amplitude_mode_expected": bool(expected_log_mode),
        "A_t_stored": float(derived["A_t"]),
        "phase_weight_calls": phase_calls,
        "log10_OmegaGW": logs.tolist(),
    }


def main() -> int:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    if source.get("status") != "complete_fixed_background_direct_nested_grid_screen":
        raise RuntimeError("Current-source direct nested-grid screen is incomplete")
    for rel, expected in source["source_sha256"].items():
        if sha256(ROOT / rel) != expected:
            raise RuntimeError(f"Input source hash differs for {rel}")
    count_key = str(FREQUENCY_COUNT)
    grid_data = source["frequency_grids"][count_key]
    freqs = np.asarray(grid_data["log10_hz"], dtype=np.float64)
    if freqs.size != FREQUENCY_COUNT or digest(freqs) != grid_data["sha256_f64le"]:
        raise RuntimeError("937-point grid hash mismatch")
    frozen_dn = float(source["fast_self_consistent"][count_key]["root_DN_eff"])
    params = {key: float(value) for key, value in source["parameters"].items()}

    enabled = run_case(params, freqs, frozen_dn, True)
    disabled = run_case(params, freqs, frozen_dn, False)
    ref_logs = np.asarray(source["reference_map_at_fast_root"][count_key]["log10_OmegaGW"], dtype=np.float64)
    enabled_logs = np.asarray(enabled["log10_OmegaGW"], dtype=np.float64)
    disabled_logs = np.asarray(disabled["log10_OmegaGW"], dtype=np.float64)
    output = {
        "schema_version": 1,
        "status": "complete_frozen_root_phase_tail_937_ab",
        "scope": "Same parameter values, frozen Fast root, exact 937 direct frequencies, z_tail, and one thread; the only changed factor is whether phase-resolved terminal moments are applied. This is a numerical comparison, not a performance claim.",
        "input_artifact": INPUT.relative_to(ROOT).as_posix(),
        "parameters": params,
        "z_tail": Z_TAIL,
        "threads": THREADS,
        "frequency_count": FREQUENCY_COUNT,
        "frequency_sha256_f64le": digest(freqs),
        "frozen_DN_eff": frozen_dn,
        "cases": {
            "phase_resolved": enabled,
            "phase_averaged": disabled,
        },
        "phase_resolved_vs_phase_averaged_spectrum_dex": error_stats(enabled_logs, disabled_logs),
        "phase_resolved_vs_reference_z8_spectrum_dex": error_stats(enabled_logs, ref_logs),
        "phase_averaged_vs_reference_z8_spectrum_dex": error_stats(disabled_logs, ref_logs),
        "mapped_DN_gw_relative_difference_phase_resolved_vs_phase_averaged": (
            float(enabled["mapped_DN_gw"]) / float(disabled["mapped_DN_gw"]) - 1.0
        ),
        "source_sha256": {
            **source["source_sha256"],
            "input_artifact": sha256(INPUT),
            "comparison_runner": sha256(Path(__file__)),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
    }
    OUTPUT.write_text(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "active_phase_modes": enabled["phase_weight_calls"][-1]["active_mode_count"] if enabled["phase_weight_calls"] else 0,
        "phase_spectrum_delta": output["phase_resolved_vs_phase_averaged_spectrum_dex"],
        "phase_root_delta": output["mapped_DN_gw_relative_difference_phase_resolved_vs_phase_averaged"],
        "phase_vs_reference": output["phase_resolved_vs_reference_z8_spectrum_dex"],
        "average_vs_reference": output["phase_averaged_vs_reference_z8_spectrum_dex"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
