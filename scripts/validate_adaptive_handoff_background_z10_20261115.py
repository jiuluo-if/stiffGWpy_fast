"""Validate the adaptive per-mode handoff on the low-T/high-tilt holdout vs Reference z10."""

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
from stiffgwpy_fast.reference import integrate_spectrum, spectrum_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_background_holdout_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_background_z10_20261115.json"
REFERENCE_Z8_CHECKPOINT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_background_reference_z8_20261115.json"
TAILS = (5.0, 8.0, 10.0)


def run_fast(params: dict[str, float], freqs: np.ndarray, root: float, *, z_tail: float, adaptive: bool) -> dict[str, object]:
    model = LCDM_SG(**params)
    model.cosmo_param["DN_eff"] = root
    if adaptive:
        model._fast_adaptive_handoff_below = 6.0
        model._fast_adaptive_handoff_to = 8.0
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
    try:
        result = model.SGWB_iter(
            engine="fast", fallback=False, tol=0.99,
            z_tail=z_tail, threads=1, accuracy_mode="fast",
            freq_grid="goal", physical_guard=True,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    elapsed = time.perf_counter() - started
    if result is None:
        raise RuntimeError(f"Fast tail={z_tail}, adaptive={adaptive} failed: {model.fast_failure_reason}")
    return {
        "log10_OmegaGW": np.asarray(model.log10OmegaGW, dtype=np.float64),
        "Ogw": np.asarray(model.Ogw_today, dtype=np.float64),
        "Oj": np.asarray(model.Oj_today, dtype=np.float64),
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "runtime_s_single_unpaired": elapsed,
        "selected_mode_count": getattr(model, "_fast_adaptive_handoff_selected", None),
        "selected_mode_fraction": getattr(model, "_fast_adaptive_handoff_fraction", None),
    }


def reference_run(model: LCDM_SG, freqs: np.ndarray, root: float, log_amp: float, nt: float, z_tail: float) -> dict[str, object]:
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        Ogw, Oj, Opgw, used_tail = spectrum_reference(
            model, freqs, root, z_tail=z_tail, rtol=1e-9,
            workers=1, backend="cartesian",
        )
        base = np.where(used_tail, 3.0 * Opgw, Ogw - Oj)
        if np.any(~np.isfinite(base)) or np.any(base <= 0.0):
            raise FloatingPointError("Reference unit spectrum is not positive finite")
        log_power = log_amp + nt * (math.log(10.0) * freqs - math.log(gp.f_piv))
        log_integrand = np.log(base) + log_power
        shift = float(np.max(log_integrand))
        normalized = np.exp(log_integrand - shift)
        g2_scaled, qerr, ierr = integrate_spectrum(freqs, normalized, np.zeros_like(normalized))
    h = float(model.derived_param["h"])
    map_scale = gp.Neff0 / (gp.Omega_nh2 / (h * h))
    mapped = map_scale * math.exp(shift + math.log(g2_scaled))
    return {
        "log10_OmegaGW": (log_integrand / gp.ln10).tolist(),
        "mapped_DN_gw": float(mapped),
        "map_relative_residual": float(mapped / root - 1.0),
        "runtime_s": time.perf_counter() - started,
        "z_tail": z_tail,
        "rtol": 1e-9,
        "used_tail_count": int(np.count_nonzero(used_tail)),
        "normalized_nodes_underflowed_to_zero": int(np.count_nonzero(normalized == 0.0)),
        "quadrature_error_relative": None if qerr is None else float(qerr / g2_scaled),
        "interpolation_error_relative": None if ierr is None else float(ierr / g2_scaled),
        "warnings": [str(item.message) for item in caught],
    }


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    params = data["parameters"]
    freqs = np.asarray(data["frequency_log10_hz"], dtype=np.float64)
    root = float(data["frozen_DN_eff"])
    log_amp = math.log(params["r"]) + math.log(params["A_s"])
    fast = {
        "z5": run_fast(params, freqs, root, z_tail=5.0, adaptive=False),
        "z8": run_fast(params, freqs, root, z_tail=8.0, adaptive=False),
        "adaptive": run_fast(params, freqs, root, z_tail=5.0, adaptive=True),
    }
    unit_params = {**params, "r": 1.0, "A_s": 1.0, "n_t": 0.0, "DN_eff": 0.0}
    ref_model = LCDM_SG(**unit_params)
    ref_model.gen_expansion()
    expected_grid_hash = hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest()
    if REFERENCE_Z8_CHECKPOINT.exists():
        checkpoint = json.loads(REFERENCE_Z8_CHECKPOINT.read_text(encoding="utf-8"))
        if (
            checkpoint["frequency_count"] != int(freqs.size)
            or checkpoint["frequency_sha256_f64le"] != expected_grid_hash
            or checkpoint["frozen_DN_eff"] != root
            or checkpoint["parameters"] != params
        ):
            raise AssertionError("saved Reference z8 checkpoint does not match this holdout")
    else:
        reference_z8 = reference_run(ref_model, freqs, root, log_amp, params["n_t"], 8.0)
        checkpoint = {
            "schema_version": 1,
            "status": "complete_background_holdout_reference_z8_checkpoint",
            "parameters": params,
            "frequency_count": int(freqs.size),
            "frequency_sha256_f64le": expected_grid_hash,
            "frozen_DN_eff": root,
            "reference": reference_z8,
            "source_sha256": {
                "reference.py": hashlib.sha256((ROOT / "stiffgwpy_fast/reference.py").read_bytes()).hexdigest(),
                "diagnostic_script": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            },
        }
        REFERENCE_Z8_CHECKPOINT.write_text(json.dumps(checkpoint, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    reference = {"8.0": checkpoint["reference"]}
    reference["10.0"] = reference_run(ref_model, freqs, root, log_amp, params["n_t"], 10.0)
    ref10 = reference["10.0"]
    errors = {}
    for name in fast:
        delta = fast[name]["log10_OmegaGW"] - ref10["log10_OmegaGW"]
        errors[name] = {
            "max_abs_dex": float(np.max(np.abs(delta))),
            "p95_abs_dex": float(np.quantile(np.abs(delta), 0.95)),
            "max_error_log10_hz": float(freqs[np.argmax(np.abs(delta))]),
        }
    adaptive_map = gp.Neff0 * gp.ln10 * float(np.sum(FS._pchip_integrals_vectorized(
        freqs, fast["adaptive"]["Ogw"] - fast["adaptive"]["Oj"]
    ))) / (gp.Omega_nh2 / float(ref_model.derived_param["h"]) ** 2)
    ref_shift = (
        np.asarray(reference["8.0"]["log10_OmegaGW"], dtype=np.float64)
        - np.asarray(ref10["log10_OmegaGW"], dtype=np.float64)
    )
    output = {
        "schema_version": 1,
        "status": "complete_adaptive_per_mode_reference_z10_background_holdout",
        "scope": "Actual private per-mode kernel at low-T/high-tilt cr0 holdout; direct Reference z8 and z10 on the same native 221-node grid and frozen 235-root.",
        "parameters": params,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
        "frozen_DN_eff": root,
        "adaptive_thresholds": {"z_today_below": 6.0, "defer_to_z": 8.0},
        "reference_z8_to_z10": {
            "max_abs_dex": float(np.max(np.abs(ref_shift))),
            "p95_abs_dex": float(np.quantile(np.abs(ref_shift), 0.95)),
            "map_relative_change": float(reference["8.0"]["mapped_DN_gw"] / ref10["mapped_DN_gw"] - 1.0),
        },
        "reference": reference,
        "fast": {
            name: {key: value for key, value in row.items() if key not in ("Ogw", "Oj", "log10_OmegaGW")}
            for name, row in fast.items()
        },
        "spectrum_error_vs_reference_z10": errors,
        "adaptive_map": {
            "mapped_DN_gw": float(adaptive_map),
            "relative_map_error_vs_reference_z10": float(adaptive_map / ref10["mapped_DN_gw"] - 1.0),
            "relative_fixed_root_residual": float(adaptive_map / root - 1.0),
        },
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
                "input_json": "docs/research/results/adaptive_handoff_background_holdout_20261115.json",
                "diagnostic_script": "scripts/validate_adaptive_handoff_background_z10_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "reference_z8_to_z10": output["reference_z8_to_z10"],
        "errors": errors,
        "adaptive_map": output["adaptive_map"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
