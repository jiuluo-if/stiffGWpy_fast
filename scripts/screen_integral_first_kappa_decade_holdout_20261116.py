"""Test integral-first scaling at a new kappa decade with direct Reference tails."""

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

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude, construct_fast_frequency_grid  # noqa: E402
from stiffgwpy_fast.reference import spectrum_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 5e-324,
    "n_t": 17.45,
    "cr": 0.0,
    "T_re": 10.0,
    "DN_re": 10.0,
    "kappa10": 1e-2,
}
OUTPUT = ROOT / "docs" / "research" / "results" / "integral_first_kappa_decade_holdout_20261116.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def grid_hash(freqs: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest()


def write_checkpoint(payload: dict[str, object]) -> None:
    payload["source_sha256"] = {
        name: sha(ROOT / path)
        for name, path in {
            "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
            "freq_adaptive.py": "stiffgwpy_fast/freq_adaptive.py",
            "exact_background.py": "stiffgwpy_fast/exact_background.py",
            "reference.py": "stiffgwpy_fast/reference.py",
            "stiff_SGWB.py": "stiffgwpy_fast/stiff_SGWB.py",
            "holdout_grid_input": "docs/research/results/adaptive_handoff_background_holdout_20261115.json",
            "reference_z10_input": "docs/research/results/adaptive_handoff_background_z10_20261115.json",
            "script": "scripts/screen_integral_first_kappa_decade_holdout_20261116.py",
        }.items()
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def refresh_metadata_hashes() -> int:
    """Refresh provenance after a metadata-only script edit without recomputing Reference."""
    if not OUTPUT.exists():
        raise FileNotFoundError(f"no completed checkpoint at {OUTPUT}")
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("status") != "complete_current_fast_kappa_holdout":
        raise ValueError("refuse metadata refresh on an incomplete Reference checkpoint")
    write_checkpoint(payload)
    print(json.dumps({"output": str(OUTPUT), "status": payload["status"], "updated_hashes": payload["source_sha256"]}, indent=2))
    return 0


def run_fast(freqs: np.ndarray, *, dn_eff: float, z_tail: float, tol: float) -> dict[str, object]:
    model = LCDM_SG(**{**PARAMS, "DN_eff": float(dn_eff)})
    model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = float(dn_eff)
        return original_gen(*args, **kwargs)

    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    model.gen_expansion = fixed_background_gen
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast", fallback=False, tol=tol, z_tail=z_tail,
            threads=1, accuracy_mode="fast", freq_grid="goal",
            physical_guard=True,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    spectrum = getattr(model, "log10OmegaGW", None)
    return {
        "returned": result is model,
        "input_DN_eff": float(dn_eff),
        "root_after_call": float(model.cosmo_param["DN_eff"]),
        "mapped_DN_gw": None if not hasattr(model, "DN_gw") else float(model.DN_gw[-1]),
        "failure_reason": model.fast_failure_reason,
        "failure_detail": getattr(model, "fast_failure_detail", None),
        "runtime_s": time.perf_counter() - started,
        "log10_OmegaGW": None if spectrum is None else np.asarray(spectrum, dtype=np.float64).tolist(),
    }


def reference_case(model, freqs, log_power, dn_eff, z_tail):
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        Ogw, Oj, Opgw, used_tail = spectrum_reference(
            model, freqs, dn_eff, z_tail=z_tail, rtol=1e-9,
            workers=1, backend="cartesian",
        )
    base = np.where(used_tail, 3.0 * Opgw, Ogw - Oj)
    if not np.isfinite(base).all() or np.any(base <= 0.0):
        raise FloatingPointError("Reference unit-transfer terminal spectrum is not positive finite")
    log_integrand = np.log(base) + log_power
    shift = float(np.max(log_integrand))
    normalized = np.exp(log_integrand - shift)
    from stiffgwpy_fast.reference import integrate_spectrum

    g2, qerr, ierr = integrate_spectrum(freqs, normalized, np.zeros_like(normalized))
    h = float(model.derived_param["h"])
    log_map = math.log(gp.Neff0 / (gp.Omega_nh2 / h**2)) + shift + math.log(g2)
    return {
        "z_tail": z_tail,
        "rtol": 1e-9,
        "runtime_s": time.perf_counter() - started,
        "mapped_DN_gw": math.exp(log_map),
        "log10_OmegaGW": (log_integrand / math.log(10.0)).tolist(),
        "used_tail_count": int(np.count_nonzero(used_tail)),
        "normalized_nodes_underflowed_to_zero": int(np.count_nonzero(normalized == 0.0)),
        "quadrature_error_relative": None if qerr is None else float(qerr / g2),
        "interpolation_error_relative": None if ierr is None else float(ierr / g2),
        "warnings": [str(item.message) for item in caught],
    }


def main() -> int:
    if sys.argv[1:] == ["--refresh-metadata-hashes"]:
        return refresh_metadata_hashes()
    if sys.argv[1:]:
        raise ValueError("unknown command-line option")
    model = LCDM_SG(**PARAMS)
    model.gen_expansion()
    log_amp = _log_tensor_amplitude(model.cosmo_param, model.derived_param)
    freqs = np.asarray(
        construct_fast_frequency_grid(model, log_amplitude=log_amp),
        dtype=np.float64,
    )
    start = {
        "schema_version": 1,
        "status": "fast_preflight_complete_reference_running",
        "scope": "New kappa decade (10x the prior low-T/high-tilt holdout): direct current Fast root, same-root maps, and direct Reference z8/z10 on the exact Fast grid.",
        "parameters": PARAMS,
        "fast_options": {"tol": 1e-7, "threads": 1, "physical_guard": True},
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": grid_hash(freqs),
        "log_A_t": log_amp,
        "stored_A_t": float(model.derived_param["A_t"]),
        "fast_preflight": None,
        "fast_fixed_maps": None,
        "reference": {},
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
    }
    root_run = run_fast(freqs, dn_eff=0.0, z_tail=5.0, tol=1e-7)
    start["fast_preflight"] = root_run
    if not root_run["returned"]:
        start["status"] = "preflight_rejected_fast_failure"
        write_checkpoint(start)
        print(json.dumps({"output": str(OUTPUT), "fast_preflight": root_run}, indent=2))
        return 0
    root = float(root_run["root_after_call"])
    fast_fixed = {
        "z5": run_fast(freqs, dn_eff=root, z_tail=5.0, tol=0.99),
        "z8": run_fast(freqs, dn_eff=root, z_tail=8.0, tol=0.99),
    }
    start["frozen_DN_eff"] = root
    start["fast_fixed_maps"] = fast_fixed
    start["status"] = "fast_preflight_complete_reference_running"
    write_checkpoint(start)

    unit_model = LCDM_SG(**{**PARAMS, "r": 1.0, "A_s": 1.0, "n_t": 0.0, "DN_eff": 0.0})
    unit_model.gen_expansion()
    log_power = log_amp + float(PARAMS["n_t"]) * (
        math.log(10.0) * freqs - math.log(gp.f_piv)
    )
    start["status"] = "reference_z8_running"
    write_checkpoint(start)
    ref8 = reference_case(unit_model, freqs, log_power, root, 8.0)
    start["reference"]["z8"] = ref8
    start["status"] = "reference_z10_running"
    write_checkpoint(start)
    ref10 = reference_case(unit_model, freqs, log_power, root, 10.0)
    start["reference"]["z10"] = ref10
    current_z8 = np.asarray(fast_fixed["z8"]["log10_OmegaGW"], dtype=np.float64)
    ref8_spectrum = np.asarray(ref8["log10_OmegaGW"], dtype=np.float64)
    ref10_spectrum = np.asarray(ref10["log10_OmegaGW"], dtype=np.float64)
    delta10 = current_z8 - ref10_spectrum
    start["tail_comparison"] = {
        "reference_z8_to_z10_max_abs_dex": float(np.max(np.abs(ref8_spectrum - ref10_spectrum))),
        "reference_z8_to_z10_p95_abs_dex": float(np.percentile(np.abs(ref8_spectrum - ref10_spectrum), 95.0)),
        "reference_z8_to_z10_relative_map_change": ref8["mapped_DN_gw"] / ref10["mapped_DN_gw"] - 1.0,
        "fast_z8_to_reference_z10_max_abs_dex": float(np.max(np.abs(delta10))),
        "fast_z8_to_reference_z10_p95_abs_dex": float(np.percentile(np.abs(delta10), 95.0)),
        "fast_z8_to_reference_z10_worst_log10f": float(freqs[int(np.argmax(np.abs(delta10)))]),
        "fast_z8_to_reference_z10_passes_1e-3_dex": bool(np.max(np.abs(delta10)) <= 1e-3),
        "fast_z8_map_relative_error_vs_reference_z10": fast_fixed["z8"]["mapped_DN_gw"] / ref10["mapped_DN_gw"] - 1.0,
        "fixed_point_residual_relative": fast_fixed["z5"]["mapped_DN_gw"] / root - 1.0,
    }
    start["status"] = "complete_current_fast_kappa_holdout"
    write_checkpoint(start)
    print(json.dumps({
        "output": str(OUTPUT),
        "frequency_count": freqs.size,
        "root": root,
        "fixed_point_residual": start["tail_comparison"]["fixed_point_residual_relative"],
        "reference_z8_to_z10_map": start["tail_comparison"]["reference_z8_to_z10_relative_map_change"],
        "fast_z8_vs_z10_max": start["tail_comparison"]["fast_z8_to_reference_z10_max_abs_dex"],
        "fast_z8_vs_z10_p95": start["tail_comparison"]["fast_z8_to_reference_z10_p95_abs_dex"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
