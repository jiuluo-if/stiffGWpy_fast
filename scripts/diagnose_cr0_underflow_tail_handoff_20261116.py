"""Recheck adaptive WKB handoff on a current cr=0 underflow holdout."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import diagnose_cr_positive_consistency_tail_handoff_20261116 as handoff  # noqa: E402

from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HOLDOUT_INPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_background_holdout_20261115.json"
REFERENCE_INPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_background_z10_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr0_underflow_tail_handoff_20261116.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    holdout = json.loads(HOLDOUT_INPUT.read_text(encoding="utf-8"))
    reference = json.loads(REFERENCE_INPUT.read_text(encoding="utf-8"))
    params = {key: float(value) for key, value in holdout["parameters"].items()}
    ref_sha = reference["source_sha256"]["reference.py"]
    if ref_sha != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
        raise RuntimeError("saved direct Reference source differs from the current implementation")
    freqs = np.asarray(holdout["frequency_log10_hz"], dtype=np.float64)
    if freqs.size != 221 or handoff.digest(freqs) != holdout["frequency_sha256_f64le"]:
        raise RuntimeError("saved direct Reference frequency grid changed")
    root = float(reference["frozen_DN_eff"])
    ref8 = reference["reference"]["8.0"]
    ref10 = reference["reference"]["10.0"]
    ref_spectrum = np.asarray(ref10["log10_OmegaGW"], dtype=np.float64)
    if ref_spectrum.shape != freqs.shape or not np.isfinite(ref_spectrum).all():
        raise RuntimeError("saved direct Reference z10 spectrum is incomplete")

    root_model = LCDM_SG(**params)
    original_goal = handoff.FA.goal_oriented_freqs
    original_construct = handoff.FA.construct_fast_frequency_grid
    handoff.FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    handoff.FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    try:
        root_result = root_model.SGWB_iter(
            engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
            threads=1, accuracy_mode="fast", freq_grid="goal", physical_guard=True,
        )
    finally:
        handoff.FA.goal_oriented_freqs = original_goal
        handoff.FA.construct_fast_frequency_grid = original_construct
    if root_result is not root_model:
        raise RuntimeError(f"current Fast underflow root failed: {root_model.fast_failure_reason}")
    if handoff.digest(np.asarray(root_model.f, dtype=np.float64)) != holdout["frequency_sha256_f64le"]:
        raise RuntimeError("current Fast goal grid differs from the direct Reference grid")
    current_root = float(root_model.cosmo_param["DN_eff"])
    root_relative_difference = current_root / root - 1.0

    process = __import__("psutil").Process()
    first_calls = {
        variant: handoff.run_variant(variant, params, freqs, root, process)
        for variant in handoff.VARIANTS
    }
    records = {variant: [] for variant in handoff.VARIANTS}
    for repeat in range(handoff.REPEATS):
        for variant in handoff.ORDERS[repeat % len(handoff.ORDERS)]:
            records[variant].append(handoff.run_variant(variant, params, freqs, root, process))

    summaries = {}
    for variant in handoff.VARIANTS:
        summary = handoff.summarize(records[variant], ref_spectrum, freqs)
        summary["max_abs_spectrum_error_vs_reference_z10_dex"] = summary.pop(
            "max_abs_spectrum_error_vs_reference_z12_dex"
        )
        summary["p95_abs_spectrum_error_vs_reference_z10_dex"] = summary.pop(
            "p95_abs_spectrum_error_vs_reference_z12_dex"
        )
        summary["passes_1e-3_dex"] = bool(
            summary["max_abs_spectrum_error_vs_reference_z10_dex"] <= 1e-3
        )
        summaries[variant] = summary

    ref8_spectrum = np.asarray(ref8["log10_OmegaGW"], dtype=np.float64)
    tail_delta = ref_spectrum - ref8_spectrum
    output = {
        "schema_version": 1,
        "status": "complete_current_cr0_underflow_handoff_diagnostic",
        "scope": "Current Fast z5/z8/per-mode handoff A/B on the same minimum-subnormal cr=0 holdout, frozen root, exact archived 221-node direct-Reference grid and one thread. The current native grid builder differs from this archived grid hash, so the exact saved direct grid is forced for current Fast root and fixed-root runs.",
        "holdout_input": str(HOLDOUT_INPUT.relative_to(ROOT)).replace("\\", "/"),
        "reference_input": str(REFERENCE_INPUT.relative_to(ROOT)).replace("\\", "/"),
        "parameters": params,
        "derived_nt_used_by_solver": float(LCDM_SG(**params).derived_param["nt"]),
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": handoff.digest(freqs),
        "stored_A_t": float(LCDM_SG(**params).derived_param["A_t"]),
        "log_A_t_finite": bool(np.isfinite(handoff.FA._log_tensor_amplitude(params, LCDM_SG(**params).derived_param))),
        "reference": {
            "source_sha256_matches_current": True,
            "frozen_DN_eff": root,
            "current_fast_root_DN_eff": current_root,
            "current_fast_root_relative_difference_vs_reference_root": root_relative_difference,
            "z8_map_DN_gw": float(ref8["mapped_DN_gw"]),
            "z10_map_DN_gw": float(ref10["mapped_DN_gw"]),
            "z8_to_z10_map_relative_change": float(ref10["mapped_DN_gw"] / ref8["mapped_DN_gw"] - 1.0),
            "z8_to_z10_max_abs_dex": float(np.max(np.abs(tail_delta))),
            "z8_to_z10_p95_abs_dex": float(np.percentile(np.abs(tail_delta), 95.0)),
            "z8_warnings": ref8["warnings"],
            "z10_warnings": ref10["warnings"],
            "z10_normalized_nodes_underflowed_to_zero": int(ref10["normalized_nodes_underflowed_to_zero"]),
        },
        "adaptive_rule": {
            "formula": "log10(1 + exp(-z_today) + exp(-2*z_today))",
            "target_ripple_dex": handoff.RIPPLE_TARGET_DEX,
            "z_today_threshold": handoff.ripple_threshold(handoff.RIPPLE_TARGET_DEX),
            "defer_to_z_tail": handoff.DEFER_TO,
            "phase_resolved_rows_remain_on_original_handoff": True,
        },
        "runtime_protocol": {
            "threads": 1,
            "warmup_calls_per_variant": 1,
            "paired_repeats_per_variant": handoff.REPEATS,
            "order_cycle": [list(order) for order in handoff.ORDERS],
            "module_import_time_excluded": True,
            "first_solver_call_may_include_jit": True,
        },
        "first_solver_call_after_module_import": {
            name: {key: value for key, value in record.items() if key not in ("log10_OmegaGW", "handoff_eps", "z_today")}
            for name, record in first_calls.items()
        },
        "variants": summaries,
        "adaptive_vs_z8_speedup_p50": summaries["z8"]["p50_runtime_s"] / summaries["adaptive"]["p50_runtime_s"],
        "adaptive_vs_z5_p50_runtime_ratio": summaries["adaptive"]["p50_runtime_s"] / summaries["z5"]["p50_runtime_s"],
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "numba": __import__("numba").__version__,
            "psutil": __import__("psutil").__version__,
        },
        "source_sha256": {
            name: sha(ROOT / path)
            for name, path in {
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "freq_adaptive.py": "stiffgwpy_fast/freq_adaptive.py",
                "stiff_SGWB.py": "stiffgwpy_fast/stiff_SGWB.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "holdout_input": "docs/research/results/adaptive_handoff_background_holdout_20261115.json",
                "reference_input": "docs/research/results/adaptive_handoff_background_z10_20261115.json",
                "handoff_runner": "scripts/diagnose_cr_positive_consistency_tail_handoff_20261116.py",
                "diagnostic_script": "scripts/diagnose_cr0_underflow_tail_handoff_20261116.py",
            }.items()
        },
        "spectra": {
            "reference_z10_log10_OmegaGW": ref_spectrum.tolist(),
            "steady_state_first_log10_OmegaGW": {
                name: np.asarray(records[name][0]["log10_OmegaGW"], dtype=np.float64).tolist()
                for name in handoff.VARIANTS
            },
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "current_root_vs_reference": root_relative_difference,
        "adaptive_vs_z8_speedup_p50": output["adaptive_vs_z8_speedup_p50"],
        "variants": {
            name: {
                "p50_s": summary["p50_runtime_s"],
                "p95_s": summary["p95_runtime_s"],
                "max_dex": summary["max_abs_spectrum_error_vs_reference_z10_dex"],
                "p95_dex": summary["p95_abs_spectrum_error_vs_reference_z10_dex"],
                "selected": summary["adaptive_selected_count"],
            }
            for name, summary in summaries.items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
