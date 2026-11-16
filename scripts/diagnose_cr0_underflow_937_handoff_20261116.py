"""Test ripple-aware handoff on the direct 937-node underflow holdout."""

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

GRID_INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
HOLDOUT_INPUT = ROOT / "docs" / "research" / "results" / "integral_first_blue_tilt_holdout_20261116.json"
REFERENCE_INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_z10_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr0_underflow_937_handoff_20261116.json"
GRID_KEY = "937"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    grid_data = json.loads(GRID_INPUT.read_text(encoding="utf-8"))
    holdout = json.loads(HOLDOUT_INPUT.read_text(encoding="utf-8"))
    reference = json.loads(REFERENCE_INPUT.read_text(encoding="utf-8"))
    params = {key: float(value) for key, value in holdout["parameters"].items()}
    freqs = np.asarray(grid_data["frequency_grids"][GRID_KEY]["log10_hz"], dtype=np.float64)
    grid_sha = handoff.digest(freqs)
    expected_sha = grid_data["frequency_grids"][GRID_KEY]["sha256_f64le"]
    if freqs.size != 937 or grid_sha != expected_sha:
        raise RuntimeError("canonical direct 937-point grid hash mismatch")
    if grid_sha != holdout["frequency_grids"][GRID_KEY]["sha256_f64le"]:
        raise RuntimeError("current holdout 937-point grid differs from the canonical nested grid")
    if reference["source_sha256"]["reference.py"] != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
        raise RuntimeError("direct Reference source hash differs from the current Reference")

    ref_spectrum = np.asarray(reference["reference"]["log10_OmegaGW"], dtype=np.float64)
    if ref_spectrum.shape != freqs.shape or not np.isfinite(ref_spectrum).all():
        raise RuntimeError("direct Reference z10 spectrum is incomplete")
    root = float(reference["frozen_DN_eff"])

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
        raise RuntimeError(f"current Fast 937-point root failed: {root_model.fast_failure_reason}")
    current_root = float(root_model.cosmo_param["DN_eff"])
    current_root_vs_reference = current_root / root - 1.0

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

    ref8_spectrum = np.asarray(reference["reference_z8"]["log10_OmegaGW"], dtype=np.float64)
    direct_tail_delta = ref_spectrum - ref8_spectrum
    output = {
        "schema_version": 1,
        "status": "complete_current_cr0_underflow_937_handoff_diagnostic",
        "scope": "Current Fast z5/z8/per-mode handoff A/B on the direct 937-node minimum-subnormal cr=0 holdout. The exact canonical nested grid, archived direct Reference z8/z10 source, current Fast root, and fixed root/grid/thread are checked independently.",
        "input_artifacts": {
            "grid": str(GRID_INPUT.relative_to(ROOT)).replace("\\", "/"),
            "current_fast_holdout": str(HOLDOUT_INPUT.relative_to(ROOT)).replace("\\", "/"),
            "direct_reference": str(REFERENCE_INPUT.relative_to(ROOT)).replace("\\", "/"),
        },
        "parameters": params,
        "derived_nt_used_by_solver": float(root_model.derived_param["nt"]),
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": grid_sha,
        "stored_A_t": float(root_model.derived_param["A_t"]),
        "log_A_t_finite": bool(np.isfinite(handoff.FA._log_tensor_amplitude(params, root_model.derived_param))),
        "reference": {
            "source_sha256_matches_current": True,
            "frozen_DN_eff": root,
            "current_fast_self_consistent_root": current_root,
            "current_fast_root_relative_difference_vs_reference": current_root_vs_reference,
            "z8_map_DN_gw": float(reference["reference_z8"]["mapped_DN_gw"]),
            "z10_map_DN_gw": float(reference["reference"]["mapped_DN_gw"]),
            "z8_to_z10_map_relative_change": float(reference["reference"]["mapped_DN_gw"] / reference["reference_z8"]["mapped_DN_gw"] - 1.0),
            "z8_to_z10_max_abs_dex": float(np.max(np.abs(direct_tail_delta))),
            "z8_to_z10_p95_abs_dex": float(np.percentile(np.abs(direct_tail_delta), 95.0)),
            "z8_warnings": reference["reference_z8"]["runtime_warnings"],
            "z10_warnings": reference["reference"]["runtime_warnings"],
            "z10_normalized_nodes_underflowed_to_zero": int(reference["reference"]["normalized_nodes_underflowed_to_zero"]),
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
                "grid_input": "docs/research/results/underflow_grid_235_469_937_20261115.json",
                "current_fast_holdout_input": "docs/research/results/integral_first_blue_tilt_holdout_20261116.json",
                "direct_reference_input": "docs/research/results/underflow_min_subnormal_reference_z10_20261115.json",
                "handoff_runner": "scripts/diagnose_cr_positive_consistency_tail_handoff_20261116.py",
                "diagnostic_script": "scripts/diagnose_cr0_underflow_937_handoff_20261116.py",
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
        "current_root_vs_reference": current_root_vs_reference,
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
