"""Recheck adaptive WKB handoff on a current cr>0 normal-amplitude control."""

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

INPUT = ROOT / "docs" / "research" / "results" / "underflow_tail_depth_consistency_branch_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_normal_control_handoff_20261116.json"
VARIANTS = handoff.VARIANTS
ORDERS = handoff.ORDERS
REPEATS = handoff.REPEATS


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    params = {key: float(value) for key, value in source["parameters"].items()}
    saved_root = float(source["fast_root_DN_eff"])
    ref8 = source["reference_fixed_root_physical_amplitude"]["8"]
    ref10 = source["reference_fixed_root_physical_amplitude"]["10"]
    ref10_spectrum = np.asarray(ref10["log10_OmegaGW"], dtype=np.float64)
    current_reference_sha = sha(ROOT / "stiffgwpy_fast" / "reference.py")
    if source["source_sha256"]["reference.py"] != current_reference_sha:
        raise RuntimeError("saved direct Reference source differs from the current implementation")

    root_model = LCDM_SG(**params)
    root_result = root_model.SGWB_iter(
        engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
        threads=1, accuracy_mode="fast", freq_grid="goal", physical_guard=True,
    )
    if root_result is not root_model:
        raise RuntimeError(f"current Fast normal-amplitude control failed: {root_model.fast_failure_reason}")
    freqs = np.asarray(root_model.f, dtype=np.float64)
    if freqs.size != source["frequency_count"] or handoff.digest(freqs) != source["frequency_grid_sha256_f64le"]:
        raise RuntimeError("current Fast grid differs from the archived direct Reference grid")
    ref_spectrum = np.asarray(ref10["log10_OmegaGW"], dtype=np.float64)
    if ref_spectrum.shape != freqs.shape or not np.isfinite(ref_spectrum).all():
        raise RuntimeError("saved direct Reference z10 spectrum is incomplete")

    process = __import__("psutil").Process()
    first_calls = {
        variant: handoff.run_variant(variant, params, freqs, saved_root, process)
        for variant in VARIANTS
    }
    records = {variant: [] for variant in VARIANTS}
    for repeat in range(REPEATS):
        for variant in ORDERS[repeat % len(ORDERS)]:
            records[variant].append(
                handoff.run_variant(variant, params, freqs, saved_root, process)
            )

    summaries = {}
    for variant in VARIANTS:
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
    direct_tail_delta = ref10_spectrum - ref8_spectrum
    boundary = float(source["derived_nt"])
    output = {
        "schema_version": 1,
        "status": "complete_current_cr_positive_normal_control_handoff_diagnostic",
        "scope": "Current Fast z5/z8/per-mode handoff A/B on the same physical-amplitude cr>0 control, frozen root, exact archived 76-node grid and one thread. Direct Reference z8/z10 is reused only after source and grid hashes match the current pipeline.",
        "input_artifact": str(INPUT.relative_to(ROOT)).replace("\\", "/"),
        "parameters": params,
        "derived_nt_used_by_solver": boundary,
        "frozen_DN_eff": saved_root,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": handoff.digest(freqs),
        "current_fast_root_DN_eff": float(root_model.cosmo_param["DN_eff"]),
        "current_fast_root_relative_difference_vs_saved": float(root_model.cosmo_param["DN_eff"] / saved_root - 1.0),
        "direct_reference": {
            "source_sha256_matches_current": True,
            "fixed_DN_eff": saved_root,
            "z8_map_DN_gw": float(ref8["mapped_DN_gw"]),
            "z10_map_DN_gw": float(ref10["mapped_DN_gw"]),
            "z8_to_z10_map_relative_change": float(ref10["mapped_DN_gw"] / ref8["mapped_DN_gw"] - 1.0),
            "z8_to_z10_max_abs_dex": float(np.max(np.abs(direct_tail_delta))),
            "z8_to_z10_p95_abs_dex": float(np.percentile(np.abs(direct_tail_delta), 95.0)),
            "z8_warnings": ref8["warnings"],
            "z10_warnings": ref10["warnings"],
        },
        "adaptive_rule": {
            "formula": "log10(1 + exp(-z_today) + exp(-2*z_today))",
            "target_ripple_dex": handoff.RIPPLE_TARGET_DEX,
            "z_today_threshold": handoff.ripple_threshold(handoff.RIPPLE_TARGET_DEX),
            "defer_to_z_tail": handoff.DEFER_TO,
        },
        "runtime_protocol": {
            "threads": 1,
            "warmup_calls_per_variant": 1,
            "paired_repeats_per_variant": REPEATS,
            "order_cycle": [list(order) for order in ORDERS],
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
                "reference_input": "docs/research/results/underflow_tail_depth_consistency_branch_20261115.json",
                "handoff_runner": "scripts/diagnose_cr_positive_consistency_tail_handoff_20261116.py",
                "diagnostic_script": "scripts/diagnose_cr_positive_normal_control_handoff_20261116.py",
            }.items()
        },
        "spectra": {
            "reference_z10_log10_OmegaGW": ref_spectrum.tolist(),
            "steady_state_first_log10_OmegaGW": {
                name: np.asarray(records[name][0]["log10_OmegaGW"], dtype=np.float64).tolist()
                for name in VARIANTS
            },
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "actual_derived_nt": boundary,
        "adaptive_vs_z8_speedup_p50": output["adaptive_vs_z8_speedup_p50"],
        "variants": {
            name: {
                "p50_s": summary["p50_runtime_s"],
                "p95_s": summary["p95_runtime_s"],
                "max_dex": summary["max_abs_spectrum_error_vs_reference_z10_dex"],
                "p95_dex": summary["p95_abs_spectrum_error_vs_reference_z10_dex"],
            }
            for name, summary in summaries.items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
