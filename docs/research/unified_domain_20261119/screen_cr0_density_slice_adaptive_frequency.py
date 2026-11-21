"""Paired adaptive/dense Fast screen on the prior CR0 density-tilt slice."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HERE = Path(__file__).resolve().parent
BASELINE = HERE / "screen_cr0_density_13p875m_to13p99609375m_nt0302734375_slice.json"
OUTPUT = HERE / "screen_cr0_density_slice_adaptive_frequency.json"


def main() -> None:
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    rows = []
    for sample in baseline["cases"]:
        params = sample["input"]["params"]
        old = sample["baseline_fast"]

        goal_model = LCDM_SG(**params)
        started = time.perf_counter()
        goal_result = SGWB_iter_fast(
            goal_model,
            config=resolve_config("fast", threads=2, freq_grid="goal"),
            tol=1.0e-6,
            physical_guard=True,
        )
        goal_seconds = time.perf_counter() - started

        adaptive_model = LCDM_SG(**params)
        started = time.perf_counter()
        adaptive_result = SGWB_iter_fast(
            adaptive_model,
            config=resolve_config("fast", threads=2, freq_grid="adaptive"),
            tol=1.0e-6,
            physical_guard=True,
            freq_grid_target=1.0e-6,
            freq_grid_max_points=1500,
        )
        adaptive_seconds = time.perf_counter() - started
        row = {
            "case_id": sample["input"]["case_id"],
            "family": sample["input"]["family"],
            "params": params,
            "old_fast_success": old["success"],
            "old_fast_failure_reason": old["failure_reason"],
            "current_goal_success": goal_result is goal_model,
            "current_goal_failure_reason": goal_model.fast_failure_reason,
            "current_goal_seconds": goal_seconds,
            "current_goal_frequency_count": int(goal_model.f.size)
            if goal_result is goal_model
            else None,
            "current_goal_root": float(goal_model.DN_gw[-1])
            if goal_result is goal_model
            else None,
            "adaptive_fast_success": adaptive_result is adaptive_model,
            "adaptive_fast_failure_reason": adaptive_model.fast_failure_reason,
            "adaptive_seconds": adaptive_seconds,
        }
        if adaptive_result is not adaptive_model:
            rows.append(row)
            continue

        row.update(
            {
                "adaptive_root": float(adaptive_model.DN_gw[-1]),
                "adaptive_root_relative_correction_vs_current_goal": (
                    abs(
                        float(adaptive_model.DN_gw[-1])
                        - float(goal_model.DN_gw[-1])
                    )
                    / max(abs(float(adaptive_model.DN_gw[-1])), 1.0e-300)
                    if goal_result is goal_model
                    else None
                ),
                "adaptive_total_DN_eff": float(adaptive_model.cosmo_param["DN_eff"]),
                "frequency_count": int(adaptive_model.f.size),
                "frequency_error_estimate": float(adaptive_model.freq_grid_error),
            }
        )
        dense_model = LCDM_SG(**params)
        started = time.perf_counter()
        dense_result = SGWB_iter_fast(
            dense_model,
            config=resolve_config("fast", threads=2, freq_grid="grid_independent"),
            tol=1.0e-6,
            physical_guard=True,
            freq_res=4.6,
            eval_freqs=np.asarray(adaptive_model.f, dtype=np.float64),
        )
        dense_seconds = time.perf_counter() - started
        row["dense_seconds"] = dense_seconds
        row["dense_fast_success"] = dense_result is dense_model
        row["dense_fast_failure_reason"] = dense_model.fast_failure_reason
        if dense_result is dense_model:
            dense_indices = np.searchsorted(-dense_model.f, -adaptive_model.f)
            common_nodes = np.array_equal(dense_model.f[dense_indices], adaptive_model.f)
            row["common_native_nodes"] = bool(common_nodes)
            row["root_relative_error_vs_dense_fast"] = abs(
                float(adaptive_model.DN_gw[-1]) - float(dense_model.DN_gw[-1])
            ) / max(abs(float(dense_model.DN_gw[-1])), 1.0e-300)
            valid = np.isfinite(adaptive_model.log10OmegaGW) & np.isfinite(
                dense_model.log10OmegaGW[dense_indices]
            )
            spectrum_error = np.abs(
                adaptive_model.log10OmegaGW[valid]
                - dense_model.log10OmegaGW[dense_indices][valid]
            )
            row["spectrum_valid_count"] = int(valid.sum())
            row["max_native_spectrum_error_dex"] = (
                float(np.max(spectrum_error)) if spectrum_error.size else None
            )
        rows.append(row)

    payload = {
        "status": "complete",
        "reference_code_modified": False,
        "baseline_file": BASELINE.name,
        "baseline_source_sha256": baseline.get("candidate_source_sha256"),
        "adaptive_settings": {"freq_grid_target": 1.0e-6, "freq_grid_max_points": 1500},
        "dense_fast_settings": {"freq_res": 4.6},
        "source_sha256": {
            name: hashlib.sha256(
                (ROOT / "stiffgwpy_fast" / name).read_bytes()
            ).hexdigest()
            for name in (
                "fast_sgwb.py",
                "reference.py",
                "exact_background.py",
                "_wkb_transfer.py",
            )
        },
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
