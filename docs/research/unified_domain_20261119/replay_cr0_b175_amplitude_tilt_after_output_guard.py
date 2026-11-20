"""Replay the CR0 amplitude/tilt strip after the output-overflow guard change."""
from __future__ import annotations

import hashlib
import importlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
fast_sgwb = importlib.import_module("stiffgwpy_fast.fast_sgwb")
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG
SGWB_iter_fast = fast_sgwb.SGWB_iter_fast
resolve_config = fast_sgwb.resolve_config

INPUTS = (
    HERE / "screen_cr0_b175_amplitude_tilt_lattice.json",
    HERE / "screen_cr0_b175_amplitude_tilt_sobol.json",
)
REFERENCE = HERE / "validate_cr0_b175_amplitude_tilt_z8.json"
OUTPUT = HERE / "replay_cr0_b175_amplitude_tilt_after_output_guard.json"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    screen_rows = []
    for path in INPUTS:
        source = json.loads(path.read_text(encoding="utf-8"))
        for item in source["cases"]:
            params = item["input"]["params"]
            old = item["baseline_fast"]
            prior = item["candidate_fast"]
            model = LCDM_SG(**params)
            result = SGWB_iter_fast(
                model, config=resolve_config("fast", threads=2), tol=1e-6,
                physical_guard=True,
            )
            current_root = float(model.DN_gw[-1]) if result is model else None
            screen_rows.append({
                "case_id": item["input"]["case_id"],
                "old_failure_reason": old["failure_reason"],
                "prior_candidate_success": prior["success"],
                "current_success": result is model,
                "current_failure_reason": model.fast_failure_reason,
                "prior_root": prior["root"],
                "current_root": current_root,
                "root_abs_delta": (
                    None if current_root is None else abs(current_root-float(prior["root"]))
                ),
                "n_freq": int(len(model.f)) if result is model else None,
                "frequency_sha256": (
                    sha256_bytes(np.ascontiguousarray(model.f, dtype=np.float64).tobytes())
                    if result is model else None
                ),
                "spectrum_sha256": (
                    sha256_bytes(np.ascontiguousarray(model.log10OmegaGW, dtype=np.float64).tobytes())
                    if result is model else None
                ),
            })

    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    reference_rows = []
    for item in reference["cases"]:
        model = LCDM_SG(**item["params"])
        result = SGWB_iter_fast(
            model, config=resolve_config("fast", threads=2), tol=1e-6,
            physical_guard=True,
        )
        saved_freq = np.asarray(item["raw"]["freqs_log10_hz"], dtype=np.float64)
        saved_log = np.asarray(item["raw"]["fast_log10_OmegaGW"], dtype=np.float64)
        current_freq = np.asarray(model.f, dtype=np.float64) if result is model else np.array([])
        current_log = np.asarray(model.log10OmegaGW, dtype=np.float64) if result is model else np.array([])
        delta = float(np.max(np.abs(current_log-saved_log))) if result is model else None
        print(item["case_id"], bool(np.array_equal(current_freq, saved_freq)), delta,
              flush=True)
        reference_rows.append({
            "case_id": item["case_id"],
            "success": result is model,
            "frequency_exact": bool(np.array_equal(current_freq, saved_freq)),
            "spectrum_exact": bool(np.array_equal(current_log, saved_log)),
            "frequency_max_abs_delta": (
                float(np.max(np.abs(current_freq-saved_freq)))
                if result is model and current_freq.shape == saved_freq.shape else None
            ),
            "spectrum_max_abs_delta": delta,
            "root": float(model.DN_gw[-1]) if result is model else None,
        })

    payload = {
        "status": "complete_post_classification_guard_replay",
        "previous_domain_source_commit": "c3c6ab3fc34e705ebdea3d16739740a227d1acb7",
        "current_workspace_commit": "1cdb98b28cb780c031fea96c136677abd5fba079",
        "current_source_sha256": {
            name: file_sha(ROOT / "stiffgwpy_fast" / name)
            for name in ("exact_background.py", "fast_sgwb.py", "freq_adaptive.py", "_wkb_transfer.py", "reference.py")
        },
        "reference_code_modified": False,
        "screen_cases": screen_rows,
        "same_root_reference_cases": reference_rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "screen_cases": len(screen_rows),
        "current_successes": sum(row["current_success"] for row in screen_rows),
        "max_root_abs_delta": max(row["root_abs_delta"] or 0.0 for row in screen_rows),
        "reference_cases": len(reference_rows),
        "reference_exact_spectra": sum(row["spectrum_exact"] for row in reference_rows),
        "artifact": str(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
