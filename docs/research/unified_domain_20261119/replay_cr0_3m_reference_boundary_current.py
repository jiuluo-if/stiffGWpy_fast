"""Replay the accepted 3M CR0 Reference anchors on the current Fast source."""
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
fast = importlib.import_module("stiffgwpy_fast.fast_sgwb")
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG
SOURCE_FILES = ("exact_background.py", "fast_sgwb.py", "reference.py")
INPUT_FILES = (
    ("validate_cr0_3m_wkb_guard_z8_multi.json", ("upper_density_zero_tilt", "upper_density_negative_tilt", "independent_sobol_upper_density")),
    ("validate_cr0_high_density_extension_z8.json", ("upper_density_zero_tilt_corner", "upper_density_negative_tilt_corner")),
)


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values, dtype=np.float64).tobytes()).hexdigest()


def main() -> None:
    rows = []
    config = fast.resolve_config("fast", threads=2)
    for filename, names in INPUT_FILES:
        artifact = json.loads((HERE / filename).read_text(encoding="utf-8"))
        for saved in artifact["cases"]:
            if saved["case_id"] not in names:
                continue
            params = saved["params"]
            model = LCDM_SG(**params)
            solved = fast.SGWB_iter_fast(model, config=config, tol=1e-6, physical_guard=True)
            if solved is not model:
                rows.append({"case_id": saved["case_id"], "success": False, "reason": model.fast_failure_reason})
                continue
            old_f = np.asarray(saved["raw"]["freqs_log10_hz"], dtype=np.float64)
            old_spec = np.asarray(saved["raw"]["fast_log10_OmegaGW"], dtype=np.float64)
            new_f = np.asarray(model.f, dtype=np.float64)
            new_spec = np.asarray(model.log10OmegaGW, dtype=np.float64)
            rows.append({
                "case_id": saved["case_id"], "params": params, "success": True,
                "saved_fast_source_sha256": artifact["source_sha256"]["fast_sgwb.py"],
                "saved_reference_source_sha256": artifact["source_sha256"]["reference.py"],
                "current_root": float(model.DN_gw[-1]), "saved_root": float(saved["fast_root"]),
                "root_bitwise_equal": float(model.DN_gw[-1]) == float(saved["fast_root"]),
                "frequency_sha256_current": digest(new_f),
                "frequency_sha256_saved": digest(old_f), "frequency_bitwise_equal": np.array_equal(new_f, old_f),
                "spectrum_sha256_current": digest(new_spec),
                "spectrum_sha256_saved": digest(old_spec), "spectrum_bitwise_equal": np.array_equal(new_spec, old_spec),
                "max_spectrum_absdiff_dex": float(np.max(np.abs(new_spec-old_spec))) if new_spec.shape == old_spec.shape else None,
            })
    source_hashes = {
        name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
        for name in SOURCE_FILES
    }
    payload = {
        "status": "complete_current_fast_replay",
        "source_base_commit": "59a23aec2300f08b50124d7f35dec4b228388d37",
        "source_sha256": source_hashes,
        "reference_source_sha256": sorted({r["saved_reference_source_sha256"] for r in rows}),
        "cases": rows,
    }
    out = HERE / "replay_cr0_3m_reference_boundary_current.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"cases": len(rows), "all_exact": all(r.get("root_bitwise_equal") and r.get("frequency_bitwise_equal") and r.get("spectrum_bitwise_equal") for r in rows), "output": str(out)}, indent=2))


if __name__ == "__main__":
    main()
