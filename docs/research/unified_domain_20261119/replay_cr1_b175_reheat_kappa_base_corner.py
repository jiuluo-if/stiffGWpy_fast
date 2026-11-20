"""Recheck the shared CR1 base corner for exact reuse of its z8 Reference arrays."""
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
SOURCE = HERE / "validate_cr1_b175_square_z8.json"
OUTPUT = HERE / "replay_cr1_b175_reheat_kappa_base_corner.json"
CASE_ID = "high_square_corner"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha_array(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def main() -> None:
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = next(row for row in data["cases"] if row["case_id"] == CASE_ID)
    model = LCDM_SG(**case["params"])
    result = fast.SGWB_iter_fast(
        model, config=fast.resolve_config("fast", threads=2), tol=1e-6,
        physical_guard=True,
    )
    raw = case["raw"]
    old_freq = np.asarray(raw["freqs_log10_hz"], dtype=np.float64)
    old_log = np.asarray(raw["fast_log10_OmegaGW"], dtype=np.float64)
    current_freq = np.asarray(model.f, dtype=np.float64) if result is model else np.array([])
    current_log = np.asarray(model.log10OmegaGW, dtype=np.float64) if result is model else np.array([])
    payload = {
        "status": "complete_cr1_base_corner_replay",
        "case_id": CASE_ID,
        "params": case["params"],
        "success": result is model,
        "frequency_exact": bool(np.array_equal(current_freq, old_freq)),
        "spectrum_exact": bool(np.array_equal(current_log, old_log)),
        "root_exact": float(model.DN_gw[-1]) == float(case["fast_root"]),
        "max_spectrum_abs_delta": (
            float(np.max(np.abs(current_log-old_log)))
            if result is model and current_log.shape == old_log.shape else None
        ),
        "root": float(model.DN_gw[-1]) if result is model else None,
        "previous_reference_root_error": case["root_relative_error_proxy"],
        "previous_reference_spectrum_max_dex": case["spectrum_max_abs_dex"],
        "current_source_sha256": {
            name: sha(ROOT / "stiffgwpy_fast" / name)
            for name in ("exact_background.py", "fast_sgwb.py", "freq_adaptive.py", "reference.py")
        },
        "saved_frequency_sha256": sha_array(old_freq),
        "current_frequency_sha256": sha_array(current_freq) if result is model else None,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
