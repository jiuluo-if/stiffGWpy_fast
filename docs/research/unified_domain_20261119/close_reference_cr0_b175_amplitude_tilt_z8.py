"""Third-map z8 Reference closure for CR0 amplitude/tilt representative points."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from stiffgwpy_fast.reference import run_reference
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

HERE = Path(__file__).resolve().parent
FIRST = HERE / "validate_cr0_b175_amplitude_tilt_z8.json"
OUTPUT = HERE / "close_reference_cr0_b175_amplitude_tilt_z8.json"
CASE_IDS = (
    "low_amplitude_high_tilt",
    "high_amplitude_high_tilt",
    "independent_sobol_holdout",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    artifact = json.loads(FIRST.read_text(encoding="utf-8"))
    source_hashes = {
        name: sha256(HERE.parents[2] / "stiffgwpy_fast" / name)
        for name in ("exact_background.py", "fast_sgwb.py", "reference.py")
    }
    rows = []
    for case_id in CASE_IDS:
        case = next(row for row in artifact["cases"] if row["case_id"] == case_id)
        model = LCDM_SG(**case["params"])
        freqs = np.asarray(case["raw"]["freqs_log10_hz"], dtype=np.float64)
        first_map = float(case["reference_map_at_fast_root"])
        maps, times = [], []
        current = first_map
        for _ in range(2):
            started = time.perf_counter()
            result = run_reference(
                model, dn_eff=current, freq_subset=freqs, z_tail=8.0,
                rtol=1e-9, self_consistent=False, workers=2,
            )
            times.append(time.perf_counter() - started)
            current = float(result["DN_gw"])
            maps.append(current)
        second_map, third_map = maps
        fast_root = float(case["fast_root"])
        rows.append({
            "case_id": case_id,
            "fast_root": fast_root,
            "first_reference_map_at_fast_root": first_map,
            "second_reference_map": second_map,
            "third_reference_map": third_map,
            "first_to_second_closure_relative": abs(second_map-first_map)
            / max(abs(second_map), 1e-300),
            "second_to_third_closure_relative": abs(third_map-second_map)
            / max(abs(third_map), 1e-300),
            "fast_root_relative_to_third_map": abs(fast_root-third_map)
            / max(abs(third_map), 1e-300),
            "reference_map_seconds": times,
        })
        print(case_id, rows[-1]["second_to_third_closure_relative"],
              rows[-1]["fast_root_relative_to_third_map"], flush=True)
    payload = {
        "status": "complete_cr0_amplitude_tilt_reference_map_closure",
        "reference_code_modified": False,
        "z_tail": 8.0,
        "same_saved_frequency_grids": True,
        "source_sha256": source_hashes,
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
