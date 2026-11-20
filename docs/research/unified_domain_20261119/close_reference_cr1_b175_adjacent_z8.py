"""Close z8 Reference maps at the CR1 adjacent-box boundary and holdout."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from stiffgwpy_fast.reference import run_reference
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

HERE = Path(__file__).resolve().parent
INPUTS = (
    (HERE / "validate_cr1_b175_boundary_probe_z8.json", "boundary_diagonal_1500k"),
    (HERE / "validate_cr1_b175_boundary_probe_z8.json", "boundary_diagonal_1600k"),
    (HERE / "validate_cr1_b175_adjacent_sobol_z8.json", "adjacent_independent_sobol_upper_sum"),
)
OUTPUT = HERE / "close_reference_cr1_b175_adjacent_z8.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    source_hashes = {
        name: sha256(HERE.parents[2] / "stiffgwpy_fast" / name)
        for name in ("exact_background.py", "fast_sgwb.py", "reference.py")
    }
    prior = {}
    if OUTPUT.exists():
        saved = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if saved.get("source_sha256") == source_hashes:
            prior = {row["case_id"]: row for row in saved.get("cases", [])}
    rows = []
    for path, case_id in INPUTS:
        artifact = json.loads(path.read_text(encoding="utf-8"))
        case = next(row for row in artifact["cases"] if row["case_id"] == case_id)
        model = LCDM_SG(**case["params"])
        frequencies = np.asarray(case["raw"]["freqs_log10_hz"], dtype=np.float64)
        first_map = float(case["reference_map_at_fast_root"])
        old = prior.get(case_id, {})
        if old.get("first_reference_map_at_fast_root") == first_map:
            second_map = float(old["second_reference_map"])
            third_map = float(old["third_reference_map"])
            times = old["reference_map_seconds"]
        else:
            maps = []
            times = []
            next_input = first_map
            for _ in range(2):
                started = time.perf_counter()
                result = run_reference(
                    model, dn_eff=next_input, freq_subset=frequencies,
                    z_tail=8.0, rtol=1e-9, self_consistent=False, workers=2,
                )
                times.append(time.perf_counter() - started)
                next_input = float(result["DN_gw"])
                maps.append(next_input)
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
        "status": "complete_cr1_adjacent_reference_map_closure",
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
