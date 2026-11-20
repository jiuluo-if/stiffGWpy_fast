"""Third z8 Reference maps for CR1 root-closure checks at b175f0b."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from stiffgwpy_fast.reference import run_reference
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

HERE = Path(__file__).resolve().parent
FIRST = HERE / "validate_cr1_b175_square_z8.json"
OUTPUT = HERE / "close_reference_cr1_b175_square_z8.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    first = json.loads(FIRST.read_text(encoding="utf-8"))
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
    for case in first["cases"]:
        model = LCDM_SG(**case["params"])
        freqs = np.asarray(case["raw"]["freqs_log10_hz"], dtype=np.float64)
        first_map = float(case["reference_map_at_fast_root"])
        old_row = prior.get(case["case_id"], {})
        if old_row.get("first_reference_map_at_fast_root") == first_map:
            second_map = float(old_row["second_reference_map"])
            first_to_second = old_row.get(
                "first_to_second_closure_relative",
                old_row.get("second_map_closure_relative"),
            )
            first_map_seconds = old_row.get(
                "first_map_seconds", old_row.get("second_map_seconds")
            )
        else:
            started = time.perf_counter()
            second = run_reference(
                model, dn_eff=first_map, freq_subset=freqs, z_tail=8.0,
                rtol=1e-9, self_consistent=False, workers=2,
            )
            first_map_seconds = time.perf_counter() - started
            second_map = float(second["DN_gw"])
            first_to_second = abs(second_map-first_map) / max(abs(second_map), 1e-300)
        started = time.perf_counter()
        third = run_reference(
            model, dn_eff=second_map, freq_subset=freqs, z_tail=8.0,
            rtol=1e-9, self_consistent=False, workers=2,
        )
        third_map_seconds = time.perf_counter() - started
        third_map = float(third["DN_gw"])
        fast_root = float(case["fast_root"])
        rows.append({
            "case_id": case["case_id"],
            "fast_root": fast_root,
            "first_reference_map_at_fast_root": first_map,
            "second_reference_map": second_map,
            "third_reference_map": third_map,
            "first_to_second_closure_relative": first_to_second,
            "second_to_third_closure_relative": abs(third_map-second_map)
            / max(abs(third_map), 1e-300),
            "fast_root_relative_to_third_map": abs(fast_root-third_map)
            / max(abs(third_map), 1e-300),
            "first_map_seconds": first_map_seconds,
            "third_map_seconds": third_map_seconds,
        })
        print(case["case_id"], rows[-1]["second_to_third_closure_relative"],
              rows[-1]["fast_root_relative_to_third_map"], flush=True)
    payload = {
        "status": "complete_cr1_second_reference_maps",
        "reference_code_modified": False,
        "z_tail": 8.0,
        "same_frequency_grid": True,
        "source_sha256": source_hashes,
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
