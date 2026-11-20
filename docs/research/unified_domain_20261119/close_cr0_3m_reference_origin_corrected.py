"""Correct the CR0 3M Reference closure using input DN_eff as fixed-point origin."""
from __future__ import annotations

import hashlib
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
reference = importlib.import_module("stiffgwpy_fast.reference").run_reference
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG
FAST_REPLAY = json.loads((HERE / "replay_cr0_3m_reference_boundary_current.json").read_text(encoding="utf-8"))
FIRST = json.loads((HERE / "validate_cr0_3m_wkb_guard_z8_multi.json").read_text(encoding="utf-8"))
OUTPUT = HERE / "close_cr0_3m_reference_origin_corrected.json"
CASE_IDS = (
    "upper_density_zero_tilt",
    "upper_density_negative_tilt",
    "independent_sobol_upper_density",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reference_total_from_map(input_dn_eff_origin: float, mapped_dn_gw: float) -> float:
    """Form the next total DN_eff from the original input plus Reference map."""
    return float(input_dn_eff_origin + mapped_dn_gw)


def main() -> None:
    rows = []
    for case_id in CASE_IDS:
        saved = next(r for r in FIRST["cases"] if r["case_id"] == case_id)
        replay = next(r for r in FAST_REPLAY["cases"] if r["case_id"] == case_id)
        if not all(replay[k] for k in ("root_bitwise_equal", "frequency_bitwise_equal", "spectrum_bitwise_equal")):
            raise RuntimeError(f"current Fast replay changed for {case_id}; cannot reuse saved Reference map")
        params = saved["params"]
        origin = float(params["DN_eff"])
        fast_root = float(saved["fast_root"])
        total_at_fast_root = _reference_total_from_map(origin, fast_root)
        map0 = float(saved["reference_map_history"][0]["mapped_dn_gw"])
        current = map0
        model = LCDM_SG(**params)
        frequencies = np.asarray(saved["raw"]["freqs_log10_hz"], dtype=np.float64)
        maps = [{
            "input_total_dn_eff": total_at_fast_root,
            "mapped_dn_gw": map0,
            "origin_correct": True,
            "source": "saved_reference_map_at_fast_root",
        }]
        timings = []
        for _ in range(2):
            total_dn_eff = _reference_total_from_map(origin, current)
            started = time.perf_counter()
            result = reference(
                model, dn_eff=total_dn_eff, freq_subset=frequencies,
                z_tail=8.0, rtol=1e-9, self_consistent=False, workers=2,
            )
            elapsed = time.perf_counter() - started
            mapped = float(result["DN_gw"])
            maps.append({"input_total_dn_eff": total_dn_eff, "mapped_dn_gw": mapped})
            timings.append(elapsed)
            current = mapped
        residual = abs(maps[-1]["mapped_dn_gw"] - maps[-2]["mapped_dn_gw"]) / max(abs(maps[-1]["mapped_dn_gw"]), 1e-300)
        root_error = abs(current - fast_root) / max(abs(current), 1e-300)
        rows.append({
            "case_id": case_id, "params": params, "input_dn_eff_origin": origin,
            "fast_root": fast_root, "fast_total_dn_eff": total_at_fast_root,
            "map_history": maps, "map_closure_relative": residual,
            "fast_root_relative_to_corrected_reference": root_error,
            "additional_reference_seconds": timings,
            "reference_source_sha256": FIRST["source_sha256"]["reference.py"],
            "current_fast_replay_exact": True,
        })
    payload = {
        "status": "complete_correct_origin_reference_closure",
        "method": "DN_eff_total = input_DN_eff_origin + DN_gw; iterate at origin + previous map",
        "reference_code_modified": False,
        "source_sha256": {
            name: sha256(ROOT / "stiffgwpy_fast" / name)
            for name in ("exact_background.py", "fast_sgwb.py", "_wkb_transfer.py", "reference.py")
        },
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    for row in rows:
        print(row["case_id"], row["fast_root_relative_to_corrected_reference"], row["map_closure_relative"], flush=True)
    print("saved", OUTPUT)


if __name__ == "__main__":
    main()
