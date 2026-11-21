"""Third-map closure for the tight lower-temperature high-r corner."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
reference = __import__("stiffgwpy_fast.reference", fromlist=["run_reference"]).run_reference
FIRST = HERE / "validate_cr1_r_upperband_edge_probes_z8.json"
OUTPUT = HERE / "close_reference_cr1_r_upperband_z8.json"
CASE_ID = "T1875_high_r"


def main() -> None:
    artifact = json.loads(FIRST.read_text(encoding="utf-8"))
    case = next(row for row in artifact["cases"] if row["case_id"] == CASE_ID)
    model = __import__("stiffgwpy_fast.stiff_SGWB", fromlist=["LCDM_SG"]).LCDM_SG(**case["params"])
    frequencies = np.asarray(case["raw"]["freqs_log10_hz"], dtype=np.float64)
    first = float(case["reference_map_at_fast_root"])
    current = first
    maps, timings = [], []
    for _ in range(2):
        started = time.perf_counter()
        result = reference(model, dn_eff=current, freq_subset=frequencies,
                           z_tail=8.0, rtol=1e-9,
                           self_consistent=False, workers=2)
        timings.append(time.perf_counter()-started)
        current = float(result["DN_gw"])
        maps.append(current)
    second, third = maps
    fast_root = float(case["fast_root"])
    row = {
        "case_id":CASE_ID,"fast_root":fast_root,
        "first_reference_map_at_fast_root":first,
        "second_reference_map":second,"third_reference_map":third,
        "first_to_second_closure_relative":abs(second-first)/abs(second),
        "second_to_third_closure_relative":abs(third-second)/abs(third),
        "fast_root_relative_to_third_map":abs(fast_root-third)/abs(third),
        "reference_map_seconds":timings,
    }
    payload={
        "status":"complete_cr1_upperband_map_closure",
        "reference_code_modified":False,"z_tail":8.0,"same_saved_frequency_grid":True,
        "source_sha256":{name:hashlib.sha256((ROOT/"stiffgwpy_fast"/name).read_bytes()).hexdigest()
                         for name in ("exact_background.py","fast_sgwb.py","freq_adaptive.py","reference.py")},
        "case":row,
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(row,indent=2))


if __name__ == "__main__":
    main()
