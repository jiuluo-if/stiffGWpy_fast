"""Repeat the pinned old Fast at its unstable high-r holdout point."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("cr1_old_repeatability_runner", HERE / "screen_h2_identity_region.py")
if spec is None or spec.loader is None:
    raise RuntimeError("could not load pinned-old-Fast replay helper")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)
SOURCE_SCREEN = HERE / "screen_cr1_r_upperband_candidate_lattice.json"
OUTPUT = HERE / "audit_cr1_r_upperband_old_repeatability_20.json"
REPEATS = 20


def main() -> None:
    paired = json.loads(SOURCE_SCREEN.read_text(encoding="utf-8"))
    case = next(row for row in paired["cases"]
                if row["input"]["case_id"] == "holdout_00")
    params = case["input"]["params"]
    cases = [{"case_id":f"repeat_{i:02d}","params":params} for i in range(REPEATS)]
    rows = screen._run_baseline(cases)
    successes = sum(bool(row.get("success")) for row in rows)
    payload = {
        "status":"complete_old_fast_same_input_repeatability_audit",
        "baseline_commit":screen.BASE_SHA,
        "case_id":"known_old_success_high_r_holdout",
        "params":params,
        "repeats":REPEATS,
        "successes":successes,
        "failures":REPEATS-successes,
        "source_screen_old_outcome":case["baseline_fast"],
        "rows":rows,
        "source_sha256":{
            "screen_h2_identity_region.py":hashlib.sha256((HERE/"screen_h2_identity_region.py").read_bytes()).hexdigest(),
        },
        "interpretation":"This is a point-level repeatability audit, not a domain success rate. Any repeated old failure remains counted separately from physical invalidity.",
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:payload[k] for k in ("baseline_commit","params","repeats","successes","failures")},indent=2))


if __name__ == "__main__":
    main()
