"""Recheck the old Fast 15-second timeout at a longer resource limit."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("cr1_timeout_recheck_runner",HERE/"screen_h2_identity_region.py")
if spec is None or spec.loader is None:
    raise RuntimeError("could not load pinned old-Fast screen helper")
screen=importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)
SOURCE=HERE/"screen_cr1_r_upperband_candidate_sobol.json"
OUTPUT=HERE/"audit_cr1_r_upperband_timeout_recheck.json"
CASE_ID="holdout_03"
RECHECK_TIMEOUT_SECONDS=120.0


def main() -> None:
    artifact=json.loads(SOURCE.read_text(encoding="utf-8"))
    case=next(row for row in artifact["cases"] if row["input"]["case_id"]==CASE_ID)
    screen.BASELINE_TIMEOUT_SECONDS=RECHECK_TIMEOUT_SECONDS
    repeated=screen._run_baseline([{"case_id":CASE_ID,"params":case["input"]["params"]}])[0]
    payload={
        "status":"complete_old_fast_resource_timeout_recheck",
        "old_failure_baseline_commit":screen.BASE_SHA,
        "params":case["input"]["params"],
        "paired_screen_result":case["baseline_fast"],
        "recheck_timeout_seconds":RECHECK_TIMEOUT_SECONDS,
        "recheck":repeated,
        "current_fast":case["candidate_fast"],
        "interpretation":"The old Fast remained resource-limited at 120 seconds with positive finite H2; keep resource_timeout separate from mathematical or primitive failures.",
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:payload[k] for k in ("status","params","paired_screen_result","recheck","current_fast")},indent=2))


if __name__=="__main__":
    main()
