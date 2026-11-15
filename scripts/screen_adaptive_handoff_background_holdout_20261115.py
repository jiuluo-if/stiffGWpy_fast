"""Run the adaptive handoff screen on a high-temperature/high-kappa holdout."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import screen_adaptive_handoff_tilt_holdout_20261115 as screen  # noqa: E402

OUTPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_background_holdout_20261115.json"
REJECTED_OUTPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_physical_preflight_rejection_20261115.json"


def main() -> int:
    if OUTPUT.exists():
        rejected = json.loads(OUTPUT.read_text(encoding="utf-8"))
        rejected["status"] = "preflight_rejected_physical_guard_not_a_valid_holdout"
        rejected["rejection"] = {
            "reason": "Exploratory root with physical_guard=False was DN_eff=1.382750485e9, above the production physical guard threshold 5; discard its candidate accuracy metrics.",
            "used_for_acceptance": False,
        }
        rejected["source_sha256"]["background_holdout_wrapper"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        REJECTED_OUTPUT.write_text(json.dumps(rejected, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    params = {
        **screen.PARAMS,
        "n_t": 17.45,
        "T_re": 10.0,
        "kappa10": 1e-3,
    }
    screen.PARAMS = params
    screen.OUTPUT = OUTPUT
    screen.main()
    data = json.loads(OUTPUT.read_text(encoding="utf-8"))
    data["holdout_label"] = "n_t=17.45, T_re=10, kappa10=1e-3, r=5e-324"
    data["source_sha256"]["background_holdout_wrapper"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    OUTPUT.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
