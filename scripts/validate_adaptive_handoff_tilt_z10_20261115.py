"""Run the z8/z10 per-mode validation on the independent n_t=17.30 holdout."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import validate_adaptive_handoff_background_z10_20261115 as validation  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_tilt_holdout_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_tilt_z10_20261115.json"
Z8_CHECKPOINT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_tilt_reference_z8_20261115.json"


def main() -> int:
    validation.INPUT = INPUT
    validation.OUTPUT = OUTPUT
    validation.REFERENCE_Z8_CHECKPOINT = Z8_CHECKPOINT
    validation.main()
    data = json.loads(OUTPUT.read_text(encoding="utf-8"))
    data["status"] = "complete_adaptive_handoff_tilt_holdout_reference_z10"
    data["scope"] = "Independent n_t=17.30, default T_re/kappa10, minimum-subnormal r holdout with direct Reference z8 and z10 tails on the saved native grid."
    data["holdout_label"] = "n_t=17.30, r=5e-324, T_re=2000, kappa10=.01"
    data["source_sha256"]["input_json"] = hashlib.sha256(INPUT.read_bytes()).hexdigest()
    data["source_sha256"]["input_json_path"] = INPUT.relative_to(ROOT).as_posix()
    data["source_sha256"]["holdout_wrapper"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    OUTPUT.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
