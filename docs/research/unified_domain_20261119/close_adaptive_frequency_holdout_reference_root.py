"""Check Reference-map closure at the adaptive Fast CR0 holdout root."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HERE = Path(__file__).resolve().parent
INPUT = HERE / "validate_adaptive_frequency_holdout_reference.json"
OUTPUT = HERE / "close_adaptive_frequency_holdout_reference_root.json"


def main() -> None:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    case = source["case"]
    model = LCDM_SG(**case["params"])
    started = time.perf_counter()
    mapped_again = run_reference(
        model,
        dn_eff=float(case["reference_mapped_root"]),
        freq_subset=np.asarray(case["frequencies"], dtype=np.float64),
        z_tail=8.0,
        rtol=1.0e-9,
        self_consistent=False,
        workers=2,
    )
    seconds = time.perf_counter() - started
    second_root = float(mapped_again["DN_gw"])
    result = {
        "status": "complete",
        "reference_code_modified": False,
        "input_file": INPUT.name,
        "fast_root": float(case["fast_root"]),
        "reference_first_map": float(case["reference_mapped_root"]),
        "reference_second_map": second_root,
        "reference_closure_relative_error": abs(
            second_root - float(case["reference_mapped_root"])
        ) / max(abs(second_root), 1.0e-300),
        "adaptive_fast_root_relative_error_to_second_map": abs(
            second_root - float(case["fast_root"])
        ) / max(abs(second_root), 1.0e-300),
        "reference_seconds": seconds,
        "source_sha256": {
            name: hashlib.sha256(
                (ROOT / "stiffgwpy_fast" / name).read_bytes()
            ).hexdigest()
            for name in (
                "fast_sgwb.py",
                "reference.py",
                "exact_background.py",
                "_wkb_transfer.py",
            )
        },
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
