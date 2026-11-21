"""Verify current Fast support exactly matches the saved Reference inputs."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HERE = Path(__file__).resolve().parent
INPUTS = (
    "validate_adaptive_frequency_cr0_density_reference.json",
    "validate_adaptive_frequency_holdout_reference.json",
)
OUTPUT = HERE / "verify_adaptive_reference_reuse.json"


def main() -> None:
    rows = []
    for filename in INPUTS:
        saved = json.loads((HERE / filename).read_text(encoding="utf-8"))
        expected = saved["case"]
        model = LCDM_SG(**expected["params"])
        result = SGWB_iter_fast(
            model,
            config=resolve_config("fast", threads=2, freq_grid="adaptive"),
            tol=1.0e-6,
            physical_guard=True,
            freq_grid_target=1.0e-6,
            freq_grid_max_points=1500,
        )
        frequencies_equal = np.array_equal(
            np.asarray(model.f), np.asarray(expected["frequencies"])
        )
        root_equal = float(model.DN_gw[-1]) == float(expected["fast_root"])
        total_equal = float(model.cosmo_param["DN_eff"]) == float(
            expected["fast_total_DN_eff"]
        )
        rows.append(
            {
                "reference_file": filename,
                "fast_success": result is model,
                "frequency_support_equal": bool(frequencies_equal),
                "fast_root_equal": root_equal,
                "fast_total_DN_eff_equal": total_equal,
                "frequency_count": int(model.f.size),
                "reference_status": expected["reference_status"],
                "reference_source_sha256": saved["source_sha256"]["reference.py"],
            }
        )

    payload = {
        "status": "complete",
        "reference_rerun": False,
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
        "checks": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
