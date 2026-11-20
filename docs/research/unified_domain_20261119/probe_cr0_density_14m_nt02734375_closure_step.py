"""Check a second Reference closure map at the -0.02734375 14M corner."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HERE = Path(__file__).resolve().parent
PROBE = json.loads((HERE / "probe_cr0_density_14m_nt02734375_reference.json").read_text(encoding="utf-8"))
CASE = next(row for row in PROBE["cases"] if row["case_id"] == "upper_density_edge")
PARAMS = CASE["params"]
DN_EFF = float(CASE["reference_mapped_root"])
OUTPUT = HERE / "probe_cr0_density_14m_nt02734375_closure_step.json"


def main() -> None:
    model = LCDM_SG(**PARAMS)
    result = SGWB_iter_fast(model, config=resolve_config("fast", threads=2), tol=1.0e-6, physical_guard=True)
    if result is not model:
        raise RuntimeError(f"Fast unexpectedly failed: {model.fast_failure_reason}")
    started = time.perf_counter()
    try:
        ref = run_reference(
            model, dn_eff=DN_EFF, freq_subset=np.asarray(model.f, dtype=np.float64),
            z_tail=8.0, rtol=1.0e-9, self_consistent=False, workers=2,
        )
        outcome = {
            "status": "success", "mapped_dn_gw": float(ref["DN_gw"]),
            "seconds": time.perf_counter() - started,
        }
    except Exception as exc:
        outcome = {
            "status": "reference_integration_failure", "exception_type": type(exc).__name__,
            "message": str(exc), "seconds": time.perf_counter() - started,
        }
    payload = {
        "status": "complete_extra_closure_step", "params": PARAMS,
        "input_dn_eff": DN_EFF,
        "first_map_root_residual": CASE["root_relative_error"],
        "reference_settings": {"z_tail": 8.0, "rtol": 1.0e-9, "self_consistent": False, "workers": 2},
        "reference_code_modified": False, "closure_step": outcome,
        "source_sha256": {
            name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
            for name in ("fast_sgwb.py", "reference.py", "exact_background.py", "_wkb_transfer.py")
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
