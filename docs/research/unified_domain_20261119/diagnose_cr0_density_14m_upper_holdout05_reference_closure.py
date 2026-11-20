"""Recheck the extra Reference fixed-point closure map that failed once."""
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
PROBE = json.loads(
    (HERE / "probe_cr0_density_14m_upper_holdout05_reference_map.json").read_text(encoding="utf-8")
)
PARAMS = PROBE["case"]["params"]
DN_EFF = float(PROBE["case"]["reference_mapped_root"])
OUTPUT = HERE / "diagnose_cr0_density_14m_upper_holdout05_reference_closure.json"


def main() -> None:
    model = LCDM_SG(**PARAMS)
    solved = SGWB_iter_fast(model, config=resolve_config("fast", threads=2), tol=1.0e-6, physical_guard=True)
    if solved is not model:
        raise RuntimeError(f"Fast unexpectedly failed: {model.fast_failure_reason}")
    started = time.perf_counter()
    try:
        result = run_reference(
            model, dn_eff=DN_EFF, freq_subset=np.asarray(model.f, dtype=np.float64),
            z_tail=8.0, rtol=1.0e-9, self_consistent=False, workers=2,
        )
        outcome = {
            "status": "success", "mapped_dn_gw": float(result["DN_gw"]),
            "seconds": time.perf_counter() - started,
        }
    except Exception as exc:
        outcome = {
            "status": "reference_integration_failure", "exception_type": type(exc).__name__,
            "message": str(exc), "seconds": time.perf_counter() - started,
        }
    output = {
        "status": "reproduced_reference_closure_step_diagnostic",
        "input_params": PARAMS, "dn_eff": DN_EFF,
        "reference_settings": {"z_tail": 8.0, "rtol": 1.0e-9, "self_consistent": False, "workers": 2},
        "prior_fast_root_map": {
            "root_relative_error": PROBE["case"]["root_relative_error_at_fast_input"],
            "fixed_point_residual": PROBE["case"]["fixed_point_residual_at_fast_input"],
            "spectrum_max_abs_dex": PROBE["case"]["spectrum_max_abs_dex"],
        },
        "closure_step": outcome, "reference_code_modified": False,
        "source_sha256": {
            name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
            for name in ("fast_sgwb.py", "reference.py", "exact_background.py", "_wkb_transfer.py")
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
