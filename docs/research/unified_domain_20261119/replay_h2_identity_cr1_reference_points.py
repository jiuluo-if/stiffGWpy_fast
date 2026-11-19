"""Current-source Fast replay for the cr=1 z8 Reference inputs."""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "validate_h2_identity_cr1_density_z8.json"
OUTPUT = HERE / "replay_h2_identity_cr1_reference_points_current_source.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    artifact = json.loads(SOURCE.read_text(encoding="utf-8"))
    config = resolve_config("fast", threads=2)
    rows = []
    for case in artifact["cases"]:
        model = LCDM_SG(**case["params"])
        started = time.perf_counter()
        result = SGWB_iter_fast(
            model, config=config, tol=1e-6, physical_guard=True
        )
        elapsed = time.perf_counter() - started
        raw = case["raw"]
        log_spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
        stored_log = np.asarray(raw["fast_log10_OmegaGW"], dtype=np.float64)
        row = {
            "case_id": case["case_id"],
            "success": result is model,
            "seconds": elapsed,
            "root": float(model.DN_gw[-1]),
            "stored_root": float(case["fast_root"]),
            "root_exact": float(model.DN_gw[-1]) == float(case["fast_root"]),
            "frequency_exact": bool(np.array_equal(model.f, raw["freqs_log10_hz"])),
            "spectrum_exact": bool(np.array_equal(log_spectrum, stored_log)),
            "max_spectrum_delta": float(np.max(np.abs(log_spectrum - stored_log))),
            "root_reference_proxy": float(case["root_relative_error_proxy"]),
            "spectrum_reference_max_dex": float(case["spectrum_max_abs_dex"]),
        }
        rows.append(row)
        print(row["case_id"], row["seconds"], row["root_exact"],
              row["spectrum_exact"], flush=True)
    payload = {
        "status": "complete",
        "reference_artifact": SOURCE.name,
        "reference_code_modified": False,
        "source_sha256": {
            rel: _sha256(HERE.parents[2] / rel)
            for rel in ("stiffgwpy_fast/exact_background.py", "stiffgwpy_fast/fast_sgwb.py")
        },
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
