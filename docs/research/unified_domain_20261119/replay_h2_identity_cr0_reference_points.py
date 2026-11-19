"""Re-run Fast at archived CR0 Reference inputs and verify exact lineage."""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "replay_h2_identity_cr0_reference_points_current_source.json"
ARTIFACTS = (
    "validate_h2_identity_region_z8.json",
    "validate_h2_identity_density_slice_z8.json",
    "validate_h2_identity_control_z8.json",
    "validate_h2_identity_resource_boundary_z8.json",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    rows = []
    config = resolve_config("fast", threads=2)
    for filename in ARTIFACTS:
        artifact = json.loads((HERE / filename).read_text(encoding="utf-8"))
        for case in artifact.get("cases", []):
            if not case.get("fast_success"):
                continue
            model = LCDM_SG(**case["params"])
            started = time.perf_counter()
            result = SGWB_iter_fast(
                model, config=config, tol=1e-6, physical_guard=True
            )
            elapsed = time.perf_counter() - started
            raw = case["raw"]
            actual_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
            saved_log = np.asarray(raw["fast_log10_OmegaGW"], dtype=np.float64)
            actual_freqs = np.asarray(model.f, dtype=np.float64)
            saved_freqs = np.asarray(raw["freqs_log10_hz"], dtype=np.float64)
            rows.append({
                "artifact": filename,
                "case_id": case["case_id"],
                "success": result is model,
                "replay_seconds": elapsed,
                "root_exact": float(model.DN_gw[-1]) == float(case["fast_root"]),
                "frequency_grid_exact": bool(np.array_equal(actual_freqs, saved_freqs)),
                "spectrum_exact": bool(np.array_equal(actual_log, saved_log)),
                "max_abs_spectrum_delta": float(np.max(np.abs(actual_log - saved_log))),
            })
            print(rows[-1]["artifact"], rows[-1]["case_id"],
                  rows[-1]["replay_seconds"], rows[-1]["root_exact"],
                  rows[-1]["spectrum_exact"], flush=True)
    payload = {
        "status": "complete",
        "base_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "replay_config": {"h": config.h, "threads": 2, "z_tail": config.z_tail,
                          "frequency_grid": "goal", "tol": 1e-6},
        "source_sha256": {
            rel: _sha256(HERE.parents[2] / rel)
            for rel in ("stiffgwpy_fast/exact_background.py", "stiffgwpy_fast/fast_sgwb.py")
        },
        "rows": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "rows": len(rows),
        "all_success": all(row["success"] for row in rows),
        "all_roots_exact": all(row["root_exact"] for row in rows),
        "all_grids_exact": all(row["frequency_grid_exact"] for row in rows),
        "all_spectra_exact": all(row["spectrum_exact"] for row in rows),
        "artifact": str(OUTPUT),
    }))


if __name__ == "__main__":
    main()
