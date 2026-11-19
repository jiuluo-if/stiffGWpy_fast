"""Fast-only replay against archived direct Reference checks."""

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
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HERE = Path(__file__).resolve().parent
INPUT_FILES = (
    "validate_h2_identity_cr1_density_z8.json",
    "validate_h2_identity_cr1_low_rectangle_z8.json",
    "validate_h2_identity_cr1_mid_rectangle_z8.json",
    "validate_h2_identity_cr1_upper_boundary_z8.json",
    "validate_h2_identity_cr1_600_square_z8.json",
    "validate_h2_identity_cr1_1300_boundary_z8.json",
)
OUTPUT = HERE / "replay_exact_quarter_candidate_reference_points.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    cases = []
    for filename in INPUT_FILES:
        artifact = json.loads((HERE / filename).read_text(encoding="utf-8"))
        for row in artifact["cases"]:
            if "root_relative_error_proxy" not in row:
                continue
            if row["root_relative_error_proxy"] > 2e-4:
                continue
            cases.append((filename, row))
    results = []
    for source_file, row in cases:
        model = LCDM_SG(**row["params"])
        start = time.perf_counter()
        result = SGWB_iter_fast(
            model, config=resolve_config("fast", threads=2), tol=1e-6,
            physical_guard=True,
        )
        seconds = time.perf_counter()-start
        if result is not model:
            results.append({"source":source_file,"case_id":row["case_id"],"success":False,"failure_reason":model.fast_failure_reason,"seconds":seconds})
            continue
        raw=row["raw"]
        f=np.asarray(model.f,dtype=np.float64);rf=np.asarray(raw["freqs_log10_hz"],dtype=np.float64)
        fast_log=np.asarray(model.log10OmegaGW,dtype=np.float64);ref_log=np.asarray(raw["reference_log10_OmegaGW"],dtype=np.float64)
        if len(f)!=len(rf):
            grid_delta=None
        else:
            grid_delta=float(np.max(np.abs(f-rf)))
        errors=np.abs(fast_log-ref_log) if len(fast_log)==len(ref_log) else np.asarray([np.inf])
        refmap=float(row["reference_map_at_fast_root"])
        results.append({
            "source":source_file,"case_id":row["case_id"],"params":row["params"],
            "success":True,"fast_root":float(model.DN_gw[-1]),
            "root_relative_proxy_vs_saved_reference_map":float(abs(float(model.DN_gw[-1])-refmap)/max(abs(refmap),1e-300)),
            "spectrum_max_abs_dex_vs_saved_reference":float(np.max(errors)),
            "spectrum_p95_abs_dex_vs_saved_reference":float(np.quantile(errors,0.95)),
            "frequency_grid_max_abs_shift":grid_delta,"frequency_grid_bitwise_equal":bool(grid_delta==0.0),
            "seconds":seconds,
        })
    payload={
        "status":"candidate_replay_only_reference_maps_reused",
        "candidate_source_sha256":_sha256(ROOT/"stiffgwpy_fast"/"fast_sgwb.py"),
        "reference_source_sha256":_sha256(ROOT/"stiffgwpy_fast"/"reference.py"),
        "results":results,
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"count":len(results),"all_success":all(r["success"] for r in results),"max_root_proxy":max(r.get("root_relative_proxy_vs_saved_reference_map",0) for r in results),"max_spectrum_dex":max(r.get("spectrum_max_abs_dex_vs_saved_reference",0) for r in results),"non_bitwise_grids":sum(not r.get("frequency_grid_bitwise_equal",False) for r in results),"results":results},indent=2))


if __name__ == "__main__":
    main()
