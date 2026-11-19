"""Direct z8 closure checks for representative CR1 transition-patch points."""

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
OUTPUT = HERE / "validate_exact_quarter_reheat_kappa_patch.json"
CASE_IDS = ("transition_00", "transition_02", "transition_06",
            "transition_08", "transition_holdout")


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    screen = json.loads(
        (HERE / "screen_exact_quarter_reheat_kappa_patch.json").read_text(encoding="utf-8")
    )
    pairs = {row["input"]["case_id"]: row for row in screen["paired_rows"]}
    results = []
    if OUTPUT.exists():
        previous = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if previous.get("source_sha256", {}).get("reference.py") == _sha(
            ROOT / "stiffgwpy_fast" / "reference.py"
        ):
            results.extend(previous.get("results", []))
    completed_ids = {row["case_id"] for row in results}
    for case_id in CASE_IDS:
        if case_id in completed_ids:
            continue
        pair = pairs[case_id]
        params = pair["input"]["params"]
        model = LCDM_SG(**params)
        if SGWB_iter_fast(
            model, config=resolve_config("fast", threads=2), tol=1e-6,
            physical_guard=True,
        ) is not model:
            raise RuntimeError(f"Fast failed at {case_id}: {model.fast_failure_reason}")
        fast_root = float(model.DN_gw[-1])
        freqs = np.asarray(model.f, dtype=np.float64)
        fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)

        start = time.perf_counter()
        first = run_reference(
            model, dn_eff=fast_root, freq_subset=freqs, z_tail=8.0,
            rtol=1e-9, self_consistent=False, workers=2,
        )
        first_seconds = time.perf_counter() - start
        ref_freqs = np.asarray(first["freqs"], dtype=np.float64)
        if not np.array_equal(freqs, ref_freqs):
            raise RuntimeError(f"Reference frequency grid differs at {case_id}")
        first_map = float(first["DN_gw"])
        ref_log = np.asarray(first["log10OmegaGW"], dtype=np.float64)
        ref_ogw = np.asarray(first["Ogw"], dtype=np.float64)
        ref_oj = np.asarray(first["Oj"], dtype=np.float64)
        valid = (ref_ogw - ref_oj > 0.0) & np.isfinite(fast_log) & np.isfinite(ref_log)
        spectrum_error = np.abs(fast_log[valid] - ref_log[valid])

        start = time.perf_counter()
        second = run_reference(
            model, dn_eff=first_map, freq_subset=freqs, z_tail=8.0,
            rtol=1e-9, self_consistent=False, workers=2,
        )
        second_seconds = time.perf_counter() - start
        second_map = float(second["DN_gw"])
        results.append({
            "case_id": case_id, "params": params,
            "old_fast_failure": pair["baseline_fast"]["failure_reason"],
            "candidate_fast_success": pair["candidate_fast"]["success"],
            "fast_root": fast_root, "reference_map_at_fast_root": first_map,
            "reference_map_at_first_map": second_map,
            "root_relative_to_closed_reference": abs(fast_root-second_map)/max(abs(second_map),1e-300),
            "reference_closure_residual_relative": abs(second_map-first_map)/max(abs(second_map),1e-300),
            "spectrum_max_abs_dex": float(np.max(spectrum_error)),
            "spectrum_p95_abs_dex": float(np.quantile(spectrum_error,0.95)),
            "valid_spectrum_nodes": int(valid.sum()), "frequency_grid_bitwise_equal": True,
            "first_reference_seconds": first_seconds,
            "closure_reference_seconds": second_seconds,
            "freqs_log10_hz": freqs.tolist(),
            "fast_log10_OmegaGW": fast_log.tolist(),
            "reference_log10_OmegaGW": ref_log.tolist(),
            "reference_Ogw": ref_ogw.tolist(), "reference_Oj": ref_oj.tolist(),
        })
        partial = {
            "status": "in_progress_transition_patch_reference_checks",
            "reference_code_modified": False, "reference_tail_z": 8.0,
            "source_sha256": {
                name: _sha(ROOT / "stiffgwpy_fast" / name)
                for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
            },
            "results": results,
        }
        OUTPUT.write_text(json.dumps(partial, indent=2)+"\n", encoding="utf-8")
    payload = {
        "status": "complete_reheat_kappa_patch_reference_checks",
        "reference_code_modified": False, "reference_tail_z": 8.0,
        "source_sha256": {
            name: _sha(ROOT / "stiffgwpy_fast" / name)
            for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
        },
        "results": results,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "results": [{key: value for key, value in row.items()
                     if key not in {"freqs_log10_hz", "fast_log10_OmegaGW",
                                    "reference_log10_OmegaGW", "reference_Ogw", "reference_Oj"}}
                    for row in results],
    }, indent=2))


if __name__ == "__main__":
    main()
