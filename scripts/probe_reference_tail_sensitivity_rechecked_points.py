"""Measure frozen-tail sensitivity on the rechecked difficult-point Fast grids."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402

from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    "high_r_red": (
        "docs/fast_cr1_red_extreme_fullgrid_prufer_20261105.json",
        "docs/fast_cr1_red_extreme_reference_root_defaulttail_20261105.json",
    ),
    "ninf_222": (
        "docs/fast_ninf222_fullgrid_prufer_20261105.json",
        "docs/fast_ninf222_reference_root_defaulttail_20261105.json",
    ),
    "positive_cr_low_temperature": (
        "docs/fast_cr1_lowT_fullgrid_prufer_20261105.json",
        "docs/fast_cr1_lowT_reference_root_defaulttail_20261105.json",
    ),
}
OUTPUT = ROOT / "docs" / "reference_tail_sensitivity_rechecked_20261105.json"
DEFAULTTAIL_SPECTRA = ROOT / "docs" / "defaulttail_rechecked_spectra_20261105.json"
TAIL_DEPTHS = (8.0, 10.0)
RTOL = 1e-9


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    rows: dict[str, dict] = {}
    warnings_seen: list[str] = []
    hashes: dict[str, dict] = {}
    defaulttail_bytes = DEFAULTTAIL_SPECTRA.read_bytes()
    defaulttail = json.loads(defaulttail_bytes)
    started_all = time.perf_counter()
    for name, (parent_rel, root_rel) in CASES.items():
        parent_path, root_path = ROOT / parent_rel, ROOT / root_rel
        parent_bytes, root_bytes = parent_path.read_bytes(), root_path.read_bytes()
        parent, root = json.loads(parent_bytes), json.loads(root_bytes)
        fast = parent["fast"]
        freqs = np.asarray(fast["frequencies_log10_hz"], dtype=np.float64)
        freq_hash = hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest()
        if freq_hash != root["frequency_sha256"]:
            raise RuntimeError(f"{name}: frequency hash mismatch")
        params = parent.get("parameters", parent.get("parameters_input"))
        if params is None:
            raise RuntimeError(f"{name}: missing input parameters")
        model = LCDM_SG(**params)
        fast_root = float(fast["DN_eff_final"])
        fast_logs = np.log10(np.maximum(np.asarray(fast["omega_gw"], dtype=np.float64), np.finfo(float).tiny))
        z5_fast_root_logs = np.asarray(defaulttail["cases"][name]["same_fast_root_reference_log10_omega"], dtype=np.float64)
        if z5_fast_root_logs.shape != fast_logs.shape:
            raise RuntimeError(f"{name}: frozen z5 spectrum and Fast arrays have different lengths")
        by_tail: dict[str, dict] = {}
        for z_tail in TAIL_DEPTHS:
            started = time.perf_counter()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = R.run_reference(
                    model,
                    dn_eff=fast_root,
                    freq_subset=freqs,
                    self_consistent=False,
                    z_tail=z_tail,
                    rtol=RTOL,
                    workers=1,
                    physical_guard=False,
                    backend="prufer",
                )
            warnings_seen.extend(f"{name}, z_tail={z_tail}: {item.message}" for item in caught)
            if not np.array_equal(freqs, np.asarray(result["freqs"], dtype=np.float64)):
                raise RuntimeError(f"{name}, z_tail={z_tail}: frequency nodes changed")
            ref_integrand = np.asarray(result["Ogw"], dtype=np.float64) - np.asarray(result["Oj"], dtype=np.float64)
            positive = (ref_integrand > 0.0) & (np.asarray(fast["omega_gw"]) > 0.0)
            fast_root_map = float(fast["DN_gw"])
            map_value = float(result["DN_gw"])
            abs_dex = np.abs(fast_logs[positive] - np.log10(ref_integrand[positive]))
            ref_log_values = np.log10(np.maximum(ref_integrand, np.finfo(float).tiny))
            by_tail[str(int(z_tail))] = {
                "reference_DN_gw_at_fast_root": map_value,
                "map_relative_difference_from_fast": abs(fast_root_map - map_value) / max(abs(map_value), np.finfo(float).tiny),
                "reference_fixed_point_residual_at_fast_root": fast_root - float(params["DN_eff"]) - map_value,
                "reference_used_tail_count": int(np.count_nonzero(result["used_tail"])),
                "positive_spectrum_nodes": int(np.count_nonzero(positive)),
                "fast_spectrum_max_abs_error_dex": float(np.max(abs_dex)) if abs_dex.size else None,
                "fast_spectrum_p95_abs_error_dex": float(np.percentile(abs_dex, 95)) if abs_dex.size else None,
                "spectrum_change_vs_frozen_z5_max_abs_dex": float(np.max(np.abs(ref_log_values - z5_fast_root_logs))),
                "reference_log10_omega": [float(x) if np.isfinite(x) else None for x in ref_log_values],
                "reference_used_tail": np.asarray(result["used_tail"], dtype=bool).tolist(),
                "runtime_s": time.perf_counter() - started,
            }
        rows[name] = {
            "parameters": params,
            "frequency_count": int(freqs.size),
            "frequency_sha256": freq_hash,
            "fast_root": fast_root,
            "frozen_z5_reference_root": float(root["reference_root_DN_eff"]),
            "fast_to_frozen_z5_root_relative_error": float(root["fast_to_reference_root_relative_error"]),
            "reference_maps": by_tail,
        }
        hashes[name] = {
            "fast_parent_sha256": hashlib.sha256(parent_bytes).hexdigest(),
            "defaulttail_root_sha256": hashlib.sha256(root_bytes).hexdigest(),
            "defaulttail_spectrum_sha256": hashlib.sha256(defaulttail_bytes).hexdigest(),
        }

    record = {
        "schema_version": 1,
        "experiment": "Reference frozen-tail sensitivity at difficult-point Fast roots",
        "status": "complete",
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        "reference": {"backend": "prufer", "z_tail_diagnostic_values": list(TAIL_DEPTHS), "frozen_default_z_tail": 5.0, "rtol": RTOL},
        "cases": rows,
        "warnings": warnings_seen,
        "total_runtime_s": time.perf_counter() - started_all,
        "resources": telemetry(workers=1, threads=1),
        "source_hashes": {"reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"), "runner.py": _sha(Path(__file__).resolve())},
        "input_hashes": hashes,
        "limitations": [
            "z_tail=8 and 10 are read-only sensitivity probes and do not replace frozen default z_tail=5 acceptance records.",
            "Only three named points and their existing Fast frequency grids are covered; no continuous-region claim is made.",
            "Reference fixed-point maps are evaluated at the Fast roots, not iterated at each diagnostic tail depth.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": record["status"],
        "cases": {name: {depth: {key: vals[key] for key in ("map_relative_difference_from_fast", "fast_spectrum_max_abs_error_dex", "fast_spectrum_p95_abs_error_dex")} for depth, vals in row["reference_maps"].items()} for name, row in rows.items()},
        "runtime_s": record["total_runtime_s"],
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
