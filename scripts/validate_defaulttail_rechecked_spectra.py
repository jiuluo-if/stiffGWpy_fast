"""Compare rechecked Fast spectra to frozen-default Prüfer maps at their roots."""

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
        "parameters",
    ),
    "ninf_222": (
        "docs/fast_ninf222_fullgrid_prufer_20261105.json",
        "docs/fast_ninf222_reference_root_defaulttail_20261105.json",
        "parameters",
    ),
    "positive_cr_low_temperature": (
        "docs/fast_cr1_lowT_fullgrid_prufer_20261105.json",
        "docs/fast_cr1_lowT_reference_root_defaulttail_20261105.json",
        "parameters_input",
    ),
}
OUTPUT = ROOT / "docs" / "defaulttail_rechecked_spectra_20261105.json"
Z_TAIL = 5.0
RTOL = 1e-9


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    rows = {}
    warnings_seen: list[str] = []
    started_all = time.perf_counter()
    input_hashes = {}
    for name, (parent_rel, root_rel, params_key) in CASES.items():
        parent_path = ROOT / parent_rel
        root_path = ROOT / root_rel
        parent_bytes, root_bytes = parent_path.read_bytes(), root_path.read_bytes()
        parent, root = json.loads(parent_bytes), json.loads(root_bytes)
        fast = parent["fast"]
        freqs = np.asarray(fast["frequencies_log10_hz"], dtype=np.float64)
        expected_freq_hash = hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest()
        if expected_freq_hash != root["frequency_sha256"]:
            raise RuntimeError(f"{name}: root and Fast frequency grids differ")
        if root["reference"]["z_tail"] != Z_TAIL:
            raise RuntimeError(f"{name}: root record is not on frozen default z_tail=5")
        params = parent[params_key]
        model = LCDM_SG(**params)
        ref_root = float(root["reference_root_DN_eff"])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = R.run_reference(
                model,
                dn_eff=ref_root,
                freq_subset=freqs,
                self_consistent=False,
                z_tail=Z_TAIL,
                rtol=RTOL,
                workers=1,
                physical_guard=False,
                backend="prufer",
            )
        warnings_seen.extend(f"{name}: {item.message}" for item in caught)
        ref_freqs = np.asarray(result["freqs"], dtype=np.float64)
        if not np.array_equal(freqs, ref_freqs):
            raise RuntimeError(f"{name}: Reference changed frequency nodes")
        fast_omega = np.asarray(fast["omega_gw"], dtype=np.float64)
        ref_integrand = np.asarray(result["Ogw"], dtype=np.float64) - np.asarray(result["Oj"], dtype=np.float64)
        ref_logs = np.asarray(result["log10OmegaGW"], dtype=np.float64)
        fast_logs = np.log10(np.maximum(fast_omega, np.finfo(float).tiny))
        valid_reference = ref_integrand > 0.0
        jointly_positive = valid_reference & (fast_omega > 0.0)
        abs_dex = np.abs(fast_logs[jointly_positive] - np.log10(ref_integrand[jointly_positive]))
        abs_linear = np.abs(fast_omega[valid_reference] - ref_integrand[valid_reference])
        relative_linear = abs_linear / np.maximum(np.abs(ref_integrand[valid_reference]), np.finfo(float).tiny)
        fast_root = float(fast["DN_eff_final"])
        rows[name] = {
            "parameters": params,
            "frequency_count": int(freqs.size),
            "frequency_sha256": expected_freq_hash,
            "fast_root": fast_root,
            "frozen_reference_root": ref_root,
            "fast_to_reference_root_relative_error": float(root["fast_to_reference_root_relative_error"]),
            "fixed_point_residual_at_reference_root": float(root["fixed_point_residual_at_root"]),
            "reference_DN_gw_at_root": float(result["DN_gw"]),
            "reference_map_residual_at_root": float(ref_root - params["DN_eff"] - float(result["DN_gw"])),
            "reference_positive_node_count": int(np.count_nonzero(valid_reference)),
            "joint_positive_node_count": int(np.count_nonzero(jointly_positive)),
            "fast_nonpositive_on_reference_positive_count": int(np.count_nonzero(valid_reference & (fast_omega <= 0.0))),
            "reference_nonpositive_or_floor_node_count": int(np.count_nonzero(ref_integrand <= 0.0)),
            "spectrum_max_abs_error_dex_on_positive_nodes": float(np.max(abs_dex)) if abs_dex.size else None,
            "spectrum_p95_abs_error_dex_on_positive_nodes": float(np.percentile(abs_dex, 95)) if abs_dex.size else None,
            "spectrum_max_abs_linear_error_on_positive_nodes": float(np.max(abs_linear)) if abs_linear.size else None,
            "spectrum_max_relative_linear_error_on_positive_nodes": float(np.max(relative_linear)) if relative_linear.size else None,
            "fast_tail_status_available_in_parent_artifact": False,
            "warnings": [str(item.message) for item in caught],
            "fast_log10_omega": fast_logs.tolist(),
            "reference_log10_omega": ref_logs.tolist(),
            "reference_integrand": ref_integrand.tolist(),
            "frequencies_log10_hz": freqs.tolist(),
        }
        input_hashes[name] = {
            "fast_parent_sha256": hashlib.sha256(parent_bytes).hexdigest(),
            "reference_root_sha256": hashlib.sha256(root_bytes).hexdigest(),
        }

    record = {
        "schema_version": 1,
        "experiment": "same-grid Fast versus frozen-default Prüfer spectrum at rechecked roots",
        "status": "complete",
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        "reference": {"backend": "prufer", "z_tail": Z_TAIL, "rtol": RTOL, "self_consistent_root_artifact": True},
        "cases": rows,
        "warnings": warnings_seen,
        "total_runtime_s": time.perf_counter() - started_all,
        "resources": telemetry(workers=1, threads=1),
        "source_hashes": {
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "runner.py": _sha(Path(__file__).resolve()),
        },
        "input_hashes": input_hashes,
        "limitations": [
            "Three named points only; these do not establish continuous-region coverage.",
            "Fast spectra are at the Fast roots and Reference spectra are at the frozen-default Reference roots; root differences are reported separately.",
            "Frequency grids are inherited from prior validated Fast artifacts and are not independently refined here.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": record["status"],
        "cases": {k: {key: v[key] for key in ("fast_to_reference_root_relative_error", "spectrum_max_abs_error_dex_on_positive_nodes", "spectrum_p95_abs_error_dex_on_positive_nodes", "reference_positive_node_count", "fast_nonpositive_on_reference_positive_count", "reference_map_residual_at_root")} for k, v in rows.items()},
        "runtime_s": record["total_runtime_s"],
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
