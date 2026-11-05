"""Check the z_tail=7 Fast point against the frozen default Prüfer reference."""

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

SWEEP = ROOT / "docs" / "sample2_tail_handoff_sweep_20261105.json"
ADAPTIVE_SWEEP = ROOT / "docs" / "sample2_fast_adaptive_sweep_20261105.json"
REFERENCE_ROOTS = ROOT / "docs" / "sample2_fast_adaptive_reference_20261105.json"
OUTPUT = ROOT / "docs" / "sample2_tail_handoff_frozen_reference_20261105.json"
TARGET = 1e-7
FAST_Z_TAIL = 7.0
REFERENCE_Z_TAIL = 5.0  # Frozen default in reference.run_reference.
RTOL = 1e-9


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    sweep_bytes = SWEEP.read_bytes()
    sweep = json.loads(sweep_bytes)
    fast_row = next(row for row in sweep["runs"] if row["z_tail"] == FAST_Z_TAIL)
    roots_bytes = REFERENCE_ROOTS.read_bytes()
    roots = json.loads(roots_bytes)
    adaptive_sweep_bytes = ADAPTIVE_SWEEP.read_bytes()
    adaptive_sweep = json.loads(adaptive_sweep_bytes)
    ref_row = next(row for row in roots["targets"] if row["target_dex"] == TARGET)
    freqs = np.asarray(fast_row["frequency_log10_hz"], dtype=np.float64)
    fast_logs = np.asarray(fast_row["log10_OmegaGW"], dtype=np.float64)
    if fast_row["frequency_count"] != ref_row["frequency_count"]:
        raise RuntimeError("Fast and frozen-reference root artifacts have different node counts")
    if fast_row["frequency_sha256"] != ref_row["frequency_sha256"]:
        raise RuntimeError("frequency hashes do not match")
    params = sweep["parameters"]
    if roots["sweep_sha256"] != hashlib.sha256(adaptive_sweep_bytes).hexdigest():
        raise RuntimeError("same-grid Reference artifact is not bound to this Fast sweep")
    if params != adaptive_sweep["parameters"]:
        raise RuntimeError("Fast tail sweep and Reference root use different cosmologies")

    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = R.run_reference(
            LCDM_SG(**params),
            dn_eff=float(fast_row["DN_eff_final"]),
            freq_subset=freqs,
            self_consistent=False,
            # Keep Reference at its unmodified, frozen default tail rule.
            z_tail=REFERENCE_Z_TAIL,
            rtol=RTOL,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
    runtime = time.perf_counter() - started
    ref_freqs = np.asarray(result["freqs"], dtype=np.float64)
    ref_logs = np.asarray(result["log10OmegaGW"], dtype=np.float64)
    if not np.array_equal(ref_freqs, freqs):
        raise RuntimeError("Reference did not preserve the Fast frequency nodes")
    abs_dex = np.abs(fast_logs - ref_logs)
    fast_handoff_eps = np.asarray(fast_row["handoff_eps"], dtype=np.float64)
    tail_mask = fast_handoff_eps >= 0.0
    tail_summary = {}
    for label, mask in (("tail", tail_mask), ("no_tail", ~tail_mask)):
        values = abs_dex[mask]
        tail_summary[label] = {
            "mode_count": int(np.count_nonzero(mask)),
            "max_abs_dex": float(np.max(values)) if values.size else None,
            "p95_abs_dex": float(np.percentile(values, 95)) if values.size else None,
        }
    ref_root = float(ref_row["reference_root_on_same_grid"])
    fast_root = float(fast_row["DN_eff_final"])
    ref_map = float(result["DN_gw"])
    record = {
        "schema_version": 1,
        "experiment": "sample2 z_tail7 Fast versus frozen-default Prüfer at identical frequency nodes",
        "status": "complete",
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        "parameters": params,
        "target_dex": TARGET,
        "frequency_count": int(freqs.size),
        "frequency_sha256": fast_row["frequency_sha256"],
        "fast_z_tail": FAST_Z_TAIL,
        "reference_z_tail": REFERENCE_Z_TAIL,
        "fast_root": fast_root,
        "frozen_reference_same_grid_root": ref_root,
        "fast_to_frozen_reference_root_relative_error": abs(fast_root - ref_root) / max(abs(ref_root), np.finfo(float).tiny),
        "reference_DN_gw_map_at_fast_root": ref_map,
        "fixed_point_residual_at_fast_root": fast_root - float(params["DN_eff"]) - ref_map,
        "fast_to_frozen_reference_DN_gw_map_relative_error": abs(float(fast_row["DN_gw_final"]) - ref_map) / max(abs(ref_map), np.finfo(float).tiny),
        "max_abs_spectrum_error_dex_at_same_DN_eff": float(np.max(abs_dex)),
        "p95_abs_spectrum_error_dex_at_same_DN_eff": float(np.percentile(abs_dex, 95)),
        "spectrum_error_by_fast_tail_status": tail_summary,
        "frequencies_log10_hz": freqs.tolist(),
        "fast_log10_omega": fast_logs.tolist(),
        "reference_log10_omega": ref_logs.tolist(),
        "fast_handoff_eps": fast_handoff_eps.tolist(),
        "reference_runtime_s": runtime,
        "warnings": [str(item.message) for item in caught],
        "resources": telemetry(workers=1, threads=1),
        "source_hashes": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "fast_tail_sweep": hashlib.sha256(sweep_bytes).hexdigest(),
            "fast_adaptive_sweep": hashlib.sha256(adaptive_sweep_bytes).hexdigest(),
            "same_grid_reference_roots": hashlib.sha256(roots_bytes).hexdigest(),
            "runner.py": _sha(Path(__file__).resolve()),
        },
        "limitation": "One cosmology and one adaptive frequency set; the same-grid Reference root is from the frozen z_tail=5 run. This is not a regional certificate or a paired performance benchmark.",
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "fast_root": fast_root,
        "reference_root": ref_root,
        "root_relative_error": record["fast_to_frozen_reference_root_relative_error"],
        "map_relative_error": record["fast_to_frozen_reference_DN_gw_map_relative_error"],
        "max_spectrum_error_dex": record["max_abs_spectrum_error_dex_at_same_DN_eff"],
        "reference_runtime_s": runtime,
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
