"""Compare sample2 Fast and Prüfer spectra at the same adaptive grid and root."""

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
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

import numpy as np  # noqa: E402

from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

SWEEP = ROOT / "docs" / "sample2_fast_adaptive_sweep_20261105.json"
REFERENCE = ROOT / "docs" / "sample2_reference_root_res8_tail5_gamma0_20261105.json"
OUTPUT = ROOT / "docs" / "sample2_fast_adaptive_propagation_20261105.json"
TARGET = 1e-7
RTOL = 1e-9
Z_TAIL = 5.0


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    sweep_bytes = SWEEP.read_bytes()
    sweep = json.loads(sweep_bytes)
    source_row = next(row for row in sweep["runs"] if row.get("target_dex") == TARGET)
    ref_root_data = json.loads(REFERENCE.read_text(encoding="utf-8"))
    params = sweep["parameters"]
    saved_freqs = np.asarray(source_row["frequency_log10_hz"], dtype=np.float64)
    saved_fast_logs = np.asarray(source_row["log10_OmegaGW"], dtype=np.float64)

    model = LCDM_SG(**params)
    started_fast = time.perf_counter()
    fast_result = model.SGWB_iter(
        engine="fast",
        fallback=False,
        tol=1e-6,
        z_tail=Z_TAIL,
        freq_res=1.0,
        threads=1,
        accuracy_mode="fast",
        physical_guard=False,
        freq_grid="adaptive",
        freq_grid_target=TARGET,
        freq_grid_max_points=1500,
    )
    fast_runtime = time.perf_counter() - started_fast
    if fast_result is not model:
        raise RuntimeError(f"Fast replay failed: {model.fast_failure_reason}")
    fast_freqs = np.asarray(model.f, dtype=np.float64)
    fast_logs = np.asarray(model.log10OmegaGW, dtype=np.float64)
    if not np.array_equal(fast_freqs, saved_freqs) or not np.array_equal(fast_logs, saved_fast_logs):
        raise RuntimeError("Fast replay does not reproduce the saved target-1e-7 arrays")

    dn_eff = float(model.cosmo_param["DN_eff"])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ref_started = time.perf_counter()
        ref_map = R.run_reference(
            LCDM_SG(**params),
            dn_eff=dn_eff,
            freq_subset=fast_freqs,
            self_consistent=False,
            z_tail=Z_TAIL,
            rtol=RTOL,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
        ref_runtime = time.perf_counter() - ref_started
    warnings_seen = [str(item.message) for item in caught]
    ref_freqs = np.asarray(ref_map["freqs"], dtype=np.float64)
    ref_logs = np.asarray(ref_map["log10OmegaGW"], dtype=np.float64)
    if not np.array_equal(ref_freqs, fast_freqs):
        raise RuntimeError("Reference did not preserve Fast's native frequency nodes")

    abs_dex = np.abs(fast_logs - ref_logs)
    worst = int(np.argmax(abs_dex))
    fast_tail = np.asarray(model.handoff_eps, dtype=np.float64) >= 0.0
    handoff_eps = np.asarray(model.handoff_eps, dtype=np.float64)
    ref_tail = np.asarray(ref_map["used_tail"], dtype=bool)
    same_tail = fast_tail == ref_tail
    tail_error_summary = {}
    for used_tail in (False, True):
        mask = fast_tail == used_tail
        tail_error_summary[str(used_tail).lower()] = {
            "mode_count": int(np.count_nonzero(mask)),
            "max_abs_dex": float(np.max(abs_dex[mask])) if np.any(mask) else None,
            "p95_abs_dex": float(np.percentile(abs_dex[mask], 95)) if np.any(mask) else None,
            "handoff_eps_max": float(np.max(handoff_eps[mask])) if np.any(mask & (handoff_eps >= 0.0)) else None,
        }
    local_quad_error = getattr(model, "estimated_DN_quadrature_error_rel", None)
    if local_quad_error is not None and not np.isfinite(float(local_quad_error)):
        local_quad_error = None
    record = {
        "schema_version": 1,
        "experiment": "sample2 Fast adaptive target1e-7 same-frequency fixed-background propagation check",
        "status": "complete",
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        "parameters": params,
        "target_dex": TARGET,
        "reference_root_800_point": float(ref_root_data["params"]["DN_eff"]),
        "fast_root_on_target_grid": dn_eff,
        "reference_DN_gw_map_at_fast_root": float(ref_map["DN_gw"]),
        "same_grid_relative_DN_gw_map_error": abs(float(model.DN_gw[-1]) - float(ref_map["DN_gw"])) / max(abs(float(ref_map["DN_gw"])), np.finfo(float).tiny),
        "fast_root_to_800_point_reference_root_relative_delta": abs(dn_eff - float(ref_root_data["params"]["DN_eff"])) / max(abs(float(ref_root_data["params"]["DN_eff"])), np.finfo(float).tiny),
        "fast_spectrum_max_abs_dex_error_same_nodes": float(np.max(abs_dex)),
        "fast_spectrum_p95_abs_dex_error_same_nodes": float(np.percentile(abs_dex, 95)),
        "worst_frequency_log10_hz": float(fast_freqs[worst]),
        "worst_fast_log10_omega": float(fast_logs[worst]),
        "worst_reference_log10_omega": float(ref_logs[worst]),
        "tail_semantics_mismatch_count": int(np.count_nonzero(~same_tail)),
        "fast_tail_count": int(np.count_nonzero(fast_tail)),
        "reference_tail_count": int(np.count_nonzero(ref_tail)),
        "spectrum_error_by_fast_tail_status": tail_error_summary,
        "worst_mode_index": worst,
        "worst_mode_used_tail": bool(fast_tail[worst]),
        "worst_mode_handoff_eps": float(handoff_eps[worst]),
        "fast_freq_grid_error": float(model.freq_grid_error),
        "fast_local_quadrature_error_rel": None if local_quad_error is None else float(local_quad_error),
        "fast_runtime_replay_s": fast_runtime,
        "reference_fixed_map_runtime_s": ref_runtime,
        "warnings": warnings_seen,
        "source_hashes": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "sweep_artifact": hashlib.sha256(sweep_bytes).hexdigest(),
            "runner.py": _sha(Path(__file__).resolve()),
        },
        "frequencies_log10_hz": fast_freqs.tolist(),
        "fast_log10_omega": fast_logs.tolist(),
        "reference_log10_omega": ref_logs.tolist(),
        "fast_handoff_eps": [float(x) if np.isfinite(x) else None for x in handoff_eps],
        "limitation": "Fixed-background map and spectrum comparison on one Fast adaptive grid. The independent same-grid Reference root is not solved here; this separates grid refinement from propagation error only at the Fast root.",
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "same_grid_DN_rel_error": record["same_grid_relative_DN_gw_map_error"],
        "max_spectrum_dex_error": record["fast_spectrum_max_abs_dex_error_same_nodes"],
        "fast_root_to_800_root_delta": record["fast_root_to_800_point_reference_root_relative_delta"],
        "tail_mismatches": record["tail_semantics_mismatch_count"],
        "worst_frequency": record["worst_frequency_log10_hz"],
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
