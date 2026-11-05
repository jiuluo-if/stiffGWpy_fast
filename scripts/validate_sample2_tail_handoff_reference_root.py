"""Compare the z_tail=7 Fast handoff with a same-grid self-consistent Prüfer root."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
import traceback
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

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

SWEEP = ROOT / "docs" / "sample2_tail_handoff_sweep_20261105.json"
REFERENCE_ROOT = ROOT / "docs" / "sample2_reference_root_res8_tail5_gamma0_20261105.json"
OUTPUT = ROOT / "docs" / "sample2_tail_handoff_reference_root_20261105.json"
TARGET = 1e-7
Z_TAIL = 7.0
RTOL = 1e-9
TOL = 1e-9


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def _write(record: dict) -> None:
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def main() -> int:
    sweep_bytes = SWEEP.read_bytes()
    sweep = json.loads(sweep_bytes)
    fast_row = next(row for row in sweep["runs"] if row["z_tail"] == Z_TAIL)
    reference_root_data = json.loads(REFERENCE_ROOT.read_text(encoding="utf-8"))
    freqs = np.asarray(fast_row["frequency_log10_hz"], dtype=np.float64)
    fast_logs = np.asarray(fast_row["log10_OmegaGW"], dtype=np.float64)
    params = sweep["parameters"]
    model = LCDM_SG(**params)
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    map_history: list[dict] = []
    record = {
        "schema_version": 1,
        "experiment": "sample2 adaptive target1e-7 z_tail7 same-grid Fast/Prüfer root check",
        "status": "running",
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": params,
        "target_dex": TARGET,
        "z_tail": Z_TAIL,
        "fast_root": float(fast_row["DN_eff_final"]),
        "fast_handoff_eps_max": float(fast_row["handoff_eps_max"]),
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "full_reference_root_z_tail5": float(reference_root_data["params"]["DN_eff"]),
        "map_history": map_history,
        "runner_sha256": _sha(Path(__file__).resolve()),
    }
    _write(record)
    original = R._spectrum_at_dn

    def capture_map(m, freq_values, dn_eff, z_tail, rtol, workers=None,
                    backend="cartesian", z_start=-12.0):
        started = time.perf_counter()
        result = original(m, freq_values, dn_eff, z_tail, rtol,
                          workers=workers, backend=backend, z_start=z_start)
        g2, qerr, ierr = R.integrate_spectrum(freq_values, result[0], result[1])
        mapped = float(gp.Neff0 * g2 / omega_nu)
        map_history.append({
            "DN_eff_input": float(dn_eff),
            "DN_gw_map": mapped,
            "fixed_point_residual": float(dn_eff - params["DN_eff"] - mapped),
            "quadrature_error": float(qerr),
            "interpolation_error": float(ierr),
            "runtime_s": time.perf_counter() - started,
        })
        record["status"] = "reference_root_iteration_in_progress"
        _write(record)
        return result

    R._spectrum_at_dn = capture_map
    warnings_seen: list[str] = []
    started = time.perf_counter()
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = R.run_reference(
                model,
                dn_eff=None,
                freq_subset=freqs,
                self_consistent=True,
                tol=TOL,
                max_iter=60,
                z_tail=Z_TAIL,
                rtol=RTOL,
                workers=1,
                physical_guard=False,
                backend="prufer",
            )
            warnings_seen = [str(item.message) for item in caught]
    except Exception as exc:
        record["status"] = "failed"
        record["error"] = f"{type(exc).__name__}: {exc}"
        record["traceback"] = traceback.format_exc()
        record["warnings"] = warnings_seen
        record["total_runtime_s"] = time.perf_counter() - started
        _write(record)
        raise
    finally:
        R._spectrum_at_dn = original

    ref_logs = np.asarray(result["log10OmegaGW"], dtype=np.float64)
    abs_dex = np.abs(fast_logs - ref_logs)
    fast_root = float(fast_row["DN_eff_final"])
    ref_root = float(result["DN_eff"])
    record.update({
        "status": "complete",
        "reference_root_on_same_grid": ref_root,
        "reference_DN_gw_at_root": float(result["DN_gw"]),
        "reference_fixed_point_residual_at_root": float(ref_root - params["DN_eff"] - result["DN_gw"]),
        "fast_to_same_grid_reference_root_relative_error": abs(fast_root - ref_root) / max(abs(ref_root), np.finfo(float).tiny),
        "same_grid_reference_root_to_z_tail5_800_root_relative_delta": abs(ref_root - record["full_reference_root_z_tail5"]) / max(abs(record["full_reference_root_z_tail5"]), np.finfo(float).tiny),
        "max_abs_dex_spectrum_error_same_grid": float(np.max(abs_dex)),
        "p95_abs_dex_spectrum_error_same_grid": float(np.percentile(abs_dex, 95)),
        "reference_iterations": int(result["n_iter"]),
        "final_quadrature_error": float(result["quadrature_error"]),
        "final_interpolation_error": float(result["interpolation_error"]),
        "warnings": warnings_seen,
        "total_runtime_s": time.perf_counter() - started,
        "fast_root_runtime_s": float(fast_row["runtime_s"]),
        "resources": telemetry(workers=1, threads=1),
        "source_hashes": {
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "sweep_artifact": hashlib.sha256(sweep_bytes).hexdigest(),
            "runner.py": _sha(Path(__file__).resolve()),
        },
        "frequencies_log10_hz": freqs.tolist(),
        "fast_log10_omega": fast_logs.tolist(),
        "reference_log10_omega": ref_logs.tolist(),
        "limitation": "Same-grid adaptive-root check at z_tail=7. The comparison with the 800-node z_tail=5 root changes both frequency sampling and tail depth; it is not a continuous-domain certificate.",
    })
    _write(record)
    print(json.dumps({
        "status": record["status"],
        "fast_root": fast_root,
        "reference_root": ref_root,
        "fast_to_same_grid_root_rel": record["fast_to_same_grid_reference_root_relative_error"],
        "root_to_800_z5_rel": record["same_grid_reference_root_to_z_tail5_800_root_relative_delta"],
        "max_spectrum_error_dex": record["max_abs_dex_spectrum_error_same_grid"],
        "runtime_s": record["total_runtime_s"],
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
