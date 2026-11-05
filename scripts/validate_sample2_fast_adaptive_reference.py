"""Validate built-in Fast adaptive grids with self-consistent same-grid Prüfer roots."""

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

SWEEP = ROOT / "docs" / "sample2_fast_adaptive_sweep_20261105.json"
REFERENCE_ROOT = ROOT / "docs" / "sample2_reference_root_res8_tail5_gamma0_20261105.json"
OUTPUT = ROOT / "docs" / "sample2_fast_adaptive_reference_20261105.json"
TARGETS = (1e-5, 1e-6, 1e-7)
RTOL = 1e-9
TOL = 1e-9
Z_TAIL = 5.0


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _write(record: dict) -> None:
    OUTPUT.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    sweep_bytes = SWEEP.read_bytes()
    sweep = json.loads(sweep_bytes)
    reference_root = json.loads(REFERENCE_ROOT.read_text(encoding="utf-8"))
    full_root = float(reference_root["params"]["DN_eff"])
    rows = {float(row["target_dex"]): row for row in sweep["runs"] if row.get("target_dex") is not None}
    record = {
        "schema_version": 1,
        "experiment": "same-frequency self-consistent Prüfer checks for sample2 Fast adaptive grids",
        "status": "running",
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "sweep_artifact": SWEEP.name,
        "sweep_sha256": hashlib.sha256(sweep_bytes).hexdigest(),
        "full_reference_root_artifact": REFERENCE_ROOT.name,
        "full_reference_root_sha256": _sha(REFERENCE_ROOT),
        "full_reference_root_DN_eff": full_root,
        "reference": {"backend": "prufer", "rtol": RTOL, "z_tail": Z_TAIL, "root_tolerance": TOL, "workers": 1},
        "targets": [],
        "runner_sha256": _sha(Path(__file__).resolve()),
    }
    _write(record)
    original = R._spectrum_at_dn
    model = LCDM_SG(**sweep["parameters"])
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2

    for target in TARGETS:
        row = rows[float(target)]
        freqs = np.asarray(row["frequency_log10_hz"], dtype=np.float64)
        fast_logs = np.asarray(row["log10_OmegaGW"], dtype=np.float64)
        if freqs.size != fast_logs.size or not np.isfinite(freqs).all() or not np.isfinite(fast_logs).all():
            raise RuntimeError(f"invalid Fast arrays at target {target}")
        history: list[dict] = []
        started = time.perf_counter()

        def capture_map(m, freq_values, dn_eff, z_tail, rtol, workers=None,
                        backend="cartesian", z_start=-12.0):
            map_started = time.perf_counter()
            result = original(
                m, freq_values, dn_eff, z_tail, rtol,
                workers=workers, backend=backend, z_start=z_start,
            )
            g2, qerr, ierr = R.integrate_spectrum(freq_values, result[0], result[1])
            mapped = float(gp.Neff0 * g2 / omega_nu)
            history.append({
                "DN_eff_input": float(dn_eff),
                "DN_gw_map": mapped,
                "fixed_point_residual": float(dn_eff - sweep["parameters"]["DN_eff"] - mapped),
                "quadrature_error": float(qerr),
                "interpolation_error": float(ierr),
                "runtime_s": time.perf_counter() - map_started,
            })
            record["status"] = "reference_runs_in_progress"
            _write(record)
            return result

        R._spectrum_at_dn = capture_map
        warnings_seen: list[str] = []
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                ref_result = R.run_reference(
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
            record["failed_target_dex"] = target
            record["error"] = f"{type(exc).__name__}: {exc}"
            record["traceback"] = traceback.format_exc()
            record["target_map_history"] = history
            _write(record)
            raise
        finally:
            R._spectrum_at_dn = original

        ref_logs = np.asarray(ref_result["log10OmegaGW"], dtype=np.float64)
        abs_dex = np.abs(fast_logs - ref_logs)
        ref_root = float(ref_result["DN_eff"])
        target_record = {
            "target_dex": float(target),
            "frequency_count": int(freqs.size),
            "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
            "fast_root": float(row["DN_eff_final"]),
            "reference_root_on_same_grid": ref_root,
            "reference_iterations": int(ref_result["n_iter"]),
            "fast_to_same_grid_reference_root_relative_error": abs(float(row["DN_eff_final"]) - ref_root) / max(abs(ref_root), np.finfo(float).tiny),
            "same_grid_reference_root_to_800_point_root_relative_error": abs(ref_root - full_root) / max(abs(full_root), np.finfo(float).tiny),
            "max_abs_dex_spectrum_error_at_native_nodes": float(np.max(abs_dex)),
            "p95_abs_dex_spectrum_error_at_native_nodes": float(np.percentile(abs_dex, 95)),
            "fast_freq_grid_error": float(row["freq_grid_error"]),
            "reference_fixed_point_residual_at_root": float(ref_result["DN_eff"] - sweep["parameters"]["DN_eff"] - ref_result["DN_gw"]),
            "reference_quadrature_error": float(ref_result["quadrature_error"]),
            "reference_interpolation_error": float(ref_result["interpolation_error"]),
            "reference_runtime_s": time.perf_counter() - started,
            "reference_map_history": history,
            "warnings": warnings_seen,
        }
        record["targets"].append(target_record)
        record["status"] = "running"
        _write(record)

    record.update({
        "status": "complete",
        "source_hashes": {
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "runner.py": _sha(Path(__file__).resolve()),
        },
        "limitations": [
            "Each Reference root uses the corresponding Fast frequency subset, so it is an exact-node cross-check but still inherits that grid's discretization error.",
            "The target sweep is not a same-grid paired P50/P95 performance study; Fast timings are from the prior sequential feasibility sweep.",
            "This tests the existing adaptive path and does not change production defaults or add a dispatch rule.",
        ],
        "resources": telemetry(workers=1, threads=1),
    })
    _write(record)
    print(json.dumps({
        "status": record["status"],
        "targets": [{k: row[k] for k in (
            "target_dex", "frequency_count", "fast_to_same_grid_reference_root_relative_error",
            "same_grid_reference_root_to_800_point_root_relative_error",
            "max_abs_dex_spectrum_error_at_native_nodes", "reference_runtime_s",
        )} for row in record["targets"]],
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
