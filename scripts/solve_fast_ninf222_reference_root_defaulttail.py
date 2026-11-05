"""Recheck the N_inf=222 grid against the frozen default-tail Prüfer Reference."""

from __future__ import annotations

import hashlib
import inspect
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
from stiffgwpy_fast import reference as ref  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARENT = ROOT / "docs" / "fast_ninf222_fullgrid_prufer_20261105.json"
OUTPUT = ROOT / "docs" / "fast_ninf222_reference_root_defaulttail_20261105.json"
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
    parent_bytes = PARENT.read_bytes()
    parent = json.loads(parent_bytes)
    freqs = np.asarray(parent["fast"]["frequencies_log10_hz"], dtype=np.float64)
    fast_x = float(parent["fast"]["DN_eff_final"])
    params = parent["parameters"]
    model = LCDM_SG(**params)
    original = ref._spectrum_at_dn
    map_history: list[dict] = []
    record = {
        "schema_version": 1,
        "experiment": "self-consistent Prüfer Reference root on N_inf=222 Fast grid with frozen default z_tail=5",
        "status": "running",
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": params,
        "parent_artifact": PARENT.name,
        "parent_artifact_sha256": hashlib.sha256(parent_bytes).hexdigest(),
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "reference": {
            "backend": "prufer", "rtol": RTOL, "z_tail": Z_TAIL,
            "root_tolerance": TOL, "workers": 1,
            "self_consistent": True, "physical_guard": False,
        },
        "fast_root": fast_x,
        "map_history": map_history,
        "run_reference_signature": str(inspect.signature(ref.run_reference)),
        "runner_sha256": _sha(Path(__file__).resolve()),
    }
    _write(record)
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2

    def capture_map(m, freq_values, dn_eff, z_tail, rtol, workers=None,
                    backend="cartesian", z_start=-12.0):
        started = time.perf_counter()
        result = original(
            m, freq_values, dn_eff, z_tail, rtol,
            workers=workers, backend=backend, z_start=z_start,
        )
        elapsed = time.perf_counter() - started
        g2, qerr, ierr = ref.integrate_spectrum(freq_values, result[0], result[1])
        mapped_dn = float(gp.Neff0 * g2 / omega_nu)
        map_history.append({
            "DN_eff_input": float(dn_eff),
            "DN_gw_map": mapped_dn,
            "fixed_point_residual": float(dn_eff - params["DN_eff"] - mapped_dn),
            "quadrature_error": float(qerr),
            "interpolation_error": float(ierr),
            "used_tail_fraction": float(np.mean(result[3])),
            "runtime_s": elapsed,
        })
        record["status"] = "reference_root_iteration_in_progress"
        _write(record)
        return result

    ref._spectrum_at_dn = capture_map
    warnings_seen: list[str] = []
    started = time.perf_counter()
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = ref.run_reference(
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
        ref._spectrum_at_dn = original

    ref_root = float(result["DN_eff"])
    ref_dn_gw = float(result["DN_gw"])
    input_dn = float(params["DN_eff"])
    record.update({
        "status": "complete",
        "reference_root_DN_eff": ref_root,
        "reference_DN_gw_at_root": ref_dn_gw,
        "fixed_point_residual_at_root": ref_root - input_dn - ref_dn_gw,
        "reference_iterations": int(result["n_iter"]),
        "fast_to_reference_root_relative_error": abs(fast_x - ref_root) / max(abs(ref_root), np.finfo(float).tiny),
        "final_quadrature_error": float(result["quadrature_error"]),
        "final_interpolation_error": float(result["interpolation_error"]),
        "warnings": warnings_seen,
        "total_runtime_s": time.perf_counter() - started,
        "resources": telemetry(threads=1, workers=1),
        "source_hashes": {
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "runner.py": _sha(Path(__file__).resolve()),
        },
        "limitation": (
            "Self-consistent Reference root on the validated 299-point adaptive Fast frequency set; "
            "the original z_tail=8 root artifact is a tail-sensitivity diagnostic, not the frozen Reference baseline. "
            "This is not a root enclosure or an independent frequency-convergence proof at the root."
        ),
    })
    _write(record)
    print(json.dumps({
        "status": record["status"],
        "reference_root_DN_eff": ref_root,
        "fast_relative_error": record["fast_to_reference_root_relative_error"],
        "iterations": record["reference_iterations"],
        "map_count": len(map_history),
        "runtime_s": record["total_runtime_s"],
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
