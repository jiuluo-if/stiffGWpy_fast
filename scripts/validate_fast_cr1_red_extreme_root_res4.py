"""Check the iterative Reference root against a native res4 Prüfer map."""

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

from stiffgwpy_fast import reference as ref  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

ROOT_ARTIFACT = ROOT / "docs" / "fast_cr1_red_extreme_reference_root_20261105.json"
FAST_ARTIFACT = ROOT / "docs" / "fast_cr1_red_extreme_fullgrid_prufer_20261105.json"
RES8_ARTIFACT = ROOT / "docs" / "fast_cr1_red_extreme_res8_prufer_20261105.json"
OUTPUT = ROOT / "docs" / "fast_cr1_red_extreme_root_res4_20261105.json"


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
    root_data = json.loads(ROOT_ARTIFACT.read_text(encoding="utf-8"))
    fast_data = json.loads(FAST_ARTIFACT.read_text(encoding="utf-8"))
    res8_data = json.loads(RES8_ARTIFACT.read_text(encoding="utf-8"))
    model = LCDM_SG(**root_data["parameters"])
    x = float(root_data["reference_root_DN_eff"])
    _write({
        "schema_version": 1,
        "experiment": "native res4 Prüfer map at the iterative high-r red positive-cr root",
        "status": "running",
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "reference_root_artifact": ROOT_ARTIFACT.name,
        "reference_root_artifact_sha256": _sha(ROOT_ARTIFACT),
        "fixed_DN_eff": x,
        "reference": {"backend": "prufer", "freq_res": 4.0, "z_tail": 8.0, "rtol": 1e-9, "workers": 1},
    })
    started = time.perf_counter()
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = ref.run_reference(
                model,
                dn_eff=x,
                freq_res=4.0,
                self_consistent=False,
                z_tail=8.0,
                rtol=1e-9,
                workers=1,
                physical_guard=False,
                backend="prufer",
            )
    except Exception as exc:
        _write({
            "schema_version": 1,
            "experiment": "native res4 Prüfer map at the iterative high-r red positive-cr root",
            "status": "failed",
            "reference_root_artifact": ROOT_ARTIFACT.name,
            "fixed_DN_eff": x,
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
            "runtime_s": time.perf_counter() - started,
        })
        raise
    runtime = time.perf_counter() - started
    warnings_seen = [str(item.message) for item in caught]
    freqs = np.asarray(result["freqs"], dtype=np.float64)
    omega = np.asarray(result["Ogw"], dtype=np.float64) - np.asarray(result["Oj"], dtype=np.float64)
    input_dn = float(root_data["parameters"]["DN_eff"])
    dn_map = float(result["DN_gw"])
    fast_x = float(root_data["fast_root"])
    output = {
        "schema_version": 1,
        "experiment": "native res4 Prüfer map at the iterative high-r red positive-cr root",
        "status": "complete",
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": root_data["parameters"],
        "reference_root_artifact": ROOT_ARTIFACT.name,
        "reference_root_artifact_sha256": _sha(ROOT_ARTIFACT),
        "fixed_DN_eff": x,
        "input_DN_eff": input_dn,
        "map_DN_gw": dn_map,
        "fixed_point_residual_at_iterative_root": float(x - input_dn - dn_map),
        "relative_root_residual": abs(x - input_dn - dn_map) / max(abs(x), np.finfo(float).tiny),
        "fast_root": fast_x,
        "fast_to_iterative_reference_root_relative_delta": abs(fast_x - x) / max(abs(x), np.finfo(float).tiny),
        "native_res4_at_fast_root_DN_gw": float(fast_data["reference_maps"]["native_res4"]["DN_gw"]),
        "native_res8_at_fast_root_DN_gw": float(res8_data["native_res8"]["DN_gw"]),
        "res8_map_root_anchor_delta": float(x - fast_x),
        "reference": {"backend": "prufer", "freq_res": 4.0, "z_tail": 8.0, "rtol": 1e-9, "workers": 1},
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "frequencies_log10_hz": freqs.tolist(),
        "omega_gw": omega.tolist(),
        "log10_omega_gw": np.log10(omega).tolist(),
        "quadrature_error": float(result["quadrature_error"]),
        "interpolation_error": float(result["interpolation_error"]),
        "used_tail_fraction": float(np.mean(result["used_tail"])),
        "runtime_s": runtime,
        "warnings": warnings_seen,
        "resources": telemetry(threads=1, workers=1),
        "source_hashes": {
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "runner.py": _sha(Path(__file__).resolve()),
        },
        "limitation": (
            "This is one native res4 map at the iterative root. The 618-second native res8 map "
            "was reused at the Fast root, which differs from this root by only the recorded delta; "
            "no res8 propagation was repeated at the iterative root. This is a residual/convergence "
            "screen, not a rigorous root enclosure."
        ),
    }
    _write(output)
    print(json.dumps({
        "status": output["status"],
        "DN_gw_res4_at_root": dn_map,
        "root_residual_rel": output["relative_root_residual"],
        "frequency_count": freqs.size,
        "runtime_s": runtime,
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
