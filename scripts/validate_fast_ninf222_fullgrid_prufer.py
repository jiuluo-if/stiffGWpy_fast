"""Full-grid Prüfer and native-grid checks for the N_inf=222.5 Fast point."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numba
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARENT = ROOT / "docs" / "fast_ninf122_tail8_target1e5_validation_20261105.json"
OUTPUT = ROOT / "docs" / "fast_ninf222_fullgrid_prufer_20261105.json"
TARGET = 1e-5
MAX_POINTS = 1500
Z_TAIL = 8.0
RTOL = 1e-9
NATIVE_FREQ_RES = 2.0


def _git(*args: str) -> str | None:
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(record: dict) -> None:
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": record["status"],
        "reference_maps": list(record["reference_maps"]),
        "output": str(OUTPUT),
    }), flush=True)


def _reference(model: LCDM_SG, dn_eff: float, *, freqs=None, freq_res=1.0) -> tuple[dict, float]:
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_reference(
            model,
            dn_eff=dn_eff,
            freq_res=freq_res,
            z_tail=Z_TAIL,
            rtol=RTOL,
            freq_subset=freqs,
            self_consistent=False,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
    if caught:
        raise RuntimeError(f"Prüfer Reference emitted warnings: {[str(item.message) for item in caught]}")
    return result, time.perf_counter() - started


def _map_record(result: dict, elapsed: float) -> dict:
    freqs = np.asarray(result["freqs"], dtype=float)
    spectrum = np.asarray(result["Ogw"], dtype=float) - np.asarray(result["Oj"], dtype=float)
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(spectrum)) or not np.all(spectrum > 0.0):
        raise FloatingPointError("Prüfer Reference returned invalid frequency or spectrum values")
    return {
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "frequencies_log10_hz": freqs.tolist(),
        "omega_gw": spectrum.tolist(),
        "DN_gw": float(result["DN_gw"]),
        "quadrature_error": float(result["quadrature_error"]),
        "interpolation_error": float(result["interpolation_error"]),
        "runtime_s": elapsed,
        "used_tail_count": int(np.count_nonzero(result["used_tail"])),
    }


def main() -> None:
    parent_sha = _sha(PARENT)
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    params = {**parent["parameters"], "DN_re": 200.0}
    FS.apply_accuracy_mode("fast")
    model = LCDM_SG(**params)
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=1e-6,
            z_tail=Z_TAIL,
            freq_res=1.0,
            threads=2,
            accuracy_mode="fast",
            physical_guard=False,
            freq_grid="adaptive",
            freq_grid_target=TARGET,
            freq_grid_max_points=MAX_POINTS,
            sigma_exact=True,
        )
    fast_elapsed = time.perf_counter() - started
    if out is not model or not model.SGWB_converge:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    if caught:
        raise RuntimeError(f"Fast emitted warnings: {[str(item.message) for item in caught]}")

    freqs = np.asarray(model.f, dtype=float)
    fast_spectrum = np.asarray(model.Ogw_today, dtype=float) - np.asarray(model.Oj_today, dtype=float)
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(fast_spectrum)) or not np.all(fast_spectrum > 0.0):
        raise FloatingPointError("Fast returned invalid frequencies or spectrum values")
    parent_freqs = np.asarray(parent["fast"]["frequencies_log10_hz"], dtype=float)
    if not np.array_equal(freqs, parent_freqs):
        raise AssertionError("candidate Fast frequency array differs from parent grid")
    dn_eff = float(model.cosmo_param["DN_eff"])

    record = {
        "schema_version": 1,
        "experiment": "fast_ninf222_tail8_adaptive_prufer_fullgrid_and_native_res2",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": params,
        "N_inf": float(model.derived_param["N_inf"]),
        "parent_artifact": PARENT.name,
        "parent_artifact_sha256": parent_sha,
        "equation": "fixed-background map at Fast DN_eff_final; self-consistent root sensitivity is separate",
        "fast": {
            "grid": "adaptive",
            "target_dex": TARGET,
            "max_points": MAX_POINTS,
            "sigma_exact": True,
            "z_tail": Z_TAIL,
            "threads": 2,
            "runtime_s": fast_elapsed,
            "frequency_count": int(freqs.size),
            "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
            "frequencies_log10_hz": freqs.tolist(),
            "omega_gw": fast_spectrum.tolist(),
            "DN_eff_final": dn_eff,
            "DN_gw": float(model.DN_gw[-1]),
            "frequency_grid_error": float(model.freq_grid_error),
            "warnings": [],
        },
        "reference": {
            "backend": "prufer",
            "rtol": RTOL,
            "z_tail": Z_TAIL,
            "workers": 1,
            "fixed_DN_eff": dn_eff,
            "native_freq_res": NATIVE_FREQ_RES,
        },
        "reference_maps": {},
        "resources": {
            **telemetry(threads=2),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "numba": numba.__version__,
            "fast_threads": 2,
            "reference_workers": 1,
        },
        "hashes": {
            "fast_sgwb_py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference_py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "script_py": _sha(Path(__file__)),
        },
        "status": "fast_complete_same_grid_running",
    }
    _save(record)

    print("Prüfer Reference on exact Fast frequency array", flush=True)
    exact_result, exact_elapsed = _reference(model, dn_eff, freqs=freqs, freq_res=1.0)
    exact = _map_record(exact_result, exact_elapsed)
    if exact["frequency_sha256"] != record["fast"]["frequency_sha256"]:
        raise AssertionError("same-grid Prüfer frequencies differ from Fast")
    exact_error = np.log10(fast_spectrum / np.asarray(exact["omega_gw"], dtype=float))
    exact["fast_max_abs_dex"] = float(np.max(np.abs(exact_error)))
    exact["fast_p95_abs_dex"] = float(np.quantile(np.abs(exact_error), 0.95))
    exact["fast_DN_gw_relative_error"] = abs(float(model.DN_gw[-1]) - exact["DN_gw"]) / max(
        abs(float(model.DN_gw[-1])), abs(exact["DN_gw"]), np.finfo(float).tiny
    )
    parent_ref = np.asarray(parent["reference_maps"]["adaptive_fast_grid"]["omega_gw"], dtype=float)
    exact["candidate_prufer_vs_parent_prufer_max_abs_dex"] = float(
        np.max(np.abs(np.log10(np.asarray(exact["omega_gw"], dtype=float) / parent_ref)))
    )
    record["reference_maps"]["adaptive_exact_fast_grid"] = exact
    record["status"] = "fast_complete_same_grid_complete_native_running"
    _save(record)

    print("Prüfer Reference on native freq_res=2 grid", flush=True)
    native_result, native_elapsed = _reference(model, dn_eff, freq_res=NATIVE_FREQ_RES)
    native = _map_record(native_result, native_elapsed)
    native["relative_DN_gw_delta_vs_adaptive_reference"] = abs(native["DN_gw"] - exact["DN_gw"]) / max(
        abs(native["DN_gw"]), abs(exact["DN_gw"]), np.finfo(float).tiny
    )
    native["fast_DN_gw_relative_error"] = abs(float(model.DN_gw[-1]) - native["DN_gw"]) / max(
        abs(float(model.DN_gw[-1])), abs(native["DN_gw"]), np.finfo(float).tiny
    )
    record["reference_maps"]["native_res2"] = native
    record["acceptance"] = {
        "full_spectrum_max_abs_dex_le_1e-3": exact["fast_max_abs_dex"] <= 1e-3,
        "same_grid_DN_gw_rel_le_2e-4": exact["fast_DN_gw_relative_error"] <= 2e-4,
        "native_grid_delta_le_2e-4": native["relative_DN_gw_delta_vs_adaptive_reference"] <= 2e-4,
        "fast_vs_native_DN_gw_rel_le_2e-4": native["fast_DN_gw_relative_error"] <= 2e-4,
    }
    record["status"] = "fixed_background_spectrum_and_grid_screen_complete"
    record["limitations"] = [
        "This verifies a fixed background only; it does not solve the independent Reference self-consistent root.",
        "The map checks do not prove invariance under DN_re shifts outside this tested pair of points.",
    ]
    _save(record)
    print(json.dumps({
        "status": record["status"],
        "N_inf": record["N_inf"],
        "fast_spectrum_max_abs_dex": exact["fast_max_abs_dex"],
        "fast_same_grid_DN_rel": exact["fast_DN_gw_relative_error"],
        "native_grid_delta": native["relative_DN_gw_delta_vs_adaptive_reference"],
        "fast_native_DN_rel": native["fast_DN_gw_relative_error"],
        "acceptance": record["acceptance"],
        "output": str(OUTPUT),
    }), flush=True)


if __name__ == "__main__":
    main()
