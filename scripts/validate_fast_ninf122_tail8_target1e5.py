"""Validate the N_inf>120 Fast tail8 target=1e-5 candidate with Prüfer maps."""

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

PARENT = ROOT / "docs" / "fast_ninf122_tail_sensitivity_20261105.json"
OUTPUT = ROOT / "docs" / "fast_ninf122_tail8_target1e5_validation_20261105.json"
TARGET = 1e-5
MAX_POINTS = 1500
Z_TAIL = 8.0
RTOL = 1e-9
FREQ_RES = 2.0


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(record: dict) -> None:
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": record["status"], "maps": list(record.get("reference_maps", {})), "output": str(OUTPUT)}), flush=True)


def _reference(model: LCDM_SG, dn_eff: float, *, freqs=None, freq_res=1.0):
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
        raise RuntimeError(f"Prüfer Reference emitted warnings: {[str(w.message) for w in caught]}")
    return result, time.perf_counter() - started


def _map_record(result: dict, elapsed: float) -> dict:
    freqs = np.asarray(result["freqs"], dtype=float)
    raw = np.asarray(result["Ogw"], dtype=float) - np.asarray(result["Oj"], dtype=float)
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(raw)):
        raise FloatingPointError("Reference returned non-finite frequency or spectrum values")
    return {
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "frequencies_log10_hz": freqs.tolist(),
        "omega_gw": raw.tolist(),
        "DN_gw": float(result["DN_gw"]),
        "quadrature_error": float(result["quadrature_error"]),
        "interpolation_error": float(result["interpolation_error"]),
        "used_tail_count": int(np.count_nonzero(result["used_tail"])),
        "runtime_s": elapsed,
    }


def main() -> None:
    parent_sha = _sha(PARENT)
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    params = parent["parameters"]
    FS.apply_accuracy_mode("fast")
    fast_model = LCDM_SG(**params)
    fast_start = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = fast_model.SGWB_iter(
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
    fast_elapsed = time.perf_counter() - fast_start
    if out is not fast_model:
        raise RuntimeError(f"Fast candidate failed: {fast_model.fast_failure_reason}")
    if caught:
        raise RuntimeError(f"Fast emitted warnings: {[str(w.message) for w in caught]}")
    freqs = np.asarray(fast_model.f, dtype=float)
    fast_raw = np.asarray(fast_model.Ogw_today, dtype=float) - np.asarray(fast_model.Oj_today, dtype=float)
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(fast_raw)) or not np.all(fast_raw > 0.0):
        raise FloatingPointError("Fast candidate returned invalid spectrum data")
    freq_hash = hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest()
    dn_eff = float(fast_model.cosmo_param["DN_eff"])
    budget = getattr(fast_model, "local_error_budget", None)
    record = {
        "schema_version": 1,
        "experiment": "fast_ninf122_tail8_adaptive_target1e5_prufer_fullgrid_and_res2",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": params,
        "parent_sensitivity_artifact": PARENT.name,
        "parent_sensitivity_artifact_sha256": parent_sha,
        "N_inf": float(fast_model.derived_param["N_inf"]),
        "fast": {
            "grid": "adaptive",
            "target": TARGET,
            "max_points": MAX_POINTS,
            "sigma_exact": True,
            "z_tail": Z_TAIL,
            "threads": 2,
            "runtime_s": fast_elapsed,
            "frequency_count": int(freqs.size),
            "frequency_sha256": freq_hash,
            "frequencies_log10_hz": freqs.tolist(),
            "omega_gw": fast_raw.tolist(),
            "DN_eff_final": dn_eff,
            "DN_gw": float(fast_model.DN_gw[-1]),
            "frequency_grid_error": float(fast_model.freq_grid_error),
            "local_error_status": budget.get("certification_status") if isinstance(budget, dict) else None,
            "warnings": [],
        },
        "reference": {"backend": "prufer", "rtol": RTOL, "z_tail": Z_TAIL, "fixed_DN_eff": dn_eff, "workers": 1},
        "reference_maps": {},
        "status": "adaptive_exact_grid_running",
    }
    _save(record)

    ref_model = LCDM_SG(**params)
    print(f"Prüfer Reference on exact {freqs.size}-node Fast target={TARGET:g} grid", flush=True)
    exact, elapsed = _reference(ref_model, dn_eff, freqs=freqs)
    exact_record = _map_record(exact, elapsed)
    exact_raw = np.asarray(exact_record["omega_gw"], dtype=float)
    err = np.log10(fast_raw / exact_raw)
    exact_record.update({
        "max_fast_abs_dex": float(np.max(np.abs(err))),
        "p95_fast_abs_dex": float(np.quantile(np.abs(err), 0.95)),
        "median_fast_abs_dex": float(np.median(np.abs(err))),
        "fast_DN_gw_relative_error": abs(float(fast_model.DN_gw[-1]) - exact_record["DN_gw"]) / max(abs(exact_record["DN_gw"]), np.finfo(float).tiny),
    })
    record["reference_maps"]["adaptive_fast_grid"] = exact_record
    record["status"] = "adaptive_exact_grid_complete"
    _save(record)

    print(f"Prüfer Reference native freq_res={FREQ_RES:g}", flush=True)
    record["status"] = "native_res2_running"
    _save(record)
    native, elapsed = _reference(ref_model, dn_eff, freq_res=FREQ_RES)
    native_record = _map_record(native, elapsed)
    native_record["relative_DN_gw_delta_vs_adaptive_reference"] = abs(
        native_record["DN_gw"] - exact_record["DN_gw"]
    ) / max(abs(native_record["DN_gw"]), abs(exact_record["DN_gw"]), np.finfo(float).tiny)
    native_record["fast_DN_gw_relative_error"] = abs(
        float(fast_model.DN_gw[-1]) - native_record["DN_gw"]
    ) / max(abs(float(fast_model.DN_gw[-1])), abs(native_record["DN_gw"]), np.finfo(float).tiny)
    record["reference_maps"]["native_res2"] = native_record
    record["acceptance"] = {
        "adaptive_full_spectrum_max_abs_dex_le_1e-3": exact_record["max_fast_abs_dex"] <= 1e-3,
        "adaptive_map_DN_gw_rel_le_2e-4": exact_record["fast_DN_gw_relative_error"] <= 2e-4,
        "native_res2_grid_delta_le_2e-4": native_record["relative_DN_gw_delta_vs_adaptive_reference"] <= 2e-4,
        "Fast_DN_gw_vs_native_res2_le_2e-4": native_record["fast_DN_gw_relative_error"] <= 2e-4,
    }
    record["hashes"] = {
        "fast_sgwb_py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "stiff_SGWB_py": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
        "reference_py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "driver_py": _sha(Path(__file__)),
    }
    record["resources"] = {
        **telemetry(threads=2),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "numba": numba.__version__,
        "fast_threads": 2,
        "reference_workers": 1,
    }
    record["limitations"] = [
        "Fixed-background verification only; it does not solve or certify the self-consistent Reference root.",
        "The res2 comparison is a two-grid convergence screen, not a formal quadrature error bound.",
        "This is one N_inf>120 point and does not establish a universal tail or grid dispatch threshold.",
    ]
    record["status"] = "complete"
    _save(record)
    print(json.dumps({
        "output": str(OUTPUT),
        "N_inf": record["N_inf"],
        "acceptance": record["acceptance"],
        "adaptive_reference": {k: exact_record[k] for k in (
            "frequency_count", "DN_gw", "runtime_s", "max_fast_abs_dex", "fast_DN_gw_relative_error")},
        "native_res2_reference": {k: native_record[k] for k in (
            "frequency_count", "DN_gw", "runtime_s", "relative_DN_gw_delta_vs_adaptive_reference", "fast_DN_gw_relative_error")},
    }, indent=2))


if __name__ == "__main__":
    main()
