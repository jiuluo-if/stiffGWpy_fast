"""Full-grid Reference validation of the Fast tail8 N_inf>120 candidate."""

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
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARENT = ROOT / "docs" / "fast_ninf122_tail_sensitivity_20261105.json"
CASE = "adaptive_exact_sigma_tail8"
OUTPUT = ROOT / "docs" / "fast_ninf122_tail8_fullgrid_20261105.json"
RTOL = 1e-9
Z_TAIL = 8.0
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
    print(json.dumps({"status": record["status"], "maps": list(record["reference_maps"]), "output": str(OUTPUT)}), flush=True)


def _run_reference(model: LCDM_SG, dn_eff: float, *, freqs=None, freq_res=1.0):
    start = time.perf_counter()
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
    return result, time.perf_counter() - start


def _raw(result: dict) -> np.ndarray:
    return np.asarray(result["Ogw"], dtype=float) - np.asarray(result["Oj"], dtype=float)


def _map_record(result: dict, elapsed: float) -> dict:
    freqs = np.asarray(result["freqs"], dtype=float)
    raw = _raw(result)
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
    fast = parent["results"][CASE]
    params = parent["parameters"]
    dn_eff = float(fast["DN_eff_final"])
    freqs = np.asarray(fast["frequencies_log10_hz"], dtype=float)
    fast_raw = np.asarray(fast["omega_gw"], dtype=float)
    if not (np.all(np.isfinite(freqs)) and np.all(np.isfinite(fast_raw)) and np.all(fast_raw > 0.0)):
        raise FloatingPointError("parent Fast artifact has invalid frequency/spectrum arrays")
    if hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest() != fast["frequency_sha256"]:
        raise RuntimeError("parent frequency array does not match its hash")

    model = LCDM_SG(**params)
    record = {
        "schema_version": 1,
        "experiment": "fast_ninf122_tail8_fixed_background_prufer_fullgrid_and_res2",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": params,
        "N_inf": float(model.derived_param["N_inf"]),
        "fast": {
            "parent_artifact": PARENT.name,
            "parent_artifact_sha256": parent_sha,
            "case": CASE,
            "fixed_DN_eff": dn_eff,
            "DN_gw_self_consistent_output": float(fast["DN_gw"]),
            "frequency_grid_error": float(fast["frequency_grid_error"]),
            "frequency_count": int(freqs.size),
            "frequency_sha256": fast["frequency_sha256"],
            "frequencies_log10_hz": freqs.tolist(),
            "omega_gw": fast_raw.tolist(),
        },
        "reference": {"backend": "prufer", "fixed_DN_eff": dn_eff, "rtol": RTOL, "z_tail": Z_TAIL, "workers": 1},
        "reference_maps": {},
        "status": "adaptive_fast_grid_prufer_running",
        "limitations": [
            "Fixed-background validation only; it does not solve or certify the self-consistent Reference root.",
            "The tail8 option is an explicit hard-point configuration; this does not change or certify the default tail5 profile globally.",
        ],
    }
    _save(record)

    print(f"Prüfer Reference on exact {freqs.size}-node adaptive Fast grid", flush=True)
    same, elapsed = _run_reference(model, dn_eff, freqs=freqs)
    same_record = _map_record(same, elapsed)
    same_raw = np.asarray(same_record["omega_gw"], dtype=float)
    spectrum_error = np.log10(fast_raw / same_raw)
    same_record["fast_vs_prufer_max_abs_dex"] = float(np.max(np.abs(spectrum_error)))
    same_record["fast_vs_prufer_p95_abs_dex"] = float(np.quantile(np.abs(spectrum_error), 0.95))
    same_record["fast_vs_prufer_median_abs_dex"] = float(np.median(np.abs(spectrum_error)))
    same_record["fast_DN_gw_relative_error"] = abs(float(fast["DN_gw"]) - same_record["DN_gw"]) / max(abs(same_record["DN_gw"]), np.finfo(float).tiny)
    record["reference_maps"]["adaptive_fast_grid"] = same_record
    record["status"] = "adaptive_fast_grid_prufer_complete"
    _save(record)

    print(f"Prüfer Reference native freq_res={FREQ_RES:g}", flush=True)
    dense, elapsed = _run_reference(model, dn_eff, freq_res=FREQ_RES)
    dense_record = _map_record(dense, elapsed)
    dense_record["relative_DN_gw_delta_vs_adaptive_fast_grid"] = abs(
        dense_record["DN_gw"] - same_record["DN_gw"]
    ) / max(abs(dense_record["DN_gw"]), abs(same_record["DN_gw"]), np.finfo(float).tiny)
    dense_record["fast_DN_gw_relative_error"] = abs(
        float(fast["DN_gw"]) - dense_record["DN_gw"]
    ) / max(abs(dense_record["DN_gw"]), np.finfo(float).tiny)
    record["reference_maps"]["native_freq_res_2"] = dense_record
    record["acceptance"] = {
        "same_grid_full_spectrum_max_abs_dex_le_1e-3": same_record["fast_vs_prufer_max_abs_dex"] <= 1e-3,
        "same_grid_DN_gw_relative_error_le_2e-4": same_record["fast_DN_gw_relative_error"] <= 2e-4,
        "native_res2_vs_adaptive_reference_DN_gw_delta_le_2e-4": dense_record["relative_DN_gw_delta_vs_adaptive_fast_grid"] <= 2e-4,
        "Fast_DN_gw_vs_native_res2_reference_le_2e-4": dense_record["fast_DN_gw_relative_error"] <= 2e-4,
    }
    record["hashes"] = {
        "fast_sgwb_py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "stiff_SGWB_py": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
        "reference_py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "driver_py": _sha(Path(__file__)),
    }
    record["resources"] = {
        **telemetry(threads=1),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "numba": numba.__version__,
    }
    record["status"] = "complete"
    _save(record)
    print(json.dumps({
        "output": str(OUTPUT),
        "acceptance": record["acceptance"],
        "adaptive_reference": {k: same_record[k] for k in (
            "frequency_count", "DN_gw", "runtime_s", "fast_vs_prufer_max_abs_dex", "fast_DN_gw_relative_error")},
        "native_res2_reference": {k: dense_record[k] for k in (
            "frequency_count", "DN_gw", "runtime_s", "relative_DN_gw_delta_vs_adaptive_fast_grid", "fast_DN_gw_relative_error")},
    }, indent=2))


if __name__ == "__main__":
    main()
