"""Full-grid independent Prüfer check for adaptive Fast on the positive-cr branch."""

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

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 1e-2,
    "n_t": 2.0,
    "cr": 1.0,
    "T_re": 2e3,
    "DN_re": 10.0,
    "kappa10": 1e-2,
}
TARGET = 1e-4
MAX_POINTS = 1500
Z_TAIL = 10.0
RTOL = 1e-9
EXPECTED_FREQ_SHA = "da027a66450eaedbc99aa7da4d521f434f70df6b236f9f68211ad80751bbabae"
OUTPUT = ROOT / "docs" / "adaptive_cr1_prufer_fullgrid_20261105.json"


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(record: dict) -> None:
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": record["status"],
        "completed_maps": list(record["reference_maps"]),
        "output": str(OUTPUT),
    }), flush=True)


def _reference(model: LCDM_SG, dn_eff: float, *, freqs=None, freq_res=1.0):
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        output = run_reference(
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
    return output, time.perf_counter() - started


def _spectral_record(freqs: np.ndarray, raw: np.ndarray) -> dict:
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(raw)):
        raise FloatingPointError("cannot archive non-finite frequency or spectrum values")
    return {
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "frequencies_log10_hz": freqs.tolist(),
        "omega_gw": raw.tolist(),
    }


def main() -> None:
    FS.apply_accuracy_mode("fast")
    fast = LCDM_SG(**PARAMS)
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        output = fast.SGWB_iter(
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
    if output is not fast:
        raise RuntimeError(f"Fast adaptive solve failed: {fast.fast_failure_reason}")
    if caught:
        raise RuntimeError(f"Fast emitted warnings: {[str(w.message) for w in caught]}")
    freqs = np.asarray(fast.f, dtype=float)
    freq_sha = hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest()
    if freq_sha != EXPECTED_FREQ_SHA:
        raise RuntimeError(f"Fast grid changed: {freq_sha} != {EXPECTED_FREQ_SHA}")
    fast_raw = np.asarray(fast.Ogw_today, dtype=float) - np.asarray(fast.Oj_today, dtype=float)
    if not np.all(np.isfinite(fast_raw)):
        raise FloatingPointError("Fast spectrum contains non-finite values")
    dn_eff = float(fast.cosmo_param["DN_eff"])
    fast_record = {
        **_spectral_record(freqs, fast_raw),
        "runtime_s": fast_elapsed,
        "DN_eff_final": dn_eff,
        "DN_gw": float(fast.DN_gw[-1]),
        "frequency_grid_error": float(fast.freq_grid_error),
        "derived_n_t": float(fast.derived_param["nt"]),
        "warnings": [],
    }
    record = {
        "schema_version": 1,
        "experiment": "adaptive_cr1_positive_cr_fullgrid_prufer_and_frequency_convergence",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": PARAMS,
        "fast": fast_record,
        "reference": {
            "backend": "prufer",
            "rtol": RTOL,
            "z_tail": Z_TAIL,
            "workers": 1,
            "fixed_DN_eff": dn_eff,
        },
        "reference_maps": {},
        "status": "fast_complete",
        "limitations": [
            "Fixed-background validation only; it does not independently solve the self-consistent Reference root.",
            "A separate Cartesian full-grid result and 32-bin Prüfer cross-check are archived in adaptive_standard_background_validation_20261105.json.",
        ],
    }
    _save(record)

    ref_model = LCDM_SG(**PARAMS)
    print(f"Prüfer Reference on all {len(freqs)} adaptive Fast frequencies", flush=True)
    same, elapsed = _reference(ref_model, dn_eff, freqs=freqs)
    ref_raw = np.asarray(same["Ogw"], dtype=float) - np.asarray(same["Oj"], dtype=float)
    common = np.isfinite(fast_raw) & np.isfinite(ref_raw) & (fast_raw > 0.0) & (ref_raw > 0.0)
    if not np.all(common):
        raise FloatingPointError(f"positive spectrum support differs on {int(np.count_nonzero(~common))} Fast nodes")
    error_dex = np.log10(fast_raw[common] / ref_raw[common])
    record["reference_maps"]["adaptive_fast_grid"] = {
        **_spectral_record(freqs, ref_raw),
        "DN_gw": float(same["DN_gw"]),
        "runtime_s": elapsed,
        "used_tail_count": int(np.count_nonzero(same["used_tail"])),
        "max_fast_abs_dex": float(np.max(np.abs(error_dex))),
        "p95_fast_abs_dex": float(np.quantile(np.abs(error_dex), 0.95)),
        "median_fast_abs_dex": float(np.median(np.abs(error_dex))),
        "relative_DN_gw_error_vs_fast": abs(float(fast.DN_gw[-1]) - float(same["DN_gw"])) / max(abs(float(same["DN_gw"])), np.finfo(float).tiny),
    }
    record["status"] = "adaptive_grid_prufer_complete"
    _save(record)

    print("Prüfer Reference frequency-grid convergence at freq_res=4", flush=True)
    record["status"] = "native_freq_res_4_running"
    _save(record)
    dense, elapsed = _reference(ref_model, dn_eff, freq_res=4.0)
    dense_freqs = np.asarray(dense["freqs"], dtype=float)
    dense_raw = np.asarray(dense["Ogw"], dtype=float) - np.asarray(dense["Oj"], dtype=float)
    adaptive_dn = float(same["DN_gw"])
    grid_delta = abs(float(dense["DN_gw"]) - adaptive_dn) / max(
        abs(float(dense["DN_gw"])), abs(adaptive_dn), np.finfo(float).tiny
    )
    record["reference_maps"]["native_freq_res_4"] = {
        **_spectral_record(dense_freqs, dense_raw),
        "DN_gw": float(dense["DN_gw"]),
        "runtime_s": elapsed,
        "used_tail_count": int(np.count_nonzero(dense["used_tail"])),
        "relative_DN_gw_delta_vs_adaptive_fast_grid_reference": grid_delta,
    }
    record["status"] = "native_freq_res_4_complete"
    _save(record)

    record["acceptance"] = {
        "prufer_same_grid_spectrum_max_abs_dex_le_1e-3": record["reference_maps"]["adaptive_fast_grid"]["max_fast_abs_dex"] <= 1e-3,
        "prufer_same_grid_DN_gw_relative_error_le_2e-4": record["reference_maps"]["adaptive_fast_grid"]["relative_DN_gw_error_vs_fast"] <= 2e-4,
        "native_res4_vs_adaptive_reference_DN_gw_relative_delta_le_2e-4": grid_delta <= 2e-4,
    }
    record["hashes"] = {
        "fast_sgwb_py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "stiff_SGWB_py": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
        "reference_py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "script_py": _sha(Path(__file__)),
    }
    record["resources"] = {
        **telemetry(threads=2),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "workers": 1,
    }
    record["status"] = "complete"
    _save(record)
    print(json.dumps({
        "output": str(OUTPUT),
        "acceptance": record["acceptance"],
        "adaptive_reference": {
            key: record["reference_maps"]["adaptive_fast_grid"][key]
            for key in ("frequency_count", "DN_gw", "runtime_s", "max_fast_abs_dex", "relative_DN_gw_error_vs_fast")
        },
        "freq_res4_reference": {
            key: record["reference_maps"]["native_freq_res_4"][key]
            for key in ("frequency_count", "DN_gw", "runtime_s", "relative_DN_gw_delta_vs_adaptive_fast_grid_reference")
        },
    }, indent=2))


if __name__ == "__main__":
    main()
