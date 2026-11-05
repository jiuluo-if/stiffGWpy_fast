"""Stratified Prüfer check for a high-r strongly red positive-cr Fast point."""

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
OUTPUT = ROOT / "docs" / "fast_cr1_red_extreme_modes_20261105.json"
PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 10.0,
    "n_t": 2.0,
    "cr": 1.0,
    "T_re": 1e-3,
    "DN_re": 200.0,
    "kappa10": 100.0,
}
TARGET = 1e-7
MAX_POINTS = 1500
Z_TAIL = 8.0
RTOL = 1e-9
MODE_COUNT = 12
SPECTRUM_TOL = 1e-3
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _git(*args: str) -> str | None:
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    FS.apply_accuracy_mode("fast")
    model = LCDM_SG(**PARAMS)
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
    f = np.asarray(model.f, dtype=float)
    y = np.asarray(model.Ogw_today, dtype=float) - np.asarray(model.Oj_today, dtype=float)
    if not np.all(np.isfinite(f)) or not np.all(np.isfinite(y)) or not np.all(y > 0.0):
        raise FloatingPointError("Fast returned invalid spectrum values")
    indices = np.linspace(0, f.size - 1, MODE_COUNT).astype(int)
    selected = f[indices]
    dn_eff = float(model.cosmo_param["DN_eff"])
    derived = model.derived_param

    reference_started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_reference(
            model,
            dn_eff=dn_eff,
            freq_res=1.0,
            z_tail=Z_TAIL,
            rtol=RTOL,
            freq_subset=selected,
            self_consistent=False,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
    reference_elapsed = time.perf_counter() - reference_started
    if caught:
        raise RuntimeError(f"Prüfer emitted warnings: {[str(item.message) for item in caught]}")
    rf = np.asarray(result["freqs"], dtype=float)
    ry = np.asarray(result["Ogw"], dtype=float) - np.asarray(result["Oj"], dtype=float)
    if not np.array_equal(rf, selected) or not np.all(np.isfinite(ry)) or not np.all(ry > 0.0):
        raise FloatingPointError("Prüfer selected-mode output is invalid or grid-mismatched")
    error = np.log10(y[indices] / ry)
    max_error = float(np.max(np.abs(error)))

    record = {
        "schema_version": 1,
        "experiment": "fast_cr1_high_r_red_extreme_selected_prufer_screen",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters_input": PARAMS,
        "derived": {
            "N_inf": float(derived["N_inf"]),
            "nt_effective": float(derived["nt"]),
            "DN_re_effective": float(model.cosmo_param["DN_re"]),
        },
        "fixed_DN_eff": dn_eff,
        "fast": {
            "grid": "adaptive",
            "target_dex": TARGET,
            "z_tail": Z_TAIL,
            "frequency_count": int(f.size),
            "frequency_sha256": hashlib.sha256(f.astype("<f8").tobytes()).hexdigest(),
            "frequencies_log10_hz": f.tolist(),
            "omega_gw": y.tolist(),
            "DN_eff_final": dn_eff,
            "DN_gw": float(model.DN_gw[-1]),
            "frequency_grid_error": float(model.freq_grid_error),
            "runtime_s": fast_elapsed,
        },
        "reference": {
            "backend": "prufer",
            "rtol": RTOL,
            "z_tail": Z_TAIL,
            "workers": 1,
            "selected_indices": indices.tolist(),
            "frequency_sha256": hashlib.sha256(selected.astype("<f8").tobytes()).hexdigest(),
            "frequencies_log10_hz": selected.tolist(),
            "omega_gw": ry.tolist(),
            "runtime_s": reference_elapsed,
        },
        "comparison": {
            "selected_mode_count": int(indices.size),
            "fast_vs_prufer_max_abs_dex": max_error,
            "fast_vs_prufer_p95_abs_dex": float(np.quantile(np.abs(error), 0.95)),
            "spectrum_gate_1e-3_dex": max_error <= SPECTRUM_TOL,
            "finite_positive_fast_and_reference": True,
        },
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
            "fast_sgwb_py": _sha(ROOT / "stiffgwpy_fast/fast_sgwb.py"),
            "reference_py": _sha(ROOT / "stiffgwpy_fast/reference.py"),
            "script_py": _sha(Path(__file__)),
        },
        "limitations": [
            "This is a 12-mode fixed-background preflight, not a full-spectrum or bolometric convergence result.",
            "The positive-cr branch derives nt and DN_re; input values are not independent coordinates.",
        ],
        "status": "selected_mode_screen_passed_fixed_background_only" if max_error <= SPECTRUM_TOL else "selected_mode_screen_failed",
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": record["status"], "derived": record["derived"], "comparison": record["comparison"], "output": str(OUTPUT)}), flush=True)


if __name__ == "__main__":
    main()
