"""Independently validate the opt-in adaptive Fast grid on a standard background."""

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
REFERENCE_RTOL = 1e-9
OUTPUT = ROOT / "docs" / "adaptive_standard_background_validation_20261105.json"


def _save_checkpoint(record: dict) -> None:
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "stage": record["status"],
        "completed_reference_maps": len(record.get("reference_frequency_grid_convergence", [])),
        "output": str(OUTPUT),
    }), flush=True)


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _map(model: LCDM_SG, dn_eff: float, backend: str,
         freqs: np.ndarray | None = None, freq_res: float = 1.0) -> tuple[dict, float]:
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_reference(
            model,
            dn_eff=dn_eff,
            freq_res=freq_res,
            z_tail=Z_TAIL,
            rtol=REFERENCE_RTOL,
            freq_subset=freqs,
            self_consistent=False,
            workers=1,
            physical_guard=False,
            backend=backend,
        )
    elapsed = time.perf_counter() - started
    if caught:
        raise RuntimeError(f"Reference {backend} emitted warnings: {[str(w.message) for w in caught]}")
    return result, elapsed


def main() -> None:
    started_run = time.perf_counter()
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
    fast_raw = np.asarray(fast.Ogw_today, dtype=float) - np.asarray(fast.Oj_today, dtype=float)
    positive = np.isfinite(fast_raw) & (fast_raw > 0.0)
    if np.count_nonzero(positive) < 2:
        raise FloatingPointError("Fast adaptive spectrum has fewer than two positive bins")
    dn_eff = float(fast.cosmo_param["DN_eff"])

    ref_model = LCDM_SG(**PARAMS)
    record = {
        "schema_version": 1,
        "experiment": "adaptive_fast_cr1_standard_background_same_grid_and_reference_convergence",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": PARAMS,
        "fast": {
            "grid": "adaptive",
            "target": TARGET,
            "max_points": MAX_POINTS,
            "sigma_exact": True,
            "z_tail": Z_TAIL,
            "threads": 2,
            "runtime_s": fast_elapsed,
            "frequency_count": int(freqs.size),
            "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
            "frequencies_log10_hz": freqs.tolist(),
            "omega_gw": fast_raw.tolist(),
            "DN_eff_final": dn_eff,
            "DN_gw": float(fast.DN_gw[-1]),
            "derived_n_t": float(fast.derived_param["nt"]),
            "frequency_grid_error": float(fast.freq_grid_error),
            "warnings": [],
        },
        "same_grid_reference": {},
        "reference_frequency_grid_convergence": [],
        "status": "fast_complete",
    }
    _save_checkpoint(record)

    print(f"Reference Cartesian on {len(freqs)} Fast nodes", flush=True)
    cart, cart_s = _map(ref_model, dn_eff, "cartesian", freqs=freqs)
    record["same_grid_reference"]["cartesian"] = {
        "frequencies_log10_hz": freqs.tolist(),
        "omega_gw": (np.asarray(cart["Ogw"], dtype=float) - np.asarray(cart["Oj"], dtype=float)).tolist(),
        "DN_gw": float(cart["DN_gw"]), "runtime_s": cart_s,
        "positive_bins": int(np.count_nonzero((np.asarray(cart["Ogw"]) - np.asarray(cart["Oj"])) > 0.0)),
    }
    record["status"] = "same_grid_cartesian_complete"
    _save_checkpoint(record)

    print("Prüfer independent check on stratified 32-bin subset", flush=True)
    subset_idx = np.unique(np.linspace(0, len(freqs) - 1, min(32, len(freqs)), dtype=int))
    prufer, prufer_s = _map(ref_model, dn_eff, "prufer", freqs=freqs[subset_idx])
    cart_raw = np.asarray(cart["Ogw"], dtype=float) - np.asarray(cart["Oj"], dtype=float)
    prufer_raw = np.asarray(prufer["Ogw"], dtype=float) - np.asarray(prufer["Oj"], dtype=float)
    subset_cart_raw = cart_raw[subset_idx]
    subset_fast_raw = fast_raw[subset_idx]
    common_positive = (subset_cart_raw > 0.0) & (subset_fast_raw > 0.0) & (prufer_raw > 0.0)
    if np.count_nonzero(common_positive) < 2:
        raise FloatingPointError("Independent References have fewer than two common positive bins")
    full_common_positive = positive & np.isfinite(cart_raw) & (cart_raw > 0.0)
    if np.count_nonzero(full_common_positive) < 2:
        raise FloatingPointError("Fast and Cartesian have fewer than two common positive bins")
    cartesian_error_dex = np.log10(
        fast_raw[full_common_positive] / cart_raw[full_common_positive]
    )
    cross_error_dex = np.log10(subset_cart_raw[common_positive] / prufer_raw[common_positive])
    subset_fast_error_dex = np.log10(subset_fast_raw[common_positive] / prufer_raw[common_positive])
    record["same_grid_reference"].update({
        "fixed_DN_eff": dn_eff,
        "rtol": REFERENCE_RTOL,
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "positive_common_bins": int(np.count_nonzero(full_common_positive)),
        "cartesian_vs_fast_max_abs_dex": float(np.max(np.abs(cartesian_error_dex))),
        "cartesian_vs_fast_p95_abs_dex": float(np.quantile(np.abs(cartesian_error_dex), 0.95)),
        "positive_fast_cartesian_bins": int(np.count_nonzero(positive & (cart_raw > 0.0))),
        "prufer_subset": {
            "frequency_indices": subset_idx.tolist(),
            "frequency_count": int(subset_idx.size),
            "DN_gw": float(prufer["DN_gw"]),
            "runtime_s": prufer_s,
            "common_positive_bins": int(np.count_nonzero(common_positive)),
            "max_cartesian_vs_prufer_abs_dex": float(np.max(np.abs(cross_error_dex))),
            "max_fast_vs_prufer_abs_dex": float(np.max(np.abs(subset_fast_error_dex))),
            "omega_gw": prufer_raw.tolist(),
        },
        "cartesian_relative_DN_gw_error_vs_fast": abs(
            float(fast.DN_gw[-1]) - float(cart["DN_gw"])
        ) / max(abs(float(cart["DN_gw"])), np.finfo(float).tiny),
    })
    record["status"] = "same_grid_independent_reference_complete"
    _save_checkpoint(record)

    convergence = []
    previous = None
    for resolution in (2.0, 4.0, 8.0):
        print(f"Reference Cartesian frequency convergence freq_res={resolution:g}", flush=True)
        record["status"] = f"reference_freq_res_{resolution:g}_running"
        _save_checkpoint(record)
        result, elapsed = _map(ref_model, dn_eff, "cartesian", freq_res=resolution)
        row = {
            "freq_res": resolution,
            "frequency_count": int(result["n_freq"]),
            "DN_gw": float(result["DN_gw"]),
            "runtime_s": elapsed,
        }
        if previous is not None:
            row["relative_DN_gw_delta_from_previous"] = abs(
                row["DN_gw"] - previous["DN_gw"]
            ) / max(abs(row["DN_gw"]), abs(previous["DN_gw"]), np.finfo(float).tiny)
        convergence.append(row)
        record["reference_frequency_grid_convergence"] = convergence.copy()
        record["status"] = f"reference_freq_res_{resolution:g}_complete"
        _save_checkpoint(record)
        previous = row
        if resolution >= 2 and row.get("relative_DN_gw_delta_from_previous", np.inf) <= 2e-4:
            break

    cross_max = float(record["same_grid_reference"]["prufer_subset"]["max_cartesian_vs_prufer_abs_dex"])
    cart_max = float(record["same_grid_reference"]["cartesian_vs_fast_max_abs_dex"])
    record["status"] = "complete"
    record["acceptance"] = {
        "same_grid_reference_crosscheck_max_abs_dex_le_1e-3": cross_max <= 1e-3,
        "fast_vs_cartesian_max_abs_dex_le_1e-3": cart_max <= 1e-3,
        "last_reference_grid_relative_delta_le_2e-4": convergence[-1].get("relative_DN_gw_delta_from_previous", np.inf) <= 2e-4,
        "fast_DN_gw_relative_error_le_2e-4": abs(float(fast.DN_gw[-1]) - float(cart["DN_gw"])) / max(abs(float(cart["DN_gw"])), np.finfo(float).tiny) <= 2e-4,
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
        "reference_rtol": REFERENCE_RTOL,
        "wall_elapsed_s": time.perf_counter() - started_run,
    }
    record["limitations"] = [
        "This is one standard positive-cr background and fixed-background Reference cross-check, not a self-consistent Reference root certification or global parameter-domain claim.",
        "Cartesian is compared across every adaptive Fast frequency; Prüfer is an independent 32-bin stratified cross-check.",
    ]
    _save_checkpoint(record)
    print(json.dumps({
        "output": str(OUTPUT),
        "fast_frequency_count": int(freqs.size),
        "fast_DN_gw": float(fast.DN_gw[-1]),
        "acceptance": record["acceptance"],
        "reference_frequency_grid_convergence": convergence,
    }, indent=2))


if __name__ == "__main__":
    main()
