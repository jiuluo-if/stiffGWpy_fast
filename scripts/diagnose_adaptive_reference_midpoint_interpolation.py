"""Check Reference interpolation at actual Fast midpoint frequencies.

The Fast midpoint probes use nested frequency midpoints, while archived dense
Prüfer grids are nonnested. This diagnostic solves a small set of the exact
midpoint frequencies with the unchanged Reference backend to separate
interpolant error from Fast-vs-Reference propagation error.
"""

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
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, telemetry  # noqa: E402

apply_environment()

from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = (
    "adaptive_embedded_fast_cr1_lowT_20261105.json",
    "adaptive_embedded_fast_cr1_red_extreme_20261105.json",
)
OUTPUT = ROOT / "docs" / "adaptive_reference_midpoint_interpolation_20261105.json"
RTOL = 1e-9
TOP_COUNT = 8


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _case(path: Path) -> dict:
    probe = json.loads(path.read_text(encoding="utf-8"))
    source_path = ROOT / "docs" / probe["source_artifact"]
    source = json.loads(source_path.read_text(encoding="utf-8"))
    params = source.get("parameters", source.get("parameters_input"))
    if not isinstance(params, dict):
        raise RuntimeError(f"missing parameters in {source_path.name}")
    dn_eff = float(source.get("fixed_DN_eff", source["fast"]["DN_eff_final"]))
    z_tail = float(source["reference"]["z_tail"])

    freqs = np.asarray(
        probe["actual_fast_midpoint_samples"]["frequencies_log10_hz"],
        dtype=np.float64,
    )
    fast_mid = np.asarray(probe["actual_fast_midpoint_samples"]["omega_gw"], dtype=np.float64)
    reference_grid = source["reference_maps"][probe["dense_reference_key"]]
    ref_freqs = np.asarray(reference_grid["frequencies_log10_hz"], dtype=np.float64)
    ref_spectrum = np.asarray(reference_grid["omega_gw"], dtype=np.float64)
    order = np.argsort(ref_freqs)
    dense_interp = np.asarray(
        PchipInterpolator(
            ref_freqs[order], np.log10(ref_spectrum[order]), extrapolate=False
        )(freqs),
        dtype=np.float64,
    )
    interpolated_gap = np.abs(np.log10(fast_mid) - dense_interp)
    valid = np.flatnonzero(np.isfinite(interpolated_gap))
    if valid.size != freqs.size:
        raise RuntimeError("dense Reference interpolation does not cover every Fast midpoint")

    top = valid[np.argsort(interpolated_gap[valid])[-TOP_COUNT:]]
    quantiles = np.array(
        [0, freqs.size // 4, freqs.size // 2, 3 * freqs.size // 4, freqs.size - 1],
        dtype=np.int64,
    )
    selected = np.unique(np.concatenate((top, quantiles)))
    model = LCDM_SG(**params)
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_reference(
            model,
            dn_eff=dn_eff,
            freq_subset=freqs[selected],
            z_tail=z_tail,
            rtol=RTOL,
            self_consistent=False,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
    elapsed = time.perf_counter() - started
    if caught:
        raise RuntimeError(f"Reference emitted warnings: {[str(item.message) for item in caught]}")
    exact_ref = np.asarray(result["Ogw"], dtype=np.float64) - np.asarray(
        result["Oj"], dtype=np.float64
    )
    if not np.all(np.isfinite(exact_ref)) or not np.all(exact_ref > 0.0):
        raise FloatingPointError("exact midpoint Prüfer spectrum is not finite and positive")
    fast_selected = fast_mid[selected]
    interp_selected = np.power(10.0, dense_interp[selected])
    fast_vs_exact = np.abs(np.log10(fast_selected / exact_ref))
    exact_vs_interp = np.abs(np.log10(exact_ref / interp_selected))
    return {
        "source_probe": path.name,
        "source_probe_sha256": _sha(path),
        "source_reference": source_path.name,
        "source_reference_sha256": _sha(source_path),
        "parameters": params,
        "fixed_DN_eff": dn_eff,
        "z_tail": z_tail,
        "reference_backend": "prufer",
        "reference_rtol": RTOL,
        "dense_reference_grid": probe["dense_reference_key"],
        "fast_midpoint_count": int(freqs.size),
        "exact_reference_spot_count": int(selected.size),
        "top_interpolation_gap_spot_count": int(top.size),
        "exact_spot_fast_vs_prufer_max_abs_dex": float(np.max(fast_vs_exact)),
        "exact_spot_fast_vs_prufer_p95_abs_dex": float(np.quantile(fast_vs_exact, 0.95)),
        "exact_spot_prufer_vs_interpolated_dense_max_abs_dex": float(np.max(exact_vs_interp)),
        "exact_spot_prufer_vs_interpolated_dense_p95_abs_dex": float(np.quantile(exact_vs_interp, 0.95)),
        "previous_fast_vs_interpolated_dense_max_abs_dex": float(np.max(interpolated_gap[valid])),
        "previous_fast_vs_interpolated_dense_p95_abs_dex": float(np.quantile(interpolated_gap[valid], 0.95)),
        "runtime_s": elapsed,
        "spots": {
            "frequencies_log10_hz": freqs[selected].tolist(),
            "fast_omega_gw": fast_selected.tolist(),
            "exact_prufer_omega_gw": exact_ref.tolist(),
            "interpolated_dense_prufer_omega_gw": interp_selected.tolist(),
            "fast_vs_exact_abs_dex": fast_vs_exact.tolist(),
            "exact_vs_interpolated_dense_abs_dex": exact_vs_interp.tolist(),
        },
        "limitation": "Selected-frequency spot check only; it diagnoses midpoint interpolation disagreement but does not certify the full frequency integral or root.",
    }


def main() -> None:
    record = {
        "schema_version": 1,
        "experiment": "exact_prufer_spots_at_fast_midpoints",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "cases": [_case(ROOT / "docs" / name) for name in CASES],
        "resources": telemetry(threads=1),
        "runtime_environment": platform.platform(),
        "script_sha256": _sha(Path(__file__)),
    }
    OUTPUT.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": OUTPUT.name,
        "cases": [
            {key: case[key] for key in (
                "source_probe", "exact_reference_spot_count",
                "exact_spot_fast_vs_prufer_max_abs_dex",
                "exact_spot_prufer_vs_interpolated_dense_max_abs_dex",
                "runtime_s",
            )}
            for case in record["cases"]
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
