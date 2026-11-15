"""Test the z_today handoff mask on a physical-amplitude cr>0 point."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.reference import _Background, f_hor_abs  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_tail_depth_consistency_branch_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_consistency_branch_20261115.json"
THRESHOLDS = (5.0, 6.0, 7.0, 8.0)


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    params = data["parameters"]
    root = float(data["fast_root_DN_eff"])
    fast5 = np.asarray(data["fast_fixed_root"]["5"]["log10_OmegaGW"], dtype=np.float64)
    fast8 = np.asarray(data["fast_fixed_root"]["8"]["log10_OmegaGW"], dtype=np.float64)
    ref10 = np.asarray(data["reference_fixed_root_physical_amplitude"]["10"]["log10_OmegaGW"], dtype=np.float64)

    # Reconstruct the exact native frequency grid from the production root solve.
    model = LCDM_SG(**params)
    model.gen_expansion()
    model.SGWB_iter(
        engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
        threads=1, accuracy_mode="fast", freq_grid="goal", physical_guard=True,
    )
    freqs = np.asarray(model.f, dtype=np.float64)
    grid_hash = hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest()
    if freqs.size != data["frequency_count"] or grid_hash != data["frequency_grid_sha256_f64le"]:
        raise AssertionError("reconstructed consistency frequency grid differs from saved comparison")
    model = LCDM_SG(**params)
    model.gen_expansion()
    bg = _Background(model, root)
    f_hor_today = f_hor_abs(bg, bg.N_inf)
    z_today = (freqs - f_hor_today) * gp.ln10
    h = float(model.derived_param["h"])
    map_factor = gp.Neff0 * gp.ln10 / (gp.Omega_nh2 / (h * h))
    reference_map = float(data["reference_fixed_root_physical_amplitude"]["10"]["mapped_DN_gw"])

    candidates = []
    for threshold in THRESHOLDS:
        use_z8 = z_today < threshold
        hybrid = np.where(use_z8, fast8, fast5)
        error = hybrid - ref10
        ln_spectrum = hybrid * gp.ln10
        shift = float(np.max(ln_spectrum))
        scaled = np.exp(ln_spectrum - shift)
        integral = float(np.sum(FS._pchip_integrals_vectorized(freqs, scaled)))
        mapped = map_factor * math.exp(shift) * integral
        candidates.append({
            "z_today_threshold": threshold,
            "z8_mode_count": int(np.count_nonzero(use_z8)),
            "z8_mode_fraction": float(np.mean(use_z8)),
            "max_abs_spectrum_error_vs_reference_z10_dex": float(np.max(np.abs(error))),
            "p95_abs_spectrum_error_vs_reference_z10_dex": float(np.quantile(np.abs(error), 0.95)),
            "max_error_log10_hz": float(freqs[np.argmax(np.abs(error))]),
            "mapped_DN_gw": float(mapped),
            "relative_map_error_vs_reference_z10": float(mapped / reference_map - 1.0),
            "relative_fixed_root_residual": float(mapped / root - 1.0),
        })
    output = {
        "schema_version": 1,
        "status": "complete_cr_positive_adaptive_handoff_screen_not_production",
        "scope": "Physical-amplitude cr>0 consistency branch; reuse direct same-root z5/z8 Fast and z10 Reference spectra from the recorded tail-depth run, reconstruct exact 76-node grid.",
        "parameters": params,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": grid_hash,
        "frozen_DN_eff": root,
        "f_hor_today_log10_hz": f_hor_today,
        "z_today_range": [float(np.min(z_today)), float(np.max(z_today))],
        "reference_z10_map_DN_gw": reference_map,
        "reference_z10_map_residual": float(reference_map / root - 1.0),
        "threshold_candidates": candidates,
        "performance": "not measured; hybrid spectra are composed from separate full direct z5 and z8 runs",
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            name: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            for name, path in {
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "input_json": "docs/research/results/underflow_tail_depth_consistency_branch_20261115.json",
                "diagnostic_script": "scripts/screen_adaptive_handoff_consistency_branch_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "candidates": [{
            "threshold": row["z_today_threshold"],
            "z8_fraction": row["z8_mode_fraction"],
            "max_error": row["max_abs_spectrum_error_vs_reference_z10_dex"],
            "p95_error": row["p95_abs_spectrum_error_vs_reference_z10_dex"],
            "map_relative_error": row["relative_map_error_vs_reference_z10"],
            "root_residual": row["relative_fixed_root_residual"],
        } for row in candidates],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
