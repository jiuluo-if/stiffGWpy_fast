"""Test exact-background and smaller-step hypotheses on one fixed underflow grid."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_fast_background_step_falsifier_20261115.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_variant(params, freqs, frozen_dn, reference_log, *, h, sigma_exact):
    model = LCDM_SG(**params)
    model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid

    def fixed_grid(*args, **kwargs):
        return np.asarray(freqs, dtype=np.float64).copy()

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = float(frozen_dn)
        return original_gen(*args, **kwargs)

    FA.goal_oriented_freqs = fixed_grid
    FA.construct_fast_frequency_grid = fixed_grid
    model.gen_expansion = fixed_background_gen
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=0.99,
            z_tail=5.0,
            threads=1,
            accuracy_mode="fast",
            freq_grid="goal",
            physical_guard=False,
            h=h,
            sigma_exact=sigma_exact,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    elapsed = time.perf_counter() - started
    if result is None or not np.array_equal(model.f, freqs):
        raise RuntimeError(f"Fast falsifier failed: {model.fast_failure_reason}")
    log_spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
    errors = np.abs(log_spectrum - reference_log)
    return {
        "h": h,
        "sigma_exact": sigma_exact,
        "tail_z": 5.0,
        "time_node_count": int(model.Nv.size),
        "max_time_step": float(np.max(np.diff(model.Nv))),
        "min_time_step": float(np.min(np.diff(model.Nv))),
        "runtime_s": elapsed,
        "map_DN_gw": float(model.DN_gw[-1]),
        "map_residual_at_frozen_background": float(model.DN_gw[-1] - frozen_dn),
        "spectrum_max_abs_error_dex": float(np.max(errors)),
        "spectrum_p95_abs_error_dex": float(np.quantile(errors, 0.95)),
        "log10_OmegaGW": log_spectrum.tolist(),
    }


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    params = data["parameters"]
    freqs = np.asarray(data["frequency_grids"]["235"]["log10_hz"], dtype=np.float64)
    frozen_dn = float(data["background"]["common_DN_eff"])
    reference_log = np.asarray(
        data["reference_fixed_background"]["235"]["log10_OmegaGW"],
        dtype=np.float64,
    )
    variants = (
        (0.01, False),
        (0.01, True),
        (0.005, True),
        (0.0025, True),
        (0.00125, True),
        (0.000625, True),
    )
    results = [
        _run_variant(
            params,
            freqs,
            frozen_dn,
            reference_log,
            h=h,
            sigma_exact=sigma_exact,
        )
        for h, sigma_exact in variants
    ]
    output = {
        "schema_version": 1,
        "status": "complete_single_point_background_step_falsifier",
        "parameters": params,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": data["frequency_grids"]["235"]["sha256_f64le"],
        "frozen_DN_eff": frozen_dn,
        "reference_tail_z": 8.0,
        "fast_tail_z": 5.0,
        "variants": results,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive.py": _sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "exact_background.py": _sha(ROOT / "stiffgwpy_fast" / "exact_background.py"),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "stiff_SGWB.py": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "input_json": _sha(INPUT),
            "falsifier_script": _sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "variants": [
            {
                "h": row["h"],
                "sigma_exact": row["sigma_exact"],
                "map_residual": row["map_residual_at_frozen_background"],
                "spectrum_max_p95_dex": [
                    row["spectrum_max_abs_error_dex"],
                    row["spectrum_p95_abs_error_dex"],
                ],
            }
            for row in results
        ],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
