"""Compare z_tail=5/8 at the minimum-subnormal root on direct grids."""

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

INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_grid_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_tail_depth_20261115.json"


def run_fast(freqs: np.ndarray, params: dict[str, float], frozen_dn: float) -> dict[str, object]:
    model = LCDM_SG(**params)
    model.cosmo_param["DN_eff"] = frozen_dn
    model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = frozen_dn
        return original_gen(*args, **kwargs)

    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    model.gen_expansion = fixed_background_gen
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast", fallback=False, tol=0.99, z_tail=8.0,
            threads=1, accuracy_mode="fast", freq_grid="goal",
            physical_guard=False,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    elapsed = time.perf_counter() - started
    if result is None:
        raise RuntimeError(f"Fast z8 rejected point: {model.fast_failure_reason}")
    return {
        "runtime_s": elapsed,
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "map_residual_relative_to_frozen_root": float(model.DN_gw[-1] / frozen_dn - 1.0),
        "frequency_count": int(model.f.size),
        "log10_OmegaGW": np.asarray(model.log10OmegaGW, dtype=np.float64),
    }


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    params = data["parameters"]
    cases = []
    for baseline in data["cases"]:
        freqs = np.asarray(baseline["frequency_log10_hz"], dtype=np.float64)
        root = float(baseline["frozen_DN_eff"])
        fast = run_fast(freqs, params, root)
        ref = np.asarray(baseline["reference"]["log10_OmegaGW"], dtype=np.float64)
        delta = fast["log10_OmegaGW"] - ref
        abs_delta = np.abs(delta)
        cases.append({
            "frequency_count": int(freqs.size),
            "frequency_sha256_f64le": baseline["frequency_sha256_f64le"],
            "frozen_DN_eff": root,
            "z_tail": 8.0,
            "fast_z8": {
                "runtime_s": fast["runtime_s"],
                "mapped_DN_gw": fast["mapped_DN_gw"],
                "map_residual_relative_to_frozen_root": fast["map_residual_relative_to_frozen_root"],
                "log10_OmegaGW": fast["log10_OmegaGW"].tolist(),
            },
            "reference_z8": {
                "runtime_s": baseline["reference"]["runtime_s"],
                "mapped_DN_gw": baseline["reference"]["mapped_DN_gw"],
                "runtime_warnings": baseline["reference"]["runtime_warnings"],
            },
            "spectrum_error": {
                "max_abs_dex": float(np.max(abs_delta)),
                "p95_abs_dex": float(np.quantile(abs_delta, 0.95)),
                "max_error_index": int(np.argmax(abs_delta)),
                "max_error_log10_hz": float(freqs[np.argmax(abs_delta)]),
                "max_error_signed_dex": float(delta[np.argmax(abs_delta)]),
                "signed_error_dex": delta.tolist(),
            },
        })
    output = {
        "schema_version": 1,
        "status": "complete_min_subnormal_fast_tail_depth8_direct_reference_comparison",
        "scope": "Same minimum-subnormal root, native 469/937 grids and Reference samples as z5; only Fast z_tail changes from 5 to 8.",
        "parameters": params,
        "cases": cases,
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
                "freq_adaptive.py": "stiffgwpy_fast/freq_adaptive.py",
                "input_result": "docs/research/results/underflow_min_subnormal_reference_grid_20261115.json",
                "diagnostic_script": "scripts/underflow_min_subnormal_tail_depth_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "cases": [{
            "count": row["frequency_count"],
            "fast_z8_seconds": row["fast_z8"]["runtime_s"],
            "fast_z8_map_residual": row["fast_z8"]["map_residual_relative_to_frozen_root"],
            "spectrum_max_p95_dex": [row["spectrum_error"]["max_abs_dex"], row["spectrum_error"]["p95_abs_dex"]],
            "max_error_log10_hz": row["spectrum_error"]["max_error_log10_hz"],
        } for row in cases],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
