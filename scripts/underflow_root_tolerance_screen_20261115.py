"""Measure root tolerance sensitivity and independently remap returned roots."""

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

from stiffgwpy_fast import fast_sgwb as fast_module  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BASE = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "n_t": 17.38,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 1e-2,
}
AMPLITUDES = (1e-316, 1e-318, 5e-324)
TOLERANCES = (1e-5, 1e-7, 1e-9, 1e-11)
GRID_RESULT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_root_tolerance_screen_20261115.json"


def solve(r: float, freqs: np.ndarray, *, tol: float, fixed_dn: float | None = None) -> dict[str, object]:
    model = LCDM_SG(**{**BASE, "r": r})
    if fixed_dn is not None:
        model.cosmo_param["DN_eff"] = fixed_dn
        model.gen_expansion()
    original_gen = model.gen_expansion
    original_gen_fast = fast_module.gen_fast
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    outer_background_evaluations = 0

    def count_gen_fast(*args, **kwargs):
        nonlocal outer_background_evaluations
        outer_background_evaluations += 1
        return original_gen_fast(*args, **kwargs)

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = fixed_dn
        return original_gen(*args, **kwargs)

    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    fast_module.gen_fast = count_gen_fast
    if fixed_dn is not None:
        model.gen_expansion = fixed_background_gen
    started = time.perf_counter()
    error = None
    try:
        result = model.SGWB_iter(
            engine="fast", fallback=False,
            tol=0.99 if fixed_dn is not None else tol,
            z_tail=5.0, threads=1, accuracy_mode="fast",
            freq_grid="goal", physical_guard=False,
        )
    except Exception as exc:
        result = None
        error = f"{type(exc).__name__}: {exc}"
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
        fast_module.gen_fast = original_gen_fast
    elapsed = time.perf_counter() - started
    if result is None:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}; {error}")
    return {
        "root_DN_eff": None if fixed_dn is not None else float(model.cosmo_param["DN_eff"]),
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "evaluation_background_DN_eff": (
            float(fixed_dn) if fixed_dn is not None else float(model.cosmo_param["DN_eff"])
        ),
        "runtime_s": elapsed,
        "outer_background_evaluations": outer_background_evaluations,
        "frequency_count": int(model.f.size),
        "frequency_sha256_f64le": hashlib.sha256(np.asarray(model.f, dtype="<f8").tobytes()).hexdigest(),
        "fast_failure_reason": model.fast_failure_reason,
    }


def main() -> int:
    saved = json.loads(GRID_RESULT.read_text(encoding="utf-8"))
    freqs = np.asarray(saved["frequency_grids"]["235"]["log10_hz"], dtype=np.float64)
    cases = []
    for r in AMPLITUDES:
        rows = []
        for tol in TOLERANCES:
            root_result = solve(r, freqs, tol=tol)
            root = float(root_result["root_DN_eff"])
            independent_map = solve(r, freqs, tol=tol, fixed_dn=root)
            independent_map["relative_fixed_point_residual"] = (
                independent_map["mapped_DN_gw"] / root - 1.0
            )
            rows.append({
                "outer_tol": tol,
                "self_consistent": root_result,
                "independent_map_at_returned_root": independent_map,
            })
        base_root = float(rows[0]["self_consistent"]["root_DN_eff"])
        for row in rows:
            row["relative_root_change_vs_tol_1e-5"] = (
                float(row["self_consistent"]["root_DN_eff"]) / base_root - 1.0
            )
        cases.append({"r": r, "tolerances": rows})

    output = {
        "schema_version": 1,
        "status": "complete_independent_root_tolerance_screen",
        "scope": "Three cr=0 underflow amplitudes on the same direct 235-node grid; each reported self-consistent root is independently remapped with background frozen at that root.",
        "parameters_except_r": BASE,
        "amplitudes": list(AMPLITUDES),
        "outer_tolerances": list(TOLERANCES),
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
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
                "diagnostic_script": "scripts/underflow_root_tolerance_screen_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "cases": [{
            "r": case["r"],
            "tolerances": [{
                "tol": row["outer_tol"],
                "root": row["self_consistent"]["root_DN_eff"],
                "root_change": row["relative_root_change_vs_tol_1e-5"],
                "independent_residual": row["independent_map_at_returned_root"]["relative_fixed_point_residual"],
            } for row in case["tolerances"]],
        } for case in cases],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
