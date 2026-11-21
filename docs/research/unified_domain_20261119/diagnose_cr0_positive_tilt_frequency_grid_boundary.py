"""Measure CR0 positive-tilt self-consistency sensitivity to frequency support."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = Path(__file__).with_suffix(".json")
BASE_PARAMS = {
    "Omega_bh2": 1.0e6,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 1.0e-8,
    "DN_re": 10.0,
    "kappa10": 0.01,
    "cr": 0.0,
}
SOURCE_NAMES = (
    "exact_background.py",
    "fast_sgwb.py",
    "freq_adaptive.py",
    "reference.py",
)


def _run_fast(params: dict[str, float], grid: str, freq_res: float) -> dict:
    model = LCDM_SG(**params)
    config = resolve_config("fast", threads=2, freq_grid=grid)
    started = time.perf_counter()
    result = SGWB_iter_fast(
        model,
        config=config,
        tol=1.0e-6,
        physical_guard=True,
        freq_res=freq_res,
    )
    return {
        "success": result is model,
        "failure": getattr(model, "fast_failure_reason", None),
        "root": float(model.DN_gw[-1]) if result is model else None,
        "Nf": len(getattr(model, "f", [])),
        "seconds": time.perf_counter() - started,
        "quadrature_error_rel": getattr(
            model, "estimated_DN_quadrature_error_rel", None
        ),
    }


def main() -> None:
    cases = []
    for tre in (2000.0, 3500.0, 5000.0):
        for nt in (0.05, 0.065, 0.08):
            params = {**BASE_PARAMS, "T_re": tre, "n_t": nt}
            row = {"T_re": tre, "n_t": nt}
            row["goal"] = _run_fast(params, "goal", 1.0)
            row["grid_independent"] = _run_fast(params, "grid_independent", 2.0)
            goal_root = row["goal"]["root"]
            dense_root = row["grid_independent"]["root"]
            if goal_root is not None and dense_root is not None:
                row["root_relative_grid_shift"] = abs(goal_root - dense_root) / abs(
                    dense_root
                )
            cases.append(row)

    ladder = []
    for tre in (2000.0, 5000.0):
        params = {**BASE_PARAMS, "T_re": tre, "n_t": 0.065}
        for freq_res in (1.15, 2.3, 4.6):
            row = {"T_re": tre, "n_t": 0.065, "freq_res": freq_res}
            row.update(_run_fast(params, "grid_independent", freq_res))
            ladder.append(row)

    payload = {
        "status": "complete_fast_frequency_grid_sensitivity_screen",
        "reference_code_modified": False,
        "reference_maps_run": 0,
        "grids": {
            "production": "goal freq_res=1",
            "denser": "grid_independent freq_res=2",
        },
        "source_sha256": {
            name: hashlib.sha256(
                (ROOT / "stiffgwpy_fast" / name).read_bytes()
            ).hexdigest()
            for name in SOURCE_NAMES
        },
        "cases": cases,
        "frequency_ladder_grid_independent": ladder,
        "frequency_ladder_caveat": (
            "Fast-only; same-grid Reference was run only at 204 nodes for "
            "T_re=5000, n_t=0.065; no convergence certification is implied."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
