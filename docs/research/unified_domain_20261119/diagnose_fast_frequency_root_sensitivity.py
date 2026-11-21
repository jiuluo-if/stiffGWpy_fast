"""Compare Fast self-consistency roots across frequency-support densities."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import stiffgwpy_fast  # noqa: E402
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

OUTPUT = Path(__file__).with_suffix(".json")
if Path(stiffgwpy_fast.__file__).resolve().parent != ROOT / "stiffgwpy_fast":
    raise RuntimeError("diagnostic must import stiffgwpy_fast from the repository root")
SOURCE_NAMES = (
    "exact_background.py",
    "fast_sgwb.py",
    "freq_adaptive.py",
    "reference.py",
)
COMMON = {
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
CASES = {
    "trusted_default": {
        **COMMON,
        "r": 1.0e-2,
        "n_t": 0.0,
        "cr": 1.0,
        "T_re": 2.0e3,
    },
    "trusted_lowT": {
        **COMMON,
        "r": 1.0e-2,
        "n_t": 0.0,
        "cr": 1.0,
        "T_re": 1.0e1,
    },
    "trusted_highT": {
        **COMMON,
        "r": 1.0e-2,
        "n_t": 0.0,
        "cr": 1.0,
        "T_re": 1.0e4,
    },
    "trusted_stiff": {
        **COMMON,
        "r": 1.0e-1,
        "n_t": 0.0,
        "cr": 1.0,
        "T_re": 2.0e3,
    },
    "trusted_high_kappa": {
        **COMMON,
        "r": 1.0e-2,
        "n_t": 0.0,
        "cr": 1.0,
        "T_re": 2.0e3,
        "kappa10": 1.0,
    },
    "negative_tilt_low_reheat": {
        **COMMON,
        "Omega_bh2": 1.0e6,
        "r": 1.0e-2,
        "n_t": -0.5,
        "cr": 0.0,
        "T_re": 4.0e3,
    },
    "negative_tilt_high_reheat": {
        **COMMON,
        "Omega_bh2": 1.0e6,
        "r": 1.0e-2,
        "n_t": -0.25,
        "cr": 0.0,
        "T_re": 1.0e7,
    },
    "negative_tilt_high_reheat_corner": {
        **COMMON,
        "Omega_bh2": 1.0e6,
        "r": 1.0e-2,
        "n_t": -0.5,
        "cr": 0.0,
        "T_re": 1.0e7,
    },
    "negative_tilt_high_density": {
        **COMMON,
        "Omega_bh2": 13_890_136.71875,
        "r": 1.0e-2,
        "n_t": -0.022705078125,
        "cr": 0.0,
        "T_re": 2.0e3,
    },
    "positive_tilt_reheat": {
        **COMMON,
        "Omega_bh2": 1.0e6,
        "r": 1.0e-8,
        "n_t": 0.065,
        "cr": 0.0,
        "T_re": 5.0e3,
    },
}
GRID_SPECS = (
    ("goal", 1.0),
    ("grid_independent", 1.15),
    ("grid_independent", 2.3),
    ("grid_independent", 4.6),
)


def _solve(params: dict[str, float], grid: str, freq_res: float) -> dict:
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
        "seconds_single_run": time.perf_counter() - started,
        "quadrature_error_rel": getattr(
            model, "estimated_DN_quadrature_error_rel", None
        ),
    }


def main() -> None:
    rows = []
    for case_id, params in CASES.items():
        grid_results = []
        for grid, freq_res in GRID_SPECS:
            result = _solve(params, grid, freq_res)
            grid_results.append(
                {"grid": grid, "freq_res": freq_res, **result}
            )
        valid = [row for row in grid_results if row["success"]]
        finest = max(valid, key=lambda row: row["Nf"]) if valid else None
        if finest is not None:
            for row in grid_results:
                root = row["root"]
                row["relative_to_finest_fast"] = (
                    abs(root - finest["root"]) / abs(finest["root"])
                    if root is not None
                    else None
                )
        rows.append(
            {"case_id": case_id, "params": params, "grid_results": grid_results}
        )
        print(json.dumps(rows[-1], sort_keys=True), flush=True)

    payload = {
        "status": "complete_fast_frequency_support_sensitivity_screen",
        "reference_code_modified": False,
        "reference_maps_run": 0,
        "accuracy_claim": False,
        "grid_specs": [
            {"grid": grid, "freq_res": freq_res} for grid, freq_res in GRID_SPECS
        ],
        "source_sha256": {
            name: hashlib.sha256(
                (ROOT / "stiffgwpy_fast" / name).read_bytes()
            ).hexdigest()
            for name in SOURCE_NAMES
        },
        "cases": rows,
        "caveat": (
            "Fast-only grid sensitivity screen. The finest grid is not a "
            "Reference oracle or a continuum error certificate. Timings are "
            "single runs, not performance A/B evidence."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "output": str(OUTPUT)}))


if __name__ == "__main__":
    main()
