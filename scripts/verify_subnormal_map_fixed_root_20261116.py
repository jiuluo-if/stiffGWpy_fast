"""Independently remap the tiny Fast root on its returned background."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude, construct_fast_frequency_grid  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BASE = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "n_t": 0.0,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 1e-2,
    "r": 1e-315,
}
OUTPUT = ROOT / "docs" / "research" / "results" / "subnormal_map_fixed_root_20261116.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    with (ROOT / "docs" / "research" / "results" / "subnormal_map_grid_convergence_20261115.json").open(
        encoding="utf-8"
    ) as handle:
        grid_results = json.load(handle)["cases"]
    rows = []
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    try:
        for row in grid_results:
            initial_root = float(row["fast"]["root_DN_eff"])
            bg_zero = LCDM_SG(**{**BASE, "DN_eff": 0.0})
            bg_zero.gen_expansion()
            model = LCDM_SG(**{**BASE, "DN_eff": initial_root})
            model.gen_expansion()
            bg_names = ("Nv", "sigma", "f_hor")
            background_equal = {
                name: bool(np.array_equal(getattr(model, name), getattr(bg_zero, name)))
                for name in bg_names
            }
            background_hash = {
                name: hashlib.sha256(
                    np.asarray(getattr(model, name), dtype="<f8").tobytes()
                ).hexdigest()
                for name in bg_names
            }
            log_amp = _log_tensor_amplitude(model.cosmo_param, model.derived_param)
            freqs = np.asarray(
                construct_fast_frequency_grid(
                    model, freq_res=float(row["freq_res"]), log_amplitude=log_amp
                ),
                dtype=np.float64,
            )
            FA.goal_oriented_freqs = lambda *args, _f=freqs, **kwargs: _f.copy()
            FA.construct_fast_frequency_grid = lambda *args, _f=freqs, **kwargs: _f.copy()
            grid_digest = hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest()
            result = model.SGWB_iter(
                engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
                threads=1, accuracy_mode="fast", freq_grid="goal",
                physical_guard=True,
            )
            mapped = None if result is None else float(model.DN_gw[-1])
            rows.append({
                "freq_res": row["freq_res"],
                "frequency_count": int(freqs.size),
                "frequency_sha256_f64le": grid_digest,
                "expected_frequency_sha256_f64le": row["frequency_sha256_f64le"],
                "input_root_DN_eff": initial_root,
                "fixed_root_fast_returned": result is model,
                "fixed_root_mapped_DN_gw": mapped,
                "relative_fixed_point_residual": (
                    None if mapped is None else mapped / initial_root - 1.0
                ),
                "root_converged_delta": getattr(model, "dn_converged_delta", None),
                "failure_reason": model.fast_failure_reason,
                "failure_detail": getattr(model, "fast_failure_detail", None),
                "neff_base_absorbs_root_exactly": bool(gp.Neff0 + initial_root == gp.Neff0),
                "background_arrays_bitwise_equal_to_dn0": background_equal,
                "background_array_sha256": background_hash,
            })
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    output = {
        "schema_version": 1,
        "status": "complete_direct_fast_fixed_root_remap",
        "scope": "For the cr=0,n_t=0,r=1e-315 root from freq_res=1/2, rerun Fast with DN_eff initialized to that root on the identical grid and compare the mapped root independently.",
        "parameters": BASE,
        "fast_options": {"tol": 1e-7, "z_tail": 5.0, "threads": 1},
        "cases": rows,
        "source_sha256": {
            "fast_sgwb.py": sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive.py": sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "script": sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "cases": [{
            "freq_res": x["freq_res"], "count": x["frequency_count"],
            "map": x["fixed_root_mapped_DN_gw"],
            "relative_residual": x["relative_fixed_point_residual"],
            "background_equal": x["background_arrays_bitwise_equal_to_dn0"],
            "grid_match": x["frequency_sha256_f64le"] == x["expected_frequency_sha256_f64le"],
        } for x in rows],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
