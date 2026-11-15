"""Attribute the worst 937-point minimum-subnormal error to its handoff state."""

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
sys.path.insert(0, str(ROOT / "scripts"))

import underflow_single_mode_state_attribution_20261115 as prior  # noqa: E402

from stiffgwpy_fast.global_param import f_piv  # noqa: E402
from stiffgwpy_fast.reference import solve_reference_mode  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

GRID_INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_grid_20261115.json"
TAIL_INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_tail_depth_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_handoff_20261115.json"
FREQUENCY = -16.291944809173845
ROOT_DN = 1.3543440396198841e-9


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    grid_data = json.loads(GRID_INPUT.read_text(encoding="utf-8"))
    tail_data = json.loads(TAIL_INPUT.read_text(encoding="utf-8"))
    params = grid_data["parameters"]
    grid937 = next(row for row in grid_data["cases"] if row["frequency_count"] == 937)
    tail937 = next(row for row in tail_data["cases"] if row["frequency_count"] == 937)
    freqs = np.asarray(grid937["frequency_log10_hz"], dtype=np.float64)
    index = int(np.argmin(np.abs(freqs - FREQUENCY)))
    if abs(float(freqs[index]) - FREQUENCY) > 1e-12:
        raise AssertionError("target outlier frequency not found on saved native grid")

    fast_states = {
        str(z): prior._fast_handoff_state(
            params, FREQUENCY, ROOT_DN, 0.005, False, z
        )
        for z in (5.0, 8.0)
    }
    unit_params = {**params, "r": 1.0, "A_s": 1.0, "n_t": 0.0, "DN_eff": 0.0}
    model = LCDM_SG(**unit_params)
    model.gen_expansion()
    reference_states = {
        str(z): solve_reference_mode(
            model, FREQUENCY, ROOT_DN, z_tail=z, rtol=1e-9, backend="cartesian"
        )
        for z in (5.0, 8.0, 10.0)
    }
    log_power = (
        math.log(params["r"])
        + math.log(params["A_s"])
        + params["n_t"] * (math.log(10.0) * FREQUENCY - math.log(f_piv))
    )
    reference_spectra = {}
    for key, state in reference_states.items():
        base = 3.0 * state["Opgw_today"] if state["used_tail"] else state["Ogw_today"] - state["Oj_today"]
        reference_spectra[key] = math.log10(base) + log_power / math.log(10.0)
    fast_spectra = {
        "5": float(grid937["fast"]["log10_OmegaGW"][index]),
        "8": float(tail937["fast_z8"]["log10_OmegaGW"][index]),
    }
    output = {
        "schema_version": 1,
        "status": "complete_min_subnormal_worst_mode_handoff_attribution",
        "scope": "Single 937-grid outlier mode at the minimum positive subnormal r; contrast Fast and Reference tail events under identical background and frequency.",
        "parameters": params,
        "frequency_index_937": index,
        "log10_frequency_hz": float(freqs[index]),
        "frozen_DN_eff": ROOT_DN,
        "fast_handoff_states": fast_states,
        "reference_states": {
            key: {name: state.get(name) for name in (
                "used_tail", "event_N", "eps_handoff", "n_steps", "matching_error_rel"
            )}
            for key, state in reference_states.items()
        },
        "point_spectra_log10_OmegaGW": {
            "Fast_z5": fast_spectra["5"],
            "Fast_z8": fast_spectra["8"],
            "Reference_z5": reference_spectra["5.0"],
            "Reference_z8": reference_spectra["8.0"],
            "Reference_z10": reference_spectra["10.0"],
            "Fast_z5_minus_Reference_z10_dex": fast_spectra["5"] - reference_spectra["10.0"],
            "Fast_z8_minus_Reference_z10_dex": fast_spectra["8"] - reference_spectra["10.0"],
            "Reference_z8_minus_z10_dex": reference_spectra["8.0"] - reference_spectra["10.0"],
        },
        "interpretation": "At Fast z_tail=8 the mode never reaches z=8 before today (last z=4.9907), so it stays in numerical propagation; at z_tail=5 it transfers to the terminal WKB approximation near epsilon=2.5e-3. Reference z8/z10 also do not use a tail for this mode.",
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            name: sha(ROOT / path)
            for name, path in {
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "single_mode_script": "scripts/underflow_single_mode_state_attribution_20261115.py",
                "grid_input": "docs/research/results/underflow_min_subnormal_reference_grid_20261115.json",
                "tail_input": "docs/research/results/underflow_min_subnormal_tail_depth_20261115.json",
                "diagnostic_script": "scripts/underflow_min_subnormal_handoff_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "frequency_index": index,
        "frequency_log10_hz": float(freqs[index]),
        "fast_z5_last_z": fast_states["5.0"]["last_node_before_tail_z"],
        "fast_z5_sigma": fast_states["5.0"]["grid_sigma_at_handoff"],
        "fast_z8_last_z": fast_states["8.0"]["last_node_before_tail_z"],
        "reference_tail_usage": {k: v["used_tail"] for k, v in output["reference_states"].items()},
        "point_errors_dex": output["point_spectra_log10_OmegaGW"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
