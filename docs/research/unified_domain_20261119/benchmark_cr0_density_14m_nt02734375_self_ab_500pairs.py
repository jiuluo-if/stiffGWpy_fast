"""Full persistent-process A/B for the CR0 -0.02734375 tilt extension."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "benchmark_cr0_density_13m_to13_984375m_self_ab_500pairs.py"
spec = importlib.util.spec_from_file_location("cr0_density_14m_nt02734375_ab", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the established paired full-solver benchmark")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)
benchmark.REPEATS = 500
benchmark.OUTPUT = HERE / "benchmark_cr0_density_14m_nt02734375_self_ab_tbb_500pairs.json"


def cases() -> dict[str, dict[str, float]]:
    lattice = json.loads((HERE / "screen_cr0_density_14m_nt02734375_slice.json").read_text(encoding="utf-8"))
    sobol = json.loads((HERE / "screen_cr0_density_14m_nt02734375_slice_sobol.json").read_text(encoding="utf-8"))
    upper_low_tilt = next(
        row for row in lattice["cases"] if row["input"]["case_id"] == "rectangle_06"
    )["input"]["params"]
    interior = next(
        row for row in lattice["cases"] if row["input"]["case_id"] == "rectangle_04"
    )["input"]["params"]
    independent_holdout = next(
        row for row in sobol["cases"] if row["input"]["case_id"] == "holdout_05"
    )["input"]["params"]
    return {
        "new_region_upper_density_low_tilt": upper_low_tilt,
        "new_region_interior": interior,
        "new_region_independent_sobol": independent_holdout,
        "trusted_physical_cr0_control": {
            "Omega_bh2": 0.022, "Omega_ch2": 0.12, "H0": 67.32117,
            "DN_eff": 0.0, "A_s": 2.100549e-9, "r": 0.01,
            "n_t": 0.0, "cr": 0.0, "T_re": 2000.0,
            "DN_re": 10.0, "kappa10": 0.01,
        },
    }


benchmark._cases = cases


if __name__ == "__main__":
    benchmark.main()
