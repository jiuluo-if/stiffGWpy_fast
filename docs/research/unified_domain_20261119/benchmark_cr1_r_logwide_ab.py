"""500-pair persistent-process A/B across the CR1 log-wide r rectangle."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "benchmark_cr1_upper_vs_fda_head_process_ab.py"
spec = importlib.util.spec_from_file_location("cr1_r_logwide_persistent_ab", SOURCE)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load persistent full-Fast benchmark")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)
LATTICE = json.loads((HERE / "screen_cr1_r_logwide_lattice.json").read_text(encoding="utf-8"))
SOBOL = json.loads((HERE / "screen_cr1_r_logwide_sobol.json").read_text(encoding="utf-8"))


def cases() -> dict[str, dict[str, float]]:
    lattice = {row["input"]["case_id"]: row["input"]["params"] for row in LATTICE["cases"]}
    sobol = {row["input"]["case_id"]: row["input"]["params"] for row in SOBOL["cases"]}
    return {
        "new_region_low_r_cutoff_corner": lattice["rectangle_06"],
        "new_region_high_r_corner": lattice["rectangle_02"],
        "new_region_independent_sobol": sobol["holdout_02"],
        "old_fast_success_mid_log_r": lattice["rectangle_07"],
        "trusted_physical_cr0_control": {
            "Omega_bh2": 0.022, "Omega_ch2": 0.12, "H0": 67.32117,
            "DN_eff": 0.0, "A_s": 2.100549e-9, "r": 0.01, "n_t": 0.0,
            "cr": 0.0, "T_re": 2000.0, "DN_re": 10.0, "kappa10": 0.01,
        },
    }


if __name__ == "__main__":
    benchmark.BASE_SHA = "d4ee1b9c8f9fad401c5e22c2316094a608074727"
    benchmark.REPEATS = 500
    benchmark.WARMUPS = 3
    benchmark.OUTPUT = HERE / "benchmark_cr1_r_logwide_ab_500pairs.json"
    benchmark._cases = cases
    benchmark.main()
