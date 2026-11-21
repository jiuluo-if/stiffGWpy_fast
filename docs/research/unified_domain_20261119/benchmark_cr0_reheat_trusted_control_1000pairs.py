"""Extended 1000-pair timing recheck for the near-threshold trusted CR0 control."""
from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "benchmark_cr1_upper_vs_fda_head_process_ab.py"
spec = importlib.util.spec_from_file_location("cr0_reheat_control_1000pairs", SOURCE)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load persistent full-Fast benchmark")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


def cases() -> dict[str, dict[str, float]]:
    return {
        "trusted_cr0_standard_control": {
            "Omega_bh2": 0.0223828, "Omega_ch2": 0.1201075,
            "H0": 67.32117, "DN_eff": 0.0, "A_s": 2.100549e-9,
            "r": 0.01, "n_t": 0.0, "cr": 0.0,
            "T_re": 2000.0, "DN_re": 10.0, "kappa10": 0.01,
        }
    }


if __name__ == "__main__":
    benchmark.BASE_SHA = "d4ee1b9c8f9fad401c5e22c2316094a608074727"
    benchmark.THREADING_LAYER = "tbb"
    benchmark.REPEATS = 1000
    benchmark.WARMUPS = 5
    benchmark.OUTPUT = HERE / "benchmark_cr0_reheat_trusted_control_1000pairs.json"
    benchmark._cases = cases
    benchmark.main()
