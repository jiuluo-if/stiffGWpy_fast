"""CR1 reheating/kappa full-solve A/B using the common helper benchmark."""
from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "benchmark_cr1_b175_wkb_guard_ab.py"
spec = importlib.util.spec_from_file_location("cr1_reheat_kappa_wkb_benchmark", SOURCE)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load shared WKB benchmark")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)
benchmark.PARAMS.update({
    "cr": 1.0, "Omega_bh2": 7.0e5, "Omega_ch2": 7.0e5,
    "T_re": 2000.0, "kappa10": 0.01,
})
benchmark.REPEATS = 100
benchmark.CASES = (
    ("cr1_reheat_kappa_high_corner", {
        "T_re": 2500.0, "kappa10": 0.015,
    }),
    ("cr1_reheat_kappa_sobol_holdout", {
        "T_re": 2372.528196312487, "kappa10": 0.014170661885291338,
    }),
)
benchmark.OUTPUT = HERE / "benchmark_cr1_b175_reheat_kappa_wkb_ab_100pairs.json"


if __name__ == "__main__":
    benchmark.main()
