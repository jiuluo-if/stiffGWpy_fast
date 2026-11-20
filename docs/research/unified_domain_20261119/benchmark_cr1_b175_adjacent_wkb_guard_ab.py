"""CR1 adjacent-domain WKB-helper A/B wrapper using the shared benchmark."""
from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "benchmark_cr1_b175_wkb_guard_ab.py"
spec = importlib.util.spec_from_file_location("cr1_adjacent_wkb_benchmark", SOURCE)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load shared WKB benchmark")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)
benchmark.PARAMS.update({"cr": 1.0, "Omega_ch2": 8.0e5})
benchmark.REPEATS = 100
benchmark.CASES = (
    ("adjacent_square_high_corner", {
        "Omega_bh2": 8.0e5, "Omega_ch2": 8.0e5, "n_t": 0.0,
    }),
    ("adjacent_sobol_upper_sum", {
        "Omega_bh2": 784905.9135653079,
        "Omega_ch2": 799701.8547728658,
        "n_t": 0.0,
    }),
)
benchmark.OUTPUT = HERE / "benchmark_cr1_b175_adjacent_wkb_guard_ab_100pairs.json"


if __name__ == "__main__":
    benchmark.main()
