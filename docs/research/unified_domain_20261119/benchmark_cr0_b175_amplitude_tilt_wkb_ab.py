"""CR0 positive-tilt/amplitude full-solve A/B wrapper."""
from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "benchmark_cr1_b175_wkb_guard_ab.py"
spec = importlib.util.spec_from_file_location("cr0_amplitude_tilt_wkb_benchmark", SOURCE)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load shared WKB benchmark")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)
benchmark.PARAMS.update({"cr": 0.0, "Omega_ch2": 0.1201075})
benchmark.REPEATS = 100
benchmark.CASES = (
    ("cr0_high_amplitude_high_tilt", {
        "Omega_bh2": 1.0e6, "r": 1.0e-8, "n_t": 0.1,
    }),
    ("cr0_independent_sobol_holdout", {
        "Omega_bh2": 1.0e6, "r": 6.300178881203047e-9,
        "n_t": 0.09677163888700308,
    }),
)
benchmark.OUTPUT = HERE / "benchmark_cr0_b175_amplitude_tilt_wkb_ab_100pairs.json"


if __name__ == "__main__":
    benchmark.main()
