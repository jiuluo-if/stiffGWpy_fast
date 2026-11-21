"""Paired screen for a temperature-axis extension of the accepted CR1 r band."""
from __future__ import annotations

import importlib.util
import itertools
from pathlib import Path

from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "cr1_r_temperature_extension_runner", HERE / "screen_cr1_r_upperband_candidate.py"
)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the pinned paired-screen runner")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)

T_RANGE = (2000.0, 2250.0)
R_RANGE = (0.005, 0.01)
T_GRID = (2000.0, 2125.0, 2250.0)
R_GRID = (0.005, 0.0075, 0.01)
runner.BASE = {
    "Omega_bh2": 7.0e5,
    "Omega_ch2": 7.0e5,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.0075,
    "n_t": 0.0,
    "cr": 1.0,
    "T_re": 2125.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
runner.T_RANGE = T_RANGE
runner.R_RANGE = R_RANGE


def main() -> None:
    lattice = [
        {"T_re": t, "r": ratio}
        for t, ratio in itertools.product(T_GRID, R_GRID)
    ]
    holdouts = [
        {"T_re": 2062.5, "r": 0.00625},
        {"T_re": 2187.5, "r": 0.00875},
    ]
    runner.run("screen_cr1_r_temperature_extension_lattice.json", lattice, holdouts)

    points = qmc.Sobol(d=2, scramble=True, seed=20261209).random_base2(m=3)
    sobol = [
        {
            "T_re": float(T_RANGE[0] + (T_RANGE[1] - T_RANGE[0]) * point[0]),
            "r": float(R_RANGE[0] + (R_RANGE[1] - R_RANGE[0]) * point[1]),
        }
        for point in points
    ]
    runner.run("screen_cr1_r_temperature_extension_sobol.json", [], sobol)


if __name__ == "__main__":
    main()
