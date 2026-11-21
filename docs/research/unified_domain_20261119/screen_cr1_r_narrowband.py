"""Paired screen for a data-selected narrow CR1 r strip above r=0.01."""
from __future__ import annotations

import importlib.util
import itertools
from pathlib import Path

from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "cr1_r_narrowband_runner", HERE / "screen_cr1_r_upperband_candidate.py"
)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the pinned paired-screen runner")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)

T_RANGE = (1875.0, 2250.0)
R_RANGE = (0.01, 0.01125)
T_GRID = (1875.0, 2062.5, 2250.0)
R_GRID = (0.01, 0.010625, 0.01125)
runner.BASE = {
    "Omega_bh2": 7.0e5,
    "Omega_ch2": 7.0e5,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 1.0,
    "T_re": 1937.5,
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
        {"T_re": 1968.75, "r": 0.011125},
        {"T_re": 2156.25, "r": 0.010375},
    ]
    runner.run("screen_cr1_r_narrowband_lattice.json", lattice, holdouts)

    points = qmc.Sobol(d=2, scramble=True, seed=20261208).random_base2(m=3)
    sobol = [
        {
            "T_re": float(T_RANGE[0] + (T_RANGE[1] - T_RANGE[0]) * point[0]),
            "r": float(R_RANGE[0] + (R_RANGE[1] - R_RANGE[0]) * point[1]),
        }
        for point in points
    ]
    runner.run("screen_cr1_r_narrowband_sobol.json", [], sobol)


if __name__ == "__main__":
    main()
