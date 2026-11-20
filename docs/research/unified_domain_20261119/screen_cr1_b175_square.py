from __future__ import annotations

import importlib.util
import itertools
from pathlib import Path

from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "screen_h2_identity_region.py"
spec = importlib.util.spec_from_file_location("cr1_region_screen", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load paired Fast screen runner")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)

LOW = 1.0e5
HIGH = 7.0e5
BASE = {"T_re": 2000.0, "cr": 1.0}


def run(output_name: str, rectangle, holdouts) -> None:
    screen.BASE_PARAMS.update(BASE)
    screen.RECTANGLE = rectangle
    screen.HOLDOUTS = holdouts
    screen.OUTPUT = HERE / output_name
    screen.BASELINE_PARTIAL = screen.OUTPUT.with_name(
        screen.OUTPUT.stem + "_baseline_partial.json"
    )
    screen.main()


def main() -> None:
    grid = (LOW, (LOW + HIGH) / 2.0, HIGH)
    rectangle = [
        {"Omega_bh2": ob, "Omega_ch2": oc}
        for ob, oc in itertools.product(grid, repeat=2)
    ]
    interior = [
        {"Omega_bh2": 2.5e5, "Omega_ch2": 4.0e5},
        {"Omega_bh2": 5.5e5, "Omega_ch2": 2.0e5},
    ]
    run("screen_cr1_b175_square_lattice.json", rectangle, interior)

    unit = qmc.Sobol(d=2, scramble=True, seed=20261120).random_base2(m=3)
    sobol = [
        {"Omega_bh2": float(LOW + (HIGH - LOW) * point[0]),
         "Omega_ch2": float(LOW + (HIGH - LOW) * point[1])}
        for point in unit
    ]
    run("screen_cr1_b175_square_sobol.json", [], sobol)


if __name__ == "__main__":
    main()
