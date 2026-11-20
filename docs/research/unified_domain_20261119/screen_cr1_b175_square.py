from __future__ import annotations

import argparse
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


def run(prefix: str, family: str, rectangle, holdouts) -> None:
    screen.BASE_PARAMS.update(BASE)
    screen.RECTANGLE = rectangle
    screen.HOLDOUTS = holdouts
    screen.OUTPUT = HERE / f"{prefix}_{family}.json"
    screen.BASELINE_PARTIAL = screen.OUTPUT.with_name(
        screen.OUTPUT.stem + "_baseline_partial.json"
    )
    screen.main()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--low", type=float, default=LOW)
    parser.add_argument("--high", type=float, default=HIGH)
    parser.add_argument("--seed", type=int, default=20261120)
    parser.add_argument("--prefix", default="screen_cr1_b175_square")
    args = parser.parse_args()
    if not (0.0 < args.low < args.high):
        parser.error("require 0 < --low < --high")
    grid = (args.low, (args.low + args.high) / 2.0, args.high)
    rectangle = [
        {"Omega_bh2": ob, "Omega_ch2": oc}
        for ob, oc in itertools.product(grid, repeat=2)
    ]
    span = args.high - args.low
    interior = [
        {"Omega_bh2": args.low + 0.25 * span,
         "Omega_ch2": args.low + 0.50 * span},
        {"Omega_bh2": args.low + 0.75 * span,
         "Omega_ch2": args.low + span / 6.0},
    ]
    run(args.prefix, "lattice", rectangle, interior)

    unit = qmc.Sobol(d=2, scramble=True, seed=args.seed).random_base2(m=3)
    sobol = [
        {"Omega_bh2": float(args.low + (args.high - args.low) * point[0]),
         "Omega_ch2": float(args.low + (args.high - args.low) * point[1])}
        for point in unit
    ]
    run(args.prefix, "sobol", [], sobol)


if __name__ == "__main__":
    main()
