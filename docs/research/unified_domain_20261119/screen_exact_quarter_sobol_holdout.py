"""Paired Fast screen on an independent 2D Sobol sample in the CR1 patch."""

from __future__ import annotations

import screen_exact_quarter_expanded_patch as screen
from scipy.stats import qmc

SOBOL_SEED = 20261120
SOBOL_POWER = 4


def main() -> None:
    unit = qmc.Sobol(d=2, scramble=True, seed=SOBOL_SEED).random_base2(SOBOL_POWER)
    screen.COORDS = tuple(
        (650_000.0 + 50_000.0 * float(u),
         650_000.0 + 50_000.0 * float(v))
        for u, v in unit
    )
    screen.OUTPUT = screen.HERE / "screen_exact_quarter_sobol_holdout.json"
    screen.main()


if __name__ == "__main__":
    main()
