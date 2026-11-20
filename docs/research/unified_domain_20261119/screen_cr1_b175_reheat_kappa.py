"""Paired CR1 screen over reheating temperature and kappa10."""
from __future__ import annotations

import importlib.util
import itertools
import json
from pathlib import Path

from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "screen_h2_identity_region.py"
spec = importlib.util.spec_from_file_location("cr1_reheat_kappa_paired_screen", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load pinned-baseline paired screen")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)

BASE = {
    "Omega_bh2": 7.0e5,
    "Omega_ch2": 7.0e5,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 1.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
T_GRID = (2000.0, 2250.0, 2500.0)
K_GRID = (0.01, 0.0125, 0.015)
SEED = 20261126


def run(name: str, rectangle, holdouts) -> None:
    screen.BASE_PARAMS = dict(BASE)
    screen.RECTANGLE = rectangle
    screen.HOLDOUTS = holdouts
    screen.OUTPUT = HERE / name
    screen.BASELINE_PARTIAL = screen.OUTPUT.with_name(
        screen.OUTPUT.stem + "_baseline_partial.json"
    )
    screen.main()
    data = json.loads(screen.OUTPUT.read_text(encoding="utf-8"))
    data["params_contract"]["derived_tilt_rule"] = "CR1 model-derived consistency tilt"
    data["params_contract"]["region"] = {
        "Omega_bh2": 7.0e5, "Omega_ch2": 7.0e5,
        "T_re": [2000.0, 2500.0], "kappa10": [0.01, 0.015],
    }
    screen.OUTPUT.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    rectangle = [
        {"T_re": t_re, "kappa10": kappa}
        for t_re, kappa in itertools.product(T_GRID, K_GRID)
    ]
    holdouts = [
        {"T_re": 2125.0, "kappa10": 0.01125},
        {"T_re": 2375.0, "kappa10": 0.01375},
    ]
    run("screen_cr1_b175_reheat_kappa_lattice.json", rectangle, holdouts)

    unit = qmc.Sobol(d=2, scramble=True, seed=SEED).random_base2(m=3)
    sobol = [
        {"T_re": float(2000.0 + 500.0 * point[0]),
         "kappa10": float(0.01 + 0.005 * point[1])}
        for point in unit
    ]
    run("screen_cr1_b175_reheat_kappa_sobol.json", [], sobol)


if __name__ == "__main__":
    main()
