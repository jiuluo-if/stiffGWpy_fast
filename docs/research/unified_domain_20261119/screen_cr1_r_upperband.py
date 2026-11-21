"""Paired CR1 Fast screen immediately above the accepted lower-r strip."""
from __future__ import annotations

import importlib.util
import itertools
import json
from pathlib import Path

from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("cr1_r_upperband_runner", HERE / "screen_h2_identity_region.py")
if spec is None or spec.loader is None:
    raise RuntimeError("could not load pinned-baseline paired screen")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)

BASE = {
    "Omega_bh2": 7.0e5, "Omega_ch2": 7.0e5, "H0": 67.32117,
    "DN_eff": 0.0, "A_s": 2.100549e-9, "r": 0.0075,
    "n_t": 0.0, "cr": 1.0, "T_re": 1875.0,
    "DN_re": 10.0, "kappa10": 0.01,
}
R_RANGE = (0.005, 0.01)
T_RANGE = (1750.0, 2000.0)
R_GRID = (0.005, 0.0075, 0.01)
T_GRID = (1750.0, 1875.0, 2000.0)
SEED = 20261205


def run(name: str, rectangle: list[dict[str, float]], holdouts: list[dict[str, float]]) -> None:
    screen.BASE_PARAMS = dict(BASE)
    screen.RECTANGLE = rectangle
    screen.HOLDOUTS = holdouts
    screen.OUTPUT = HERE / name
    screen.BASELINE_PARTIAL = screen.OUTPUT.with_name(screen.OUTPUT.stem + "_baseline_partial.json")
    screen.main()
    data = json.loads(screen.OUTPUT.read_text(encoding="utf-8"))
    data["params_contract"]["derived_tilt_rule"] = "n_t=-r/8 for cr>0"
    data["params_contract"]["region"] = {
        "Omega_bh2": 7.0e5, "Omega_ch2": 7.0e5,
        "T_re": list(T_RANGE), "r": list(R_RANGE), "kappa10": 0.01,
        "sobol_measure": "uniform in T_re and r",
    }
    screen.OUTPUT.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    rectangle = [{"T_re": t, "r": ratio} for t, ratio in itertools.product(T_GRID, R_GRID)]
    holdouts = [
        {"T_re": 1781.25, "r": 0.00625},
        {"T_re": 1968.75, "r": 0.00875},
    ]
    run("screen_cr1_r_upperband_lattice.json", rectangle, holdouts)
    unit = qmc.Sobol(d=2, scramble=True, seed=SEED).random_base2(m=3)
    sobol = [
        {"T_re": float(T_RANGE[0]+(T_RANGE[1]-T_RANGE[0])*point[0]),
         "r": float(R_RANGE[0]+(R_RANGE[1]-R_RANGE[0])*point[1])}
        for point in unit
    ]
    run("screen_cr1_r_upperband_sobol.json", [], sobol)


if __name__ == "__main__":
    main()
