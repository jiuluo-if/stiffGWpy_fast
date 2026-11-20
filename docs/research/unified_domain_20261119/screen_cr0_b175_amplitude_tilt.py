"""Current-source paired CR0 screen over positive tilt and tensor amplitude."""
from __future__ import annotations

import importlib.util
import itertools
import json
from pathlib import Path

from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "screen_h2_identity_region.py"
spec = importlib.util.spec_from_file_location("cr0_h2_identity_screen", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load pinned-baseline paired screen")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)

BASE = {
    "Omega_bh2": 1.0e6,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 1.0e-10,
    "n_t": 0.05,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
R_GRID = (1.0e-12, 1.0e-10, 1.0e-8)
NT_GRID = (0.05, 0.075, 0.1)
SEED = 20261123


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
    data["params_contract"]["derived_tilt_rule"] = "independent input tilt for cr=0"
    data["params_contract"]["region"] = {
        "r": [1.0e-12, 1.0e-8], "n_t": [0.05, 0.1],
        "Omega_bh2": 1.0e6, "Omega_ch2": 0.1201075,
    }
    screen.OUTPUT.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    rectangle = [
        {"r": r, "n_t": n_t}
        for r, n_t in itertools.product(R_GRID, NT_GRID)
    ]
    holdouts = [
        {"r": 1.0e-11, "n_t": 0.0875},
        {"r": 3.0e-10, "n_t": 0.0625},
    ]
    run("screen_cr0_b175_amplitude_tilt_lattice.json", rectangle, holdouts)

    unit = qmc.Sobol(d=2, scramble=True, seed=SEED).random_base2(m=3)
    sobol = [
        {"r": 10.0 ** (-12.0 + 4.0 * float(point[0])),
         "n_t": 0.05 + 0.05 * float(point[1])}
        for point in unit
    ]
    run("screen_cr0_b175_amplitude_tilt_sobol.json", [], sobol)


if __name__ == "__main__":
    main()
