"""Paired CR0 Fast screen around the physical high-kappa boundary."""
from __future__ import annotations

import importlib.util
import itertools
import json
import math
from pathlib import Path

from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "cr0_kappa_tilt_extension_runner", HERE / "screen_h2_identity_region.py"
)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the pinned paired-screen runner")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)

LOG_K_RANGE = (1.0, 2.0)
NT_RANGE = (-0.5, -0.25)
LOG_K_GRID = (1.0, 1.5, 2.0)
NT_GRID = (-0.5, -0.375, -0.25)
screen.BASE_PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": -0.375,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 10.0,
}


def run(output_name: str, rectangle: list[dict[str, float]], holdouts: list[dict[str, float]]) -> None:
    screen.RECTANGLE = rectangle
    screen.HOLDOUTS = holdouts
    screen.OUTPUT = HERE / output_name
    screen.BASELINE_PARTIAL = screen.OUTPUT.with_name(screen.OUTPUT.stem + "_baseline_partial.json")
    screen.main()
    data = json.loads(screen.OUTPUT.read_text(encoding="utf-8"))
    data["params_contract"]["derived_tilt_rule"] = "cr=0: n_t remains independent"
    data["params_contract"]["prior_scope"] = {
        "log10T_re_GeV": [-1.0, 7.0], "n_t": [-1.0, 6.0],
        "log10kappa10": [-27.0, 3.0],
    }
    data["params_contract"]["region"] = {
        "log10kappa10": list(LOG_K_RANGE), "n_t": list(NT_RANGE),
        "kappa10": [10.0, 100.0], "sobol_measure": "uniform in log10(kappa10) and n_t",
        "physical_guard": "enabled; shared total N_eff <= 5",
    }
    screen.OUTPUT.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    lattice = [
        {"kappa10": 10.0**log_k, "n_t": tilt}
        for log_k, tilt in itertools.product(LOG_K_GRID, NT_GRID)
    ]
    holdouts = [
        {"kappa10": math.sqrt(10.0 * 100.0), "n_t": -0.46875},
        {"kappa10": 10.0**1.75, "n_t": -0.28125},
    ]
    run("screen_cr0_kappa_tilt_extension_lattice.json", lattice, holdouts)

    unit = qmc.Sobol(d=2, scramble=True, seed=20261211).random_base2(m=3)
    sobol = [
        {
            "kappa10": 10.0 ** (LOG_K_RANGE[0] + (LOG_K_RANGE[1] - LOG_K_RANGE[0]) * point[0]),
            "n_t": NT_RANGE[0] + (NT_RANGE[1] - NT_RANGE[0]) * point[1],
        }
        for point in unit
    ]
    run("screen_cr0_kappa_tilt_extension_sobol.json", [], sobol)


if __name__ == "__main__":
    main()
