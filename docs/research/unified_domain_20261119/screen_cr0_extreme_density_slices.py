"""Paired current-HEAD/candidate Fast screens for two extreme-density slices."""
from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path

from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[3]
SHARED = Path(__file__).with_name("screen_h2_identity_region.py")
spec = importlib.util.spec_from_file_location("extreme_density_shared_screen", SHARED)
if spec is None or spec.loader is None:
    raise RuntimeError("cannot load paired-screen runner")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)

screen.BASE_SHA = subprocess.check_output(
    ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
).strip()
screen.BASELINE_TIMEOUT_SECONDS = 45.0
screen.BASE_PARAMS = {
    "Omega_bh2": 1.0e16,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 1.0e-16,
    "n_t": 0.0,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}


def _make_cases(slice_name: str) -> list[dict[str, object]]:
    if slice_name == "density_h0":
        rectangle = [
            {"Omega_bh2": density, "H0": h0}
            for density in (1.0e16, 1.5e16, 2.0e16)
            for h0 in (60.0, 67.5, 75.0)
        ]
        holdouts = [{"Omega_bh2": 1.25e16, "H0": 63.75},
                    {"Omega_bh2": 1.75e16, "H0": 71.25}]
        seed, dims = 20261123, ("Omega_bh2", "H0")
        bounds = ((1.0e16, 2.0e16), (60.0, 75.0))
    elif slice_name == "density_reheat":
        rectangle = [
            {"Omega_bh2": density, "T_re": tre}
            for density in (1.0e16, 1.5e16, 2.0e16)
            for tre in (1500.0, 2000.0, 2500.0)
        ]
        holdouts = [{"Omega_bh2": 1.25e16, "T_re": 1750.0},
                    {"Omega_bh2": 1.75e16, "T_re": 2250.0}]
        seed, dims = 20261124, ("Omega_bh2", "T_re")
        bounds = ((1.0e16, 2.0e16), (1500.0, 2500.0))
    else:
        raise ValueError(slice_name)

    cases = []
    for family, rows in (("rectangle", rectangle), ("holdout", holdouts)):
        for i, override in enumerate(rows):
            params = dict(screen.BASE_PARAMS)
            params.update(override)
            cases.append({"case_id": f"{family}_{i:02d}", "family": family,
                          "params": params})
    sampler = qmc.Sobol(d=2, scramble=True, seed=seed)
    points = sampler.random_base2(m=3)
    for i, point in enumerate(points):
        overrides = {
            dims[0]: bounds[0][0] + point[0] * (bounds[0][1] - bounds[0][0]),
            dims[1]: bounds[1][0] + point[1] * (bounds[1][1] - bounds[1][0]),
        }
        params = dict(screen.BASE_PARAMS)
        params.update(overrides)
        cases.append({"case_id": f"sobol_{i:02d}", "family": "sobol",
                      "params": params})
    return cases


def run(slice_name: str) -> None:
    screen._cases = lambda: _make_cases(slice_name)
    screen.OUTPUT = Path(__file__).with_name(f"screen_cr0_extreme_{slice_name}.json")
    screen.BASELINE_PARTIAL = screen.OUTPUT.with_name(
        f"screen_cr0_extreme_{slice_name}_baseline_partial.json"
    )
    screen.main()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("slice", choices=("density_h0", "density_reheat"))
    run(parser.parse_args().slice)
