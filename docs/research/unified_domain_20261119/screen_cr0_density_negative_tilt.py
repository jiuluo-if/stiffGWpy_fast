"""Paired baseline/current Fast screen on the CR0 density x tilt lattice."""

from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "screen_h2_identity_region.py"
spec = importlib.util.spec_from_file_location("screen_h2_region_cr0_tilt", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the shared Fast comparison runner")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)

screen.BASE_PARAMS.update({
    "Omega_bh2": 1.0e5,
    "Omega_ch2": 0.1201075,
    "cr": 0.0,
    "r": 0.01,
    "T_re": 2000.0,
    "kappa10": 0.01,
})
screen.RECTANGLE = [
    {"Omega_bh2": omega_b, "n_t": tilt}
    for omega_b in (1.0e5, 5.0e5, 1.0e6)
    for tilt in (-0.5, -0.25, 0.0)
]
screen.HOLDOUTS = [
    {"Omega_bh2": 3.0e5, "n_t": -0.375},
    {"Omega_bh2": 7.5e5, "n_t": -0.125},
]
screen.OUTPUT = HERE / "screen_cr0_density_negative_tilt_exploratory.json"
screen.BASELINE_PARTIAL = HERE / "screen_cr0_density_negative_tilt_baseline_partial.json"


if __name__ == "__main__":
    screen.main()
