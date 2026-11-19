"""Paired baseline/current Fast screen on the CR0 positive-tilt strip."""

from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "screen_h2_identity_region.py"
spec = importlib.util.spec_from_file_location("screen_h2_region_cr0_blue_strip", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the shared Fast comparison runner")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)

screen.BASE_PARAMS.update({
    "Omega_bh2": 1.0e6,
    "Omega_ch2": 0.1201075,
    "cr": 0.0,
    "T_re": 2000.0,
    "kappa10": 0.01,
})
screen.RECTANGLE = [
    {"r": tensor_r, "n_t": tilt}
    for tensor_r in (1.0e-12, 1.0e-10, 1.0e-8)
    for tilt in (0.05, 0.075, 0.1)
]
screen.HOLDOUTS = [
    {"r": 5.0e-11, "n_t": 0.0625},
    {"r": 3.0e-9, "n_t": 0.0875},
]
screen.OUTPUT = HERE / "screen_cr0_blue_tilt_low_strip.json"
screen.BASELINE_PARTIAL = HERE / "screen_cr0_blue_tilt_low_strip_baseline_partial.json"


if __name__ == "__main__":
    screen.main()
