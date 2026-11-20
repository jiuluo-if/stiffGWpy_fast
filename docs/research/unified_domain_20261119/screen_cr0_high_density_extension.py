"""CR0 high-density continuation screen, paired against pre-stable baseline."""

from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "screen_h2_identity_region.py"
spec = importlib.util.spec_from_file_location("screen_cr0_high_density", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the shared paired Fast screen")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)

screen.BASE_PARAMS.update(
    {
        "Omega_bh2": 1.0e6,
        "Omega_ch2": 0.1201075,
        "cr": 0.0,
        "r": 0.01,
        "T_re": 2000.0,
        "kappa10": 0.01,
    }
)
screen.RECTANGLE = [
    {"Omega_bh2": density, "n_t": tilt}
    for density in (1.0e6, 1.5e6, 2.0e6)
    for tilt in (-0.5, -0.25, 0.0)
]
screen.HOLDOUTS = [
    {"Omega_bh2": 1.25e6, "n_t": -0.375},
    {"Omega_bh2": 1.75e6, "n_t": -0.125},
]
screen.OUTPUT = HERE / "screen_cr0_high_density_extension.json"
screen.BASELINE_PARTIAL = HERE / "screen_cr0_high_density_extension_baseline_partial.json"


if __name__ == "__main__":
    screen.main()
