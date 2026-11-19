"""Seeded scrambled log-r Sobol screen inside the CR0 positive-tilt strip."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "screen_h2_identity_region.py"
SEED = 20261122
spec = importlib.util.spec_from_file_location("screen_h2_region_cr0_blue_sobol", RUNNER)
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
samples = qmc.Sobol(d=2, scramble=True, seed=SEED).random_base2(m=3)
screen.RECTANGLE = []
screen.HOLDOUTS = [
    {"r": 10.0 ** (-12.0 + float(row[0]) * 4.0),
     "n_t": 0.05 + float(row[1]) * 0.05}
    for row in samples
]
screen.OUTPUT = HERE / "screen_cr0_blue_tilt_low_strip_sobol.json"
screen.BASELINE_PARTIAL = HERE / "screen_cr0_blue_tilt_low_strip_sobol_baseline_partial.json"


if __name__ == "__main__":
    screen.main()
    payload = json.loads(screen.OUTPUT.read_text(encoding="utf-8"))
    payload["sampling_design"] = {
        "method": "scipy.stats.qmc.Sobol",
        "dimension": 2,
        "scramble": True,
        "seed": SEED,
        "n": 8,
        "parameter_box": {"log10_r": [-12.0, -8.0], "n_t": [0.05, 0.1]},
    }
    screen.OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
