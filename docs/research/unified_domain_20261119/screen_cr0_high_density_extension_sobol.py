"""Independent scrambled Sobol holdouts for the CR0 high-density extension."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "screen_h2_identity_region.py"
SEED = 20261123
spec = importlib.util.spec_from_file_location("screen_cr0_high_density_sobol", RUNNER)
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
screen.RECTANGLE = []
samples = qmc.Sobol(d=2, scramble=True, seed=SEED).random_base2(m=3)
screen.HOLDOUTS = [
    {
        "Omega_bh2": 1.0e6 + float(point[0]) * 1.0e6,
        "n_t": -0.5 + float(point[1]) * 0.5,
    }
    for point in samples
]
screen.OUTPUT = HERE / "screen_cr0_high_density_extension_sobol.json"
screen.BASELINE_PARTIAL = HERE / "screen_cr0_high_density_extension_sobol_baseline_partial.json"


if __name__ == "__main__":
    screen.main()
    payload = json.loads(screen.OUTPUT.read_text(encoding="utf-8"))
    payload["sampling_design"] = {
        "method": "scipy.stats.qmc.Sobol",
        "dimension": 2,
        "scramble": True,
        "seed": SEED,
        "n": 8,
        "parameter_box": {"Omega_bh2": [1.0e6, 2.0e6], "n_t": [-0.5, 0.0]},
    }
    screen.OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
