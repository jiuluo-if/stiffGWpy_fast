"""Independent Sobol Fast checks for the 13.9921875M CR0 box."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "screen_h2_identity_region.py"
SEED = 20261153
spec = importlib.util.spec_from_file_location("screen_cr0_density_to139921875_nt0302734375_sobol", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the shared paired Fast screen")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)
screen.BASELINE_TIMEOUT_SECONDS = 90.0
screen.BASE_PARAMS.update({
    "Omega_bh2": 13.875e6,
    "Omega_ch2": 0.1201075,
    "cr": 0.0,
    "r": 0.01,
    "T_re": 2000.0,
    "kappa10": 0.01,
})
screen.RECTANGLE = []
points = qmc.Sobol(d=2, scramble=True, seed=SEED).random_base2(m=3)
screen.HOLDOUTS = [
    {"Omega_bh2": 13.875e6 + float(p[0]) * 0.1171875e6,
     "n_t": -0.0302734375 + float(p[1]) * 0.0302734375}
    for p in points
]
screen.OUTPUT = HERE / "screen_cr0_density_13p875m_to13p9921875m_nt0302734375_slice_sobol.json"
screen.BASELINE_PARTIAL = HERE / "screen_cr0_density_13p875m_to13p9921875m_nt0302734375_slice_sobol_baseline_partial.json"


if __name__ == "__main__":
    screen.main()
    payload = json.loads(screen.OUTPUT.read_text(encoding="utf-8"))
    payload["sampling_design"] = {
        "method": "scipy.stats.qmc.Sobol", "dimension": 2, "scramble": True,
        "seed": SEED, "n": 8,
        "parameter_box": {"Omega_bh2": [13.875e6, 13.9921875e6], "n_t": [-0.0302734375, 0.0]},
    }
    helper = HERE.parents[2] / "stiffgwpy_fast" / "_wkb_transfer.py"
    payload["candidate_source_sha256"]["_wkb_transfer.py"] = hashlib.sha256(helper.read_bytes()).hexdigest()
    screen.OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
