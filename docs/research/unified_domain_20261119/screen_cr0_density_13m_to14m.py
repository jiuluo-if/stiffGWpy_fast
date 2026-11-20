"""Paired Fast screen on the CR0 density continuation from 13M to 14M."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "screen_h2_identity_region.py"
spec = importlib.util.spec_from_file_location("screen_cr0_density_13m_to14m", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the shared paired Fast screen")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)

screen.BASELINE_TIMEOUT_SECONDS = 90.0
screen.BASE_PARAMS.update({
    "Omega_bh2": 13.0e6,
    "Omega_ch2": 0.1201075,
    "cr": 0.0,
    "r": 0.01,
    "T_re": 2000.0,
    "kappa10": 0.01,
})
screen.RECTANGLE = [
    {"Omega_bh2": density, "n_t": tilt}
    for density in (13.0e6, 13.5e6, 14.0e6)
    for tilt in (-0.5, -0.25, 0.0)
]
screen.HOLDOUTS = [
    {"Omega_bh2": 13.25e6, "n_t": -0.375},
    {"Omega_bh2": 13.75e6, "n_t": -0.125},
]
screen.OUTPUT = HERE / "screen_cr0_density_13m_to14m.json"
screen.BASELINE_PARTIAL = HERE / "screen_cr0_density_13m_to14m_baseline_partial.json"


if __name__ == "__main__":
    screen.main()
    payload = json.loads(screen.OUTPUT.read_text(encoding="utf-8"))
    helper = HERE.parents[2] / "stiffgwpy_fast" / "_wkb_transfer.py"
    payload["candidate_source_sha256"]["_wkb_transfer.py"] = hashlib.sha256(helper.read_bytes()).hexdigest()
    screen.OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
