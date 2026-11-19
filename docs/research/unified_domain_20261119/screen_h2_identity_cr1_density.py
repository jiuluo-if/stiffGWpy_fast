"""Paired Fast screen on a consistency-relation (cr>0) density rectangle."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUNNER = ROOT / "docs" / "research" / "unified_domain_20261119" / "screen_h2_identity_region.py"
spec = importlib.util.spec_from_file_location("h2_region_screen_cr1", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load shared region screen runner")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)

screen.BASE_PARAMS.update({
    "Omega_bh2": 1.0e5,
    "Omega_ch2": 1.0e5,
    "T_re": 2000.0,
    "cr": 1.0,
})
screen.RECTANGLE = [
    {"Omega_bh2": ob, "Omega_ch2": oc}
    for ob in (1.0e5, 3.0e5, 1.0e6)
    for oc in (1.0e5, 3.0e5, 1.0e6)
]
screen.HOLDOUTS = [
    {"Omega_bh2": 5.0e5, "Omega_ch2": 5.0e5},
    {"Omega_bh2": 7.5e5, "Omega_ch2": 2.5e5},
]
screen.OUTPUT = Path(__file__).with_name("screen_h2_identity_cr1_density.json")
screen.BASELINE_PARTIAL = Path(__file__).with_name(
    "screen_h2_identity_cr1_density_baseline_partial.json"
)


if __name__ == "__main__":
    screen.main()
