"""Test the root-accuracy boundary at a middle CR0 positive tilt."""

from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "validate_cr0_blue_tilt_amplitude_z8.py"
spec = importlib.util.spec_from_file_location("validate_cr0_blue_tilt_nt20", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the direct Reference runner")
validation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validation)
validation.CASES = (("middle_positive_tilt", 1.0e-10, 0.2),)
validation.OUTPUT = HERE / "validate_cr0_blue_tilt_amplitude_nt20_z8.json"


if __name__ == "__main__":
    validation.main()
