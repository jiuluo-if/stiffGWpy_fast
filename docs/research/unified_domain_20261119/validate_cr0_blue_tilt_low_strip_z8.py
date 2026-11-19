"""Direct z8 checks for the lower edge and Sobol holdout of the CR0 strip."""

from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "validate_cr0_blue_tilt_amplitude_z8.py"
spec = importlib.util.spec_from_file_location("validate_cr0_blue_tilt_low_strip", RUNNER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the direct Reference runner")
validation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validation)
validation.CASES = (
    ("low_amplitude_lower_tilt_edge", 1.0e-12, 0.05),
    ("high_amplitude_lower_tilt_edge", 1.0e-8, 0.05),
    ("independent_sobol_holdout", 3.96526332334898e-9, 0.07029469762928785),
)
validation.OUTPUT = HERE / "validate_cr0_blue_tilt_low_strip_z8.json"


if __name__ == "__main__":
    validation.main()
