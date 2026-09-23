"""Contract tests for the standalone coefficient-exp recurrence."""

import os
import subprocess
import sys


def _run_isolated_check():
    code = """
import numpy as np
from scripts.benchmark_exp_recurrence_round61 import (
    _path_baseline,
    _path_candidate,
)
starts = np.array([0.2, 0.2025, 0.205], dtype=np.float64)
ends = starts + 0.0025
steps = np.full(3, 0.005, dtype=np.float64)
baseline = _path_baseline(starts, ends, steps, 0.25)
candidate = _path_candidate(starts, ends, steps, 0.25)
assert np.isfinite(candidate).all()
assert np.allclose(candidate, baseline, rtol=2e-12, atol=2e-13)
"""
    env = os.environ.copy()
    env.update({"NUMBA_NUM_THREADS": "2", "NUMBA_THREADING_LAYER": "workqueue"})
    subprocess.run([sys.executable, "-c", code], check=True, env=env)


def test_exp_recurrence_is_finite_and_close_on_a_subdivided_path():
    _run_isolated_check()
