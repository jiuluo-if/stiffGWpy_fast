"""Contract test for the event-aware tangent endpoint prototype."""

import os
import subprocess
import sys


def test_tangent_endpoint_is_finite_and_small_error_on_default_mode():
    code = """
from scripts.benchmark_outer_sensitivity_round66 import _probe
row = _probe('default', mode=20, repeats=3)
assert row['finite'], row
assert row['eligible'], row
assert row['state_relative_error'] < 1e-3, row
"""
    env = os.environ.copy()
    env.update({"NUMBA_NUM_THREADS": "2", "NUMBA_THREADING_LAYER": "workqueue"})
    subprocess.run([sys.executable, "-c", code], check=True, env=env)
