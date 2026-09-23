"""Contract test for the direct Ermakov--Kummer prototype."""

import os
import subprocess
import sys


def test_direct_phase_function_is_finite_and_stable_in_an_isolated_process():
    code = """
from scripts.benchmark_phase_function_direct_round62 import _probe
row = _probe('default')
assert row['finite']
assert row['rho_min'] > 0.0
assert row['candidate_vs_cartesian_relative_error'] > 1e-2
assert row['candidate_vs_cartesian_relative_error'] < 1e4
assert row['phase_state_stability_relative_error'] < 1e-8, row
"""
    env = os.environ.copy()
    env.update({"NUMBA_NUM_THREADS": "2", "NUMBA_THREADING_LAYER": "workqueue"})
    subprocess.run([sys.executable, "-c", code], check=True, env=env)
