"""Contract test for the standalone pairwise Magnus map."""

import os
import subprocess
import sys


def test_pair_magnus_reports_finite_candidate_in_isolated_process():
    code = """
from scripts.benchmark_pair_magnus_round63 import _probe
row = _probe('default')
assert row['finite']
assert row['native_intervals'] > 100
assert row['native_intervals'] % 2 == 0
assert row['candidate_pair_maps'] == row['native_intervals'] // 2
assert row['candidate_power_relative_error'] > 1e-7
"""
    env = os.environ.copy()
    env.update({"NUMBA_NUM_THREADS": "2", "NUMBA_THREADING_LAYER": "workqueue"})
    subprocess.run([sys.executable, "-c", code], check=True, env=env)
