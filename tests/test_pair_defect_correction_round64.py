"""Contract test for the standalone Duhamel defect-correction map."""

import os
import subprocess
import sys


def test_pair_defect_correction_is_finite_and_reduces_map_count():
    code = """
from scripts.benchmark_pair_defect_correction_round64 import _probe
row = _probe('default')
assert row['finite']
assert row['native_intervals'] % 2 == 0
assert row['candidate_pair_maps'] == row['native_intervals'] // 2
assert row['candidate_power_relative_error'] > 1e-8
"""
    env = os.environ.copy()
    env.update({"NUMBA_NUM_THREADS": "2", "NUMBA_THREADING_LAYER": "workqueue"})
    subprocess.run([sys.executable, "-c", code], check=True, env=env)
