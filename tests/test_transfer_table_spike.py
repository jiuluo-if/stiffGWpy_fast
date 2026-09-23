"""TDD contract for the combined transfer-map table feasibility screen."""

import json
import os
import subprocess
import sys


def test_transfer_table_screen_detects_insufficient_local_error():
    code = (
        "import json; "
        "from scripts.benchmark_transfer_table_spike import _local_table_gate; "
        "print(json.dumps(_local_table_gate()))"
    )
    clean_env = os.environ.copy()
    for name in (
        "NUMBA_NUM_THREADS", "NUMBA_THREADING_LAYER", "FAST_THREADS",
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        clean_env.pop(name, None)
    result = subprocess.run(
        [sys.executable, "-c", code], check=True, capture_output=True,
        text=True, env=clean_env)
    metrics = json.loads(result.stdout.splitlines()[-1])
    assert metrics["max_transfer_absolute_error"] >= 1.0e-11
