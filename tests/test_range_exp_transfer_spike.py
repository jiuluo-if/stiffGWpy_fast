"""TDD contract for the range-reduced exponential transfer spike."""

import json
import os
import subprocess
import sys


def test_range_reduced_exp_is_locally_accurate():
    code = (
        "import json; "
        "from scripts.benchmark_range_exp_transfer_spike import _local_errors; "
        "print(json.dumps(_local_errors()))"
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
    assert metrics["exp_relative_max"] < 5.0e-15
    assert metrics["transfer_absolute_max"] < 5.0e-14
