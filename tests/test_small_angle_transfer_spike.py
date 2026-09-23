"""TDD contract for the small-angle polynomial transfer spike."""

import json
import os
import subprocess
import sys


def test_small_angle_polynomial_is_locally_accurate():
    code = (
        "import json; "
        "from scripts.benchmark_small_angle_transfer_spike import _small_angle_error; "
        "print(json.dumps(_small_angle_error()))"
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
    assert json.loads(result.stdout.splitlines()[-1]) < 5.0e-15
