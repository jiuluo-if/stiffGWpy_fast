"""TDD contract for the tangent propagation feasibility screen."""

import json
import os
import subprocess
import sys


def test_tangent_step_matches_finite_difference():
    code = (
        "import json; "
        "from scripts.benchmark_tangent_outer_spike import _local_jacobian_error; "
        "print(json.dumps(_local_jacobian_error()))"
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
    assert json.loads(result.stdout.splitlines()[-1]) < 1.0e-8
