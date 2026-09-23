"""Contract test for the standalone fourth-order CF Magnus spike."""

import json
import os
import subprocess
import sys


def test_cf4_magnus_zero_threshold_is_bitwise_fallback():
    code = (
        "import json; "
        "from scripts.benchmark_cf4_magnus_spike import _run_kernel; "
        "print(json.dumps(_run_kernel('default', threshold=0.0, repeats=1, threads=2)))"
    )
    clean_env = os.environ.copy()
    for name in (
        "NUMBA_NUM_THREADS", "NUMBA_THREADING_LAYER", "FAST_THREADS",
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        clean_env.pop(name, None)
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        capture_output=True,
        text=True,
        env=clean_env,
    )
    row = json.loads(result.stdout.splitlines()[-1])
    assert row["bitwise_equal"] is True
    assert row["full_output_bitwise_equal"] is True
    assert row["accepted_blocks"] == 0
