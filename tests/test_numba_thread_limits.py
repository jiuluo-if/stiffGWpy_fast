import os
import subprocess
import sys
from pathlib import Path

import pytest

from stiffgwpy_fast import fast_sgwb


@pytest.mark.parametrize(
    ("setter", "value"),
    [
        (fast_sgwb.set_threads, 1.5),
        (fast_sgwb.set_threads, True),
        (fast_sgwb.set_col_step, 2.5),
        (fast_sgwb.set_col_step, True),
    ],
)
def test_legacy_integer_setters_reject_fractional_values(setter, value):
    with pytest.raises(ValueError, match="integer"):
        setter(value)


def test_max_threads_uses_numba_config_after_mask_is_reduced():
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["NUMBA_NUM_THREADS"] = "4"
    env.pop("FAST_THREADS", None)
    code = (
        "import numba; "
        "numba.set_num_threads(1); "
        "from stiffgwpy_fast import fast_sgwb; "
        "assert fast_sgwb.max_threads() == numba.config.NUMBA_NUM_THREADS == 4; "
        "assert fast_sgwb.get_settings()['threads'] == 1"
    )

    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=root,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
