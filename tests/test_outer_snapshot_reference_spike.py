import numpy as np

from stiffgwpy_fast import fast_sgwb as FS


def test_outer_snapshot_keeps_generation_owned_arrays_by_reference():
    nv = np.arange(4.0)
    sigma = np.linspace(1.0, 2.0, 4)
    f_hor = np.linspace(3.0, 4.0, 4)

    snapshot = FS._retain_outer_background_snapshot(nv, sigma, f_hor)

    assert snapshot == (nv, sigma, f_hor)
    for before, after in zip(snapshot, (nv, sigma, f_hor)):
        assert before is after
