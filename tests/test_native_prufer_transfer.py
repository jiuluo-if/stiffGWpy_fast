import numpy as np

from scripts.benchmark_native_prufer_transfer import (
    _cartesian_chain,
    _prufer_chain,
)


def test_native_prufer_transfer_is_finite_and_power_close():
    cart = _cartesian_chain(0.2, 3.0, 0.01, 256)
    prufer = _prufer_chain(0.2, 3.0, 0.01, 256)
    assert np.isfinite(np.asarray(prufer)).all()
    assert abs(prufer[0] * prufer[0] -
               (cart[0] * cart[0] + cart[1] * cart[1])) < 1e-8
