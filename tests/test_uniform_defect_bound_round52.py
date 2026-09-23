"""TDD contract for the Round 52 closed-form defect bound."""


def test_endpoint_defect_bound_is_zero_for_constant_coefficient():
    from scripts.benchmark_uniform_defect_bound_round52 import _defect_bound

    assert _defect_bound(0.6, 0.6, 0.005) == 0.0


def test_endpoint_defect_bound_is_symmetric_and_nonnegative():
    from scripts.benchmark_uniform_defect_bound_round52 import _defect_bound

    left = _defect_bound(-0.2, 0.7, 0.01)
    right = _defect_bound(0.7, -0.2, 0.01)
    assert left >= 0.0
    assert left == right


def test_endpoint_bound_dominates_embedded_matrix_defect_on_probes():
    from scripts.benchmark_uniform_asymptotic_defect_round51 import _uniform_block_defect
    from scripts.benchmark_uniform_defect_bound_round52 import _defect_bound

    for z0, z1, h in ((0.0, 0.01, 0.005), (0.5, 0.52, 0.005),
                      (2.0, 2.04, 0.005)):
        assert _defect_bound(z0, z1, h) >= _uniform_block_defect(z0, z1, h)
