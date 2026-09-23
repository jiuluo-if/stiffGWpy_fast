"""TDD contract for the Round 54 exact-composition Amdahl screen."""


def test_exact_composition_keeps_transcendental_count_and_adds_algebra():
    from scripts.benchmark_exact_block_schedule_round54 import _operation_model

    model = _operation_model(native_intervals=8, substeps=3, block_size=2)
    assert model["baseline_transcendentals"] == 24
    assert model["candidate_transcendentals"] == 24
    assert model["candidate_extra_matrix_products"] > 0


def test_singleton_schedule_has_no_composition_overhead():
    from scripts.benchmark_exact_block_schedule_round54 import _operation_model

    model = _operation_model(native_intervals=4, substeps=1, block_size=1)
    assert model["candidate_extra_matrix_products"] == 0
