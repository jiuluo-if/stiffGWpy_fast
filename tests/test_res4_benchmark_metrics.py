from scripts.benchmark_res4_scalar_step import batch_case_rate
from scripts.benchmark_sample11_reference_map import align_spectrum


def test_batch_throughput_counts_all_variants_in_every_pair():
    assert batch_case_rate(pair_count=15, variant_count=3, elapsed_seconds=5.0) == 9.0


def test_reference_spectrum_alignment_handles_descending_fast_frequencies():
    fast_log = align_spectrum(
        fast_freqs=[3.0, 2.0, 1.0], fast_log=[30.0, 20.0, 10.0],
        reference_freqs=[1.0, 3.0],
    )
    assert fast_log.tolist() == [10.0, 30.0]
