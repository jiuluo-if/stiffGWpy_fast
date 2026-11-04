import pytest

from scripts.benchmark_res4_scalar_step import batch_case_rate, per_mode_error_record
from scripts.benchmark_sample11_reference_map import align_spectrum


def test_batch_throughput_counts_all_variants_in_every_pair():
    assert batch_case_rate(pair_count=15, variant_count=3, elapsed_seconds=5.0) == 9.0


def test_reference_spectrum_alignment_handles_descending_fast_frequencies():
    fast_log = align_spectrum(
        fast_freqs=[3.0, 2.0, 1.0], fast_log=[30.0, 20.0, 10.0],
        reference_freqs=[1.0, 3.0],
    )
    assert fast_log.tolist() == [10.0, 30.0]


def test_per_mode_error_record_preserves_reference_and_all_variant_deltas():
    record = per_mode_error_record(
        frequencies=[-3.0, -2.0],
        reference_log=[-10.0, -20.0],
        variants={"default": [-10.1, -19.8], "candidate": [-10.0, -20.0]},
    )
    assert record["log10_f_Hz"] == [-3.0, -2.0]
    assert record["reference_log10OmegaGW"] == [-10.0, -20.0]
    assert record["absolute_error_dex"]["default"] == pytest.approx([0.1, 0.2])
    assert record["absolute_error_dex"]["candidate"] == [0.0, 0.0]
