from scripts.benchmark_res4_scalar_step import batch_case_rate


def test_batch_throughput_counts_all_variants_in_every_pair():
    assert batch_case_rate(pair_count=15, variant_count=3, elapsed_seconds=5.0) == 9.0
