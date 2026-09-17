"""Contract test for the standalone FD lookup cache spike."""

from scripts.benchmark_fd_lookup_cache import cached_fd_lookup
from stiffgwpy_fast import exact_background as EB


def test_cached_fd_lookup_returns_the_existing_interpolators():
    expected = EB._fd_from_ref()
    actual = cached_fd_lookup()
    assert actual[0] is expected[0]
    assert actual[1] is expected[1]
