"""Small ``importlib.resources`` helpers for packaged scientific data."""

from contextlib import contextmanager
from importlib import resources
from pathlib import Path
from typing import Iterator


@contextmanager
def package_path(package: str, name: str) -> Iterator[Path]:
    """Yield a filesystem path for a package resource, including zipped wheels."""
    resource = resources.files(package).joinpath(name)
    with resources.as_file(resource) as path:
        yield path
