"""Small ``importlib.resources`` helpers for packaged scientific data.

中文：统一通过 `importlib.resources` 读取 wheel 内的数据，避免假设资源一定是源码目录中的普通文件。
"""

from contextlib import contextmanager
from importlib import resources
from pathlib import Path
from typing import Iterator


@contextmanager
def package_path(package: str, name: str) -> Iterator[Path]:
    """Yield a filesystem path for a package resource, including zipped wheels.

    中文：返回的路径只在 `with` 上下文中有效；`as_file` 可能需要临时释放压缩包资源。
    """
    resource = resources.files(package).joinpath(name)
    with resources.as_file(resource) as path:
        yield path
