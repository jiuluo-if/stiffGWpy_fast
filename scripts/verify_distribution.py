"""检查 sdist 和 wheel 是否只包含发布所需内容。"""

from __future__ import annotations

import argparse
import sys
import tarfile
import zipfile
from pathlib import Path

FORBIDDEN_PARTS = {".github", "build", "dist", "docs", "scripts", "tests"}
FORBIDDEN_NAMES = {"mcmc_compare.yaml"}
REQUIRED_WHEEL_PARTS = {"stiffgwpy/__init__.py", "stiffgwpy/fast_sgwb.py"}


def _is_forbidden(name: str) -> bool:
    parts = set(name.replace("\\", "/").split("/"))
    return bool(parts & FORBIDDEN_PARTS) or any(name.endswith(item) for item in FORBIDDEN_NAMES)


def _archive_names(path: Path) -> list[str]:
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            return archive.namelist()
    if path.name.endswith(".tar.gz"):
        with tarfile.open(path, "r:gz") as archive:
            return archive.getnames()
    raise ValueError(f"不支持的发布文件：{path}")


def verify(path: Path) -> list[str]:
    names = _archive_names(path)
    errors = [f"发布文件包含禁止内容：{name}" for name in names if _is_forbidden(name)]
    normalized = {name.replace("\\", "/") for name in names}
    if path.suffix == ".whl":
        errors.extend(
            f"wheel 不应包含构建元数据：{name}"
            for name in names
            if "stiffgwpy.egg-info" in name.replace("\\", "/").split("/")
        )
        missing = sorted(REQUIRED_WHEEL_PARTS - normalized)
        errors.extend(f"wheel 缺少必要文件：{name}" for name in missing)
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description="检查 sdist 和 wheel 的发布边界")
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    archives = sorted(args.directory.glob("*.whl")) + sorted(args.directory.glob("*.tar.gz"))
    if len(archives) != 2:
        print("发布目录必须恰好包含一个 wheel 和一个 sdist。", file=sys.stderr)
        return 1
    errors = [error for archive in archives for error in verify(archive)]
    if errors:
        print("发布包校验失败：", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("发布包校验通过：docs、tests、scripts 和历史研究配置均未进入归档。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
