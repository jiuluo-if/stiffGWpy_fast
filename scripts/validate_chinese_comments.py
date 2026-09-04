"""检查提交新增的代码注释是否包含中文。"""

from __future__ import annotations

import argparse
import io
import re
import subprocess
import sys
import tokenize
from pathlib import Path

CODE_SUFFIXES = {".js", ".py", ".sh", ".ts", ".yaml", ".yml"}
MACHINE_DIRECTIVES = (
    "# noqa",
    "# pragma:",
    "# pylint:",
    "# pyright:",
    "# type:",
)
CHINESE = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]")


def _python_comment(line: str) -> str | None:
    """提取 Python 行中的注释，避免把字符串中的井号误判为注释。"""
    try:
        tokens = tokenize.generate_tokens(io.StringIO(line).readline)
        for token in tokens:
            if token.type == tokenize.COMMENT:
                return token.string
    except (IndentationError, tokenize.TokenError):
        return None
    return None


def _comment_text(line: str, suffix: str) -> str | None:
    if suffix == ".py":
        return _python_comment(line)
    marker = line.find("#")
    return line[marker:] if marker >= 0 else None


def _base_revision(requested: str | None) -> str:
    candidate = requested or ""
    if not candidate or set(candidate) == {"0"}:
        candidate = "HEAD^"
    result = subprocess.run(
        ["git", "rev-parse", "--verify", f"{candidate}^{{commit}}"],
        check=False,
        capture_output=True,
        text=True,
    )
    return candidate if result.returncode == 0 else "HEAD^"


def _changed_comments(base: str, cached: bool) -> list[tuple[str, int, str]]:
    diff_command = ["git", "diff", "--unified=0"]
    if cached:
        diff_command.append("--cached")
    else:
        diff_command.extend([base, "HEAD"])
    result = subprocess.run(
        [*diff_command, "--"],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    path = ""
    new_line = 0
    findings: list[tuple[str, int, str]] = []
    for line in result.stdout.splitlines():
        if line.startswith("+++ b/"):
            path = line[6:]
            continue
        if line.startswith("@@"):
            match = re.search(r"\+(\d+)(?:,(\d+))?", line)
            if match:
                new_line = int(match.group(1))
            continue
        if not line.startswith("+") or line.startswith("+++"):
            if line.startswith(" "):
                new_line += 1
            continue
        suffix = Path(path).suffix.lower()
        comment = _comment_text(line[1:], suffix) if suffix in CODE_SUFFIXES else None
        if comment and not any(comment.lower().startswith(item) for item in MACHINE_DIRECTIVES):
            if not CHINESE.search(comment):
                findings.append((path, new_line, comment.strip()))
        new_line += 1
    return findings


def main() -> int:
    parser = argparse.ArgumentParser(description="校验提交新增代码注释的中文要求")
    parser.add_argument("--base", help="与当前提交比较的基线提交")
    parser.add_argument("--cached", action="store_true", help="检查暂存区差异")
    args = parser.parse_args()
    base = _base_revision(args.base)
    findings = _changed_comments(base, args.cached)
    if findings:
        print("发现不符合中文注释要求的新增代码注释：", file=sys.stderr)
        for path, line, comment in findings:
            print(f"{path}:{line}: {comment}", file=sys.stderr)
        return 1
    print(f"中文注释门禁通过（基线：{base}）。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
