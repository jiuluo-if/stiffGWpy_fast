import json
import runpy
from pathlib import Path

import pytest


def test_committed_validation_manifest_matches_schema():
    root = Path(__file__).resolve().parents[1]
    validator = runpy.run_path(str(root / "scripts" / "validate_manifest.py"))
    with (root / "docs" / "validation" / "validation_manifest.json").open(
        encoding="utf-8"
    ) as handle:
        manifest = json.load(handle)
    assert validator["validate_manifest"](manifest)


def test_optional_json_returns_none_only_when_file_is_missing(tmp_path):
    builder = runpy.run_path(
        str(Path(__file__).resolve().parents[1] / "scripts" / "build_two_mode_manifest.py")
    )

    assert builder["_load_optional_json"](tmp_path / "missing.json") is None

    malformed = tmp_path / "malformed.json"
    malformed.write_text("{", encoding="utf-8")
    with pytest.raises(json.JSONDecodeError):
        builder["_load_optional_json"](malformed)
