import json
import runpy
from pathlib import Path


def test_committed_validation_manifest_matches_schema():
    root = Path(__file__).resolve().parents[1]
    validator = runpy.run_path(str(root / "scripts" / "validate_manifest.py"))
    with (root / "docs" / "validation" / "validation_manifest.json").open(
        encoding="utf-8"
    ) as handle:
        manifest = json.load(handle)
    assert validator["validate_manifest"](manifest)
