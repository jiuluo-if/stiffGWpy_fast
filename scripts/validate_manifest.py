"""Validate the committed validation manifest's stable structural contract."""

import argparse
import json
import re
from pathlib import Path

COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
PROFILE_KEYS = {
    "fast_plain_grid": "plain-grid",
    "fast_transition_refine": "transition-refine",
}
CONFIG_KEYS = {
    "h",
    "col_step",
    "z_tail",
    "freq_res",
    "transition_refine",
    "phase_max",
    "freq_grid",
    "outer_tol",
}
STATUSES = {"VERIFIED", "NOT VERIFIED", "PARTIALLY VERIFIED", "FAILED"}


def validate_manifest(manifest):
    """Raise ``ValueError`` when *manifest* violates the committed contract."""
    if not isinstance(manifest, dict):
        raise ValueError("manifest must be a JSON object")
    if manifest.get("schema_version") != 1:
        raise ValueError("unsupported schema_version; expected 1")
    commit = manifest.get("commit")
    if not isinstance(commit, str) or not COMMIT_RE.fullmatch(commit):
        raise ValueError("commit must be a 40-character lowercase git SHA")
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(manifest.get("date", ""))):
        raise ValueError("date must use YYYY-MM-DD")
    if not isinstance(manifest.get("generated_by"), str):
        raise ValueError("generated_by is required")
    semantics = manifest.get("oracle_semantics", "")
    if "continuous-sigma" not in semantics or "LSODA" not in semantics:
        raise ValueError("oracle_semantics must state the independent oracle and LSODA role")

    for key, expected_profile in PROFILE_KEYS.items():
        profile = manifest.get(key)
        if not isinstance(profile, dict):
            raise ValueError("missing profile: %s" % key)
        if profile.get("profile") != expected_profile:
            raise ValueError("%s has an unexpected profile name" % key)
        config = profile.get("config")
        if not isinstance(config, dict) or not CONFIG_KEYS.issubset(config):
            raise ValueError("%s is missing required config fields" % key)
        if config["freq_grid"] not in {"construct", "grid_independent", "adaptive"}:
            raise ValueError("%s has an invalid freq_grid" % key)
        if not isinstance(config["transition_refine"], bool):
            raise ValueError("%s transition_refine must be boolean" % key)
        counts = profile.get("status_counts")
        if not isinstance(counts, dict) or not isinstance(counts.get("total"), int):
            raise ValueError("%s status_counts.total must be an integer" % key)
        if profile.get("status") not in STATUSES:
            raise ValueError("%s has an invalid status" % key)
        if not isinstance(profile.get("accuracy"), dict):
            raise ValueError("%s accuracy section is required" % key)

    return True


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "manifest",
        nargs="?",
        type=Path,
        default=Path("docs/validation/validation_manifest.json"),
    )
    args = parser.parse_args(argv)
    with args.manifest.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    validate_manifest(manifest)
    print("validation manifest: OK", args.manifest)


if __name__ == "__main__":
    main()
