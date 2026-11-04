from pathlib import Path

import pytest


@pytest.mark.parametrize("workflow", ["ci.yml", "slow.yml"])
def test_github_actions_workflow_uses_explicit_ubuntu_images(workflow):
    path = Path(__file__).resolve().parents[1] / ".github" / "workflows" / workflow
    assert "ubuntu-latest" not in path.read_text(encoding="utf-8")
