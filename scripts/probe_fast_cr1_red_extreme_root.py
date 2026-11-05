"""Two-grid secant screen for the high-r strongly red positive-cr point."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import probe_fast_high_blue_root_secant as screen  # noqa: E402

BASELINE = ROOT / "docs" / "fast_cr1_red_extreme_fullgrid_prufer_20261105.json"
OUTPUT = ROOT / "docs" / "fast_cr1_red_extreme_root_secant_20261105.json"


def main() -> None:
    screen.BASELINE = BASELINE
    screen.OUTPUT = OUTPUT
    screen.main()
    record = json.loads(OUTPUT.read_text(encoding="utf-8"))
    record["runner_script_py"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    record["runner_script_name"] = Path(__file__).name
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
