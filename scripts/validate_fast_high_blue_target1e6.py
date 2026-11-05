"""Run the high-blue Fast validation at a tighter adaptive-grid target."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs" / "fast_high_blue_target1e6_prufer_20261105.json"
sys.path.insert(0, str(ROOT))

from scripts import validate_fast_high_blue_fullgrid_prufer as validation  # noqa: E402


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    validation.TARGET = 1e-6
    validation.OUTPUT = OUTPUT
    validation.main()

    record = json.loads(OUTPUT.read_text(encoding="utf-8"))
    freqs = np.asarray(record["fast"]["frequencies_log10_hz"], dtype=float)
    spectrum = np.asarray(record["fast"]["omega_gw"], dtype=float)
    order = np.argsort(freqs)
    x = freqs[order]
    y = np.maximum(spectrum[order], 1e-300)
    curvature = np.abs(np.asarray(PchipInterpolator(x, np.log10(y)).derivative(2)(x)))
    widths = np.diff(x)
    local_error = curvature[:-1] * widths * widths / 8.0
    weights = 0.5 * (y[:-1] + y[1:]) * widths
    fractions = local_error * weights / np.sum(weights)
    record["fast"]["frequency_grid_error_sum_estimate"] = float(np.sum(fractions))
    record["fast"]["frequency_grid_error_max_panel"] = float(np.max(fractions)) if fractions.size else 0.0
    record["fast"]["frequency_grid_error_sum_note"] = (
        "Diagnostic sum of the existing curvature-weighted panel indicator; it is not an independent error bound."
    )
    record["hashes"]["runner_script_py"] = _sha(Path(__file__))
    record["hashes"]["runner_script_name"] = Path(__file__).name
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": record["status"],
        "grid_target": record["fast"]["target_dex"],
        "frequency_count": record["fast"]["frequency_count"],
        "grid_error_max_panel": record["fast"]["frequency_grid_error_max_panel"],
        "grid_error_sum_estimate": record["fast"]["frequency_grid_error_sum_estimate"],
        "output": str(OUTPUT),
    }), flush=True)


if __name__ == "__main__":
    main()
