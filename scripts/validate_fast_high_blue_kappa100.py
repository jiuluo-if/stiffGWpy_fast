"""Validate the joint high-blue, high-stiff-density Fast corner."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import validate_fast_high_blue_fullgrid_prufer as validation  # noqa: E402

OUTPUT = ROOT / "docs" / "fast_high_blue_kappa100_prufer_20261105.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    validation.PARAMS = {**validation.PARAMS, "kappa10": 100.0}
    validation.TARGET = 1e-7
    validation.OUTPUT = OUTPUT
    validation.main()

    record = json.loads(OUTPUT.read_text(encoding="utf-8"))
    f = np.asarray(record["fast"]["frequencies_log10_hz"], dtype=float)
    y = np.asarray(record["fast"]["omega_gw"], dtype=float)
    order = np.argsort(f)
    x, y = f[order], np.maximum(y[order], 1e-300)
    c2 = np.abs(np.asarray(PchipInterpolator(x, np.log10(y)).derivative(2)(x)))
    h = np.diff(x)
    err = c2[:-1] * h * h / 8.0
    weights = 0.5 * (y[:-1] + y[1:]) * h
    frac = err * weights / np.sum(weights)
    record["fast"]["frequency_grid_error_sum_estimate"] = float(np.sum(frac))
    record["fast"]["frequency_grid_error_max_panel"] = float(np.max(frac)) if frac.size else 0.0
    record["fast"]["frequency_grid_error_sum_note"] = (
        "Diagnostic sum of the existing curvature-weighted panel indicator; not an independent error bound."
    )
    record["hashes"]["runner_script_py"] = _sha(Path(__file__))
    record["hashes"]["runner_script_name"] = Path(__file__).name
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": record["status"],
        "parameters": {k: record["parameters"][k] for k in ("r", "n_t", "kappa10", "cr")},
        "N_inf": record["derived"]["N_inf"],
        "frequency_count": record["fast"]["frequency_count"],
        "DN_gw": record["fast"]["DN_gw"],
        "acceptance": record["acceptance"],
        "output": str(OUTPUT),
    }), flush=True)


if __name__ == "__main__":
    main()
