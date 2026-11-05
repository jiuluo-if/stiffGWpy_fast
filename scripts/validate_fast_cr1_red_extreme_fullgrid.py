"""Full-grid Prüfer validation of the high-r strongly red positive-cr point."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import validate_fast_high_blue_fullgrid_prufer as validation  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

OUTPUT = ROOT / "docs" / "fast_cr1_red_extreme_fullgrid_prufer_20261105.json"
PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 10.0,
    "n_t": 2.0,
    "cr": 1.0,
    "T_re": 1e-3,
    "DN_re": 200.0,
    "kappa10": 100.0,
}


def main() -> None:
    validation.PARAMS = PARAMS
    validation.TARGET = 1e-7
    validation.OUTPUT = OUTPUT
    validation.main()

    record = json.loads(OUTPUT.read_text(encoding="utf-8"))
    model = LCDM_SG(**PARAMS)
    derived = model.derived_param
    record["derived"]["N_re"] = float(derived["N_re"])
    record["derived"]["DN_re_effective"] = float(model.cosmo_param["DN_re"])
    record["hashes"]["runner_script_py"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    record["hashes"]["runner_script_name"] = Path(__file__).name
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
