"""Close one Reference z_tail=8 map step for the archived feature-band grid."""

# The script is run by path and adds the repository root before importing the
# local checkout rather than whichever package version is installed.
# ruff: noqa: E402, I001

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.reference import run_reference
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


RESEARCH = Path(__file__).resolve().parent
INPUT = RESEARCH / "exact_matter_transfer_feature_probe.json"
PARAMS = {
    "cr": 0.0,
    "T_re": 10.0,
    "kappa10": 1e-6,
    "n_t": 1.0,
    "r": 0.01,
    "A_s": 1.9680104830424577e-7,
    "DN_re": 740.0,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    previous = json.loads(INPUT.read_text(encoding="utf-8"))
    dn = float(previous["reference_map_at_fast_root"])
    frequencies = np.asarray(previous["frequency"], dtype=float)
    start = time.perf_counter()
    mapped = run_reference(
        LCDM_SG(**PARAMS),
        dn_eff=dn,
        freq_subset=frequencies,
        z_tail=8.0,
        workers=1,
        rtol=1e-11,
    )
    elapsed = time.perf_counter() - start
    mapped_root = float(mapped["DN_gw"])
    logspec = np.log10(np.asarray(mapped["Ogw"] - mapped["Oj"], dtype=float))
    record = {
        "parameters": PARAMS,
        "z_tail": 8.0,
        "frequency": frequencies.tolist(),
        "dn_input": dn,
        "reference_map": mapped_root,
        "fixed_point_residual_relative": abs(mapped_root - dn) / max(abs(dn), 1e-300),
        "spectrum_log10omega": logspec.tolist(),
        "reference_seconds": elapsed,
        "previous_artifact_sha256": sha256(INPUT),
        "source_sha256": sha256(ROOT / "stiffgwpy_fast" / "reference.py"),
        "source_head": __import__("subprocess").check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
    }
    (RESEARCH / "feature_reference_root_z8_next_map.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in record.items() if k not in {
        "frequency", "spectrum_log10omega",
    }}, indent=2))


if __name__ == "__main__":
    main()
