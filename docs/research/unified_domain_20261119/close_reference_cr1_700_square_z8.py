"""Close one additional z8 Reference map at the high-density boundary probe."""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HERE = Path(__file__).resolve().parent
SOURCE_RESULT = HERE / "validate_h2_identity_cr1_700_square_z8.json"
OUTPUT = HERE / "close_reference_cr1_700_square_z8.json"
PARAMS = {
    "Omega_bh2": 7.0e5, "Omega_ch2": 7.0e5, "H0": 67.32117,
    "DN_eff": 0.0, "A_s": 2.100549e-9, "r": 0.01, "n_t": 0.0,
    "cr": 1.0, "T_re": 2000.0, "DN_re": 10.0, "kappa10": 0.01,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    stored = json.loads(SOURCE_RESULT.read_text(encoding="utf-8"))
    row = stored["cases"][0]
    model = LCDM_SG(**PARAMS)
    if SGWB_iter_fast(
        model, config=resolve_config("fast", threads=2), tol=1e-6,
        physical_guard=True,
    ) is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    fast_root = float(model.DN_gw[-1])
    frequencies = np.asarray(model.f, dtype=np.float64)
    if not np.array_equal(frequencies, np.asarray(row["raw"]["freqs_log10_hz"])):
        raise RuntimeError("Fast grid differs from the saved first Reference map")

    first_map = float(row["reference_map_at_fast_root"])
    start = time.perf_counter()
    closed = run_reference(
        model, dn_eff=first_map, freq_subset=frequencies, z_tail=8.0,
        rtol=1e-9, self_consistent=False, workers=2,
    )
    seconds = time.perf_counter() - start
    second_map = float(closed["DN_gw"])
    ref_signal = np.asarray(closed["Ogw"]) - np.asarray(closed["Oj"])
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    ref_log = np.asarray(closed["log10OmegaGW"], dtype=np.float64)
    valid = (ref_signal > 0.0) & np.isfinite(fast_log) & np.isfinite(ref_log)
    error = np.abs(fast_log[valid] - ref_log[valid])
    source_hashes = {
        name: _sha256(ROOT / "stiffgwpy_fast" / name)
        for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
    }
    if source_hashes["reference.py"] != stored["source_sha256"]["reference.py"]:
        raise RuntimeError("Reference source changed relative to the archived map")
    payload = {
        "params": PARAMS, "reference_code_modified": False,
        "z_tail": 8.0, "same_frequency_grid": True,
        "fast_root": fast_root, "first_reference_map_at_fast_root": first_map,
        "second_reference_map_at_first_map": second_map,
        "reference_fixed_point_residual_relative": abs(second_map-first_map)/max(abs(second_map),1e-300),
        "fast_root_relative_to_second_map": abs(fast_root-second_map)/max(abs(second_map),1e-300),
        "second_map_seconds": seconds,
        "second_map_spectrum_max_abs_dex_vs_fast": float(np.max(error)),
        "second_map_spectrum_p95_abs_dex_vs_fast": float(np.quantile(error,0.95)),
        "source_sha256": source_hashes,
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2))


if __name__ == "__main__":
    main()
