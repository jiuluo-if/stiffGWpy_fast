"""Direct z8 Reference closure at the interior of the expanded CR1 patch."""

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
OUTPUT = HERE / "validate_exact_quarter_transfer_1350_candidate.json"
PARAMS = {
    "Omega_bh2": 6.75e5, "Omega_ch2": 6.75e5, "H0": 67.32117,
    "DN_eff": 0.0, "A_s": 2.100549e-9, "r": 0.01, "n_t": 0.0,
    "cr": 1.0, "T_re": 2000.0, "DN_re": 10.0, "kappa10": 0.01,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    model = LCDM_SG(**PARAMS)
    if SGWB_iter_fast(
        model, config=resolve_config("fast", threads=2), tol=1e-6,
        physical_guard=True,
    ) is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    fast_root = float(model.DN_gw[-1])
    frequencies = np.asarray(model.f, dtype=np.float64)
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)

    start = time.perf_counter()
    first = run_reference(
        model, dn_eff=fast_root, freq_subset=frequencies, z_tail=8.0,
        rtol=1e-9, self_consistent=False, workers=2,
    )
    first_seconds = time.perf_counter() - start
    ref_freqs = np.asarray(first["freqs"], dtype=np.float64)
    if not np.array_equal(frequencies, ref_freqs):
        raise RuntimeError("direct Reference grid differs from candidate Fast grid")
    first_map = float(first["DN_gw"])
    ref_log = np.asarray(first["log10OmegaGW"], dtype=np.float64)
    ref_signal = np.asarray(first["Ogw"]) - np.asarray(first["Oj"])
    valid = (ref_signal > 0.0) & np.isfinite(fast_log) & np.isfinite(ref_log)
    spectrum_error = np.abs(fast_log[valid] - ref_log[valid])

    start = time.perf_counter()
    second = run_reference(
        model, dn_eff=first_map, freq_subset=frequencies, z_tail=8.0,
        rtol=1e-9, self_consistent=False, workers=2,
    )
    second_seconds = time.perf_counter() - start
    second_map = float(second["DN_gw"])
    hashes = {
        name: _sha256(ROOT / "stiffgwpy_fast" / name)
        for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
    }
    payload = {
        "params": PARAMS, "reference_code_modified": False,
        "fast_root": fast_root, "reference_map_at_fast_root": first_map,
        "reference_map_at_first_map": second_map,
        "reference_closure_residual_relative": abs(second_map-first_map)/max(abs(second_map),1e-300),
        "fast_root_relative_to_closed_reference": abs(fast_root-second_map)/max(abs(second_map),1e-300),
        "spectrum_max_abs_dex_at_fast_root": float(np.max(spectrum_error)),
        "spectrum_p95_abs_dex_at_fast_root": float(np.quantile(spectrum_error,0.95)),
        "valid_spectrum_nodes": int(valid.sum()), "frequency_grid_bitwise_equal": True,
        "reference_tail_z": 8.0, "first_reference_seconds": first_seconds,
        "closure_reference_seconds": second_seconds, "source_sha256": hashes,
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2))


if __name__ == "__main__":
    main()
