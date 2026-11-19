"""Validate an asymmetric Fast holdout against same-sum raw Reference data."""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

REFERENCE = HERE / "validate_h2_identity_cr1_upper_boundary_z8.json"
OUTPUT = HERE / "validate_h2_identity_cr1_600_square_holdout_symmetry.json"
PARAMS = {
    "Omega_bh2": 5.5e5,
    "Omega_ch2": 5.5e5,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 1.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    model = LCDM_SG(**PARAMS)
    start = time.perf_counter()
    result = SGWB_iter_fast(
        model, config=resolve_config("fast", threads=2), tol=1e-6,
        physical_guard=True,
    )
    fast_seconds = time.perf_counter() - start
    if result is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")

    reference_artifact = json.loads(REFERENCE.read_text(encoding="utf-8"))
    reference_row = reference_artifact["cases"][0]
    reference_raw = reference_row["raw"]
    frequencies = np.asarray(model.f, dtype=np.float64)
    reference_frequencies = np.asarray(reference_raw["freqs_log10_hz"], dtype=np.float64)
    if not np.array_equal(frequencies, reference_frequencies):
        raise RuntimeError("same-sum density frequency grids are not bitwise equal")
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    reference_log = np.asarray(reference_raw["reference_log10_OmegaGW"], dtype=np.float64)
    if not np.isfinite(fast_log).all() or not np.isfinite(reference_log).all():
        raise RuntimeError("Fast or same-sum Reference spectrum is nonfinite")
    errors = np.abs(fast_log - reference_log)
    reference_source_sha = reference_artifact["source_sha256"]["reference.py"]
    source_hashes = {
        name: _sha256(HERE.parents[2] / "stiffgwpy_fast" / name)
        for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
    }
    if source_hashes["reference.py"] != reference_source_sha:
        raise RuntimeError("Reference source hash differs from the reused raw artifact")

    payload = {
        "case_id": "square_independent_holdout_1100k",
        "params": PARAMS,
        "same_total_density_reference_input": reference_row["params"],
        "density_sum": float(PARAMS["Omega_bh2"] + PARAMS["Omega_ch2"]),
        "reference_code_modified": False,
        "fast_root": float(model.DN_gw[-1]),
        "reference_map_at_same_total_density_and_root": float(
            reference_row["reference_map_at_fast_root"]
        ),
        "root_relative_error_proxy": float(
            abs(float(model.DN_gw[-1]) - reference_row["reference_map_at_fast_root"])
            / reference_row["reference_map_at_fast_root"]
        ),
        "spectrum_max_abs_dex": float(np.max(errors)),
        "spectrum_p95_abs_dex": float(np.quantile(errors, 0.95)),
        "fast_seconds": fast_seconds,
        "frequency_arrays_bitwise_equal": True,
        "reference_tail_z": 8.0,
        "fast_log10_OmegaGW": fast_log.tolist(),
        "reference_log10_OmegaGW": reference_log.tolist(),
        "freqs_log10_hz": frequencies.tolist(),
        "reference_source_sha256": reference_source_sha,
        "source_sha256": source_hashes,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        key: value for key, value in payload.items()
        if key not in ("fast_log10_OmegaGW", "reference_log10_OmegaGW", "freqs_log10_hz")
    }, indent=2))


if __name__ == "__main__":
    main()
