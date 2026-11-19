"""Separate Fast tensor-spectrum error from frequency-quadrature error."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as GP  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "validate_h2_identity_cr1_700_square_z8.json"
OUTPUT = HERE / "diagnose_cr1_1400k_spectrum_decomposition.json"
PARAMS = {
    "Omega_bh2": 7.0e5, "Omega_ch2": 7.0e5, "H0": 67.32117,
    "DN_eff": 0.0, "A_s": 2.100549e-9, "r": 0.01, "n_t": 0.0,
    "cr": 1.0, "T_re": 2000.0, "DN_re": 10.0, "kappa10": 0.01,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    model = LCDM_SG(**PARAMS)
    if FS.SGWB_iter_fast(
        model, config=FS.resolve_config("fast", threads=2), tol=1e-6,
        physical_guard=True,
    ) is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    row = reference["cases"][0]
    raw = row["raw"]
    freqs = np.asarray(model.f, dtype=np.float64)
    if not np.array_equal(freqs, np.asarray(raw["freqs_log10_hz"], dtype=np.float64)):
        raise RuntimeError("Fast grid differs from the saved direct Reference grid")
    fast_signal = np.asarray(model.Ogw_today) - np.asarray(model.Oj_today)
    ref_signal = np.asarray(raw["reference_Ogw"]) - np.asarray(raw["reference_Oj"])
    order = np.argsort(freqs)
    x = freqs[order]
    fast_intervals = FS._pchip_integrals_vectorized(x, fast_signal[order])
    ref_intervals = FS._pchip_integrals_vectorized(x, ref_signal[order])
    scale = GP.Neff0 / (GP.Omega_nh2 / model.derived_param["h"]**2) * GP.ln10
    fast_map = float(np.sum(fast_intervals) * scale)
    ref_pchip_map = float(np.sum(ref_intervals) * scale)
    delta = (fast_intervals - ref_intervals) * scale
    top = np.argsort(-np.abs(delta))[:10]
    records = [
        {
            "logf_lower": float(x[i]), "logf_upper": float(x[i+1]),
            "fast_minus_reference_map_contribution": float(delta[i]),
        }
        for i in top
    ]
    ref_map = float(row["reference_map_at_fast_root"])
    payload = {
        "params": PARAMS, "density_sum": 1.4e6,
        "frequency_grid_bitwise_equal": True, "frequency_nodes": int(x.size),
        "fast_root": float(model.DN_gw[-1]),
        "fast_map_rebuilt_from_fast_signal_pchip": fast_map,
        "reference_map_rebuilt_from_reference_signal_pchip": ref_pchip_map,
        "reference_map_from_direct_quadrature": ref_map,
        "fast_minus_reference_pchip_map": fast_map-ref_pchip_map,
        "reference_pchip_minus_reference_quadpack": ref_pchip_map-ref_map,
        "spectrum_max_abs_dex": float(np.max(np.abs(
            np.asarray(model.log10OmegaGW)-np.asarray(raw["reference_log10_OmegaGW"])
        ))),
        "largest_interval_contributions": records,
        "source_sha256": {
            name: _sha256(Path(__file__).resolve().parents[3]/"stiffgwpy_fast"/name)
            for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
        },
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2))


if __name__ == "__main__":
    main()
