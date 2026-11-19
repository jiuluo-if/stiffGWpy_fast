"""Persist direct z8 Reference spectra for the exact-quarter CR1 patch."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HERE = Path(__file__).resolve().parent
BASE_PARAMS = {
    "H0": 67.32117, "DN_eff": 0.0, "A_s": 2.100549e-9,
    "r": 0.01, "n_t": 0.0, "cr": 1.0, "T_re": 2000.0,
    "DN_re": 10.0, "kappa10": 0.01,
}
TOTALS = {1300: (700_000.0, 600_000.0),
          1350: (675_000.0, 675_000.0),
          1400: (700_000.0, 700_000.0)}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if len(sys.argv) != 2 or int(sys.argv[1]) not in TOTALS:
        raise SystemExit("usage: capture_exact_quarter_reference_spectra.py 1300|1350|1400")
    label = int(sys.argv[1])
    ob, oc = TOTALS[label]
    params = dict(BASE_PARAMS, Omega_bh2=ob, Omega_ch2=oc)
    model = LCDM_SG(**params)
    if SGWB_iter_fast(model, config=resolve_config("fast", threads=2),
                      tol=1e-6, physical_guard=True) is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    frequencies = np.asarray(model.f, dtype=np.float64)
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    result = run_reference(
        model, dn_eff=float(model.DN_gw[-1]), freq_subset=frequencies,
        z_tail=8.0, rtol=1e-9, self_consistent=False, workers=2,
    )
    ref_freqs = np.asarray(result["freqs"], dtype=np.float64)
    if not np.array_equal(frequencies, ref_freqs):
        raise RuntimeError("Reference frequencies differ from exact-quarter Fast grid")
    ref_log = np.asarray(result["log10OmegaGW"], dtype=np.float64)
    ref_ogw = np.asarray(result["Ogw"], dtype=np.float64)
    ref_oj = np.asarray(result["Oj"], dtype=np.float64)
    valid = (ref_ogw - ref_oj > 0.0) & np.isfinite(fast_log) & np.isfinite(ref_log)
    payload = {
        "case": f"cr1_total_{label}k", "params": params,
        "reference_code_modified": False, "reference_tail_z": 8.0,
        "fast_root": float(model.DN_gw[-1]), "reference_map_at_fast_root": float(result["DN_gw"]),
        "same_frequency_grid": True,
        "spectrum_max_abs_dex": float(np.max(np.abs(fast_log[valid] - ref_log[valid]))),
        "valid_spectrum_nodes": int(valid.sum()),
        "freqs_log10_hz": frequencies.tolist(),
        "fast_log10_OmegaGW": fast_log.tolist(),
        "reference_log10_OmegaGW": ref_log.tolist(),
        "reference_Ogw": ref_ogw.tolist(), "reference_Oj": ref_oj.tolist(),
        "source_sha256": {
            name: _sha(ROOT / "stiffgwpy_fast" / name)
            for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
        },
    }
    output = HERE / f"reference_spectrum_exact_quarter_{label}k.json"
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items()
                      if key not in {"freqs_log10_hz", "fast_log10_OmegaGW",
                                     "reference_log10_OmegaGW", "reference_Ogw", "reference_Oj"}},
                     indent=2))


if __name__ == "__main__":
    main()
