"""z10 same-grid closure for the independent amplitude-tilt Sobol holdout."""

from __future__ import annotations

import hashlib
import importlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
fast = importlib.import_module("stiffgwpy_fast.fast_sgwb")
run_reference = importlib.import_module("stiffgwpy_fast.reference").run_reference
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "validate_cr0_blue_tilt_sobol_z10.json"
PARAMS = {
    "Omega_bh2": 1.0e6, "Omega_ch2": 0.1201075, "H0": 67.32117,
    "DN_eff": 0.0, "A_s": 2.100549e-9, "r": 1.676304587225033e-9,
    "n_t": 0.4469903271645308, "cr": 0.0, "T_re": 2000.0,
    "DN_re": 10.0, "kappa10": 0.01,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    model = LCDM_SG(**PARAMS)
    config = fast.resolve_config("fast", threads=2)
    fast_started = time.perf_counter()
    solved = fast.SGWB_iter_fast(
        model, config=config, tol=1.0e-6, physical_guard=True
    )
    fast_seconds = time.perf_counter() - fast_started
    if solved is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    freqs = np.asarray(model.f, dtype=np.float64)
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    fast_root = float(model.DN_gw[-1])
    root = float(model.cosmo_param["DN_eff"])
    maps = []
    first = None
    closure_residual = math.inf
    for _ in range(4):
        started = time.perf_counter()
        ref = run_reference(
            model, dn_eff=root, freq_subset=freqs, z_tail=10.0,
            rtol=1.0e-9, self_consistent=False, workers=2,
        )
        seconds = time.perf_counter() - started
        if first is None:
            first = ref
        mapped = float(ref["DN_gw"])
        maps.append({"input_root": root, "mapped_root": mapped, "seconds": seconds})
        closure_residual = abs(mapped-root)/max(abs(mapped),1e-300)
        root = mapped
        if closure_residual < 2.0e-5:
            break

    ref_log = np.asarray(first["log10OmegaGW"], dtype=np.float64)
    ref_linear = np.asarray(first["Ogw"]) - np.asarray(first["Oj"])
    valid = np.isfinite(fast_log) & np.isfinite(ref_log) & (ref_linear > 0.0)
    source = {
        name: _sha256(ROOT / "stiffgwpy_fast" / name)
        for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
    }
    payload = {
        "status": "complete_z10_root_closure",
        "reference_code_modified": False,
        "params": PARAMS,
        "fast_config": {"threads": 2, "h": float(config.h),
                        "z_tail": float(config.z_tail), "tol": 1.0e-6,
                        "frequency_grid": "goal"},
        "candidate_source_sha256": source,
        "fast_seconds": fast_seconds,
        "fast_root": fast_root,
        "closed_reference_root": root,
        "root_relative_error": abs(root-fast_root)/max(abs(root),1e-300),
        "reference_closure_residual": closure_residual,
        "map_history": maps,
        "spectrum_valid_count": int(valid.sum()),
        "spectrum_max_abs_dex": float(np.max(np.abs(fast_log[valid]-ref_log[valid]))),
        "spectrum_p95_abs_dex": float(np.quantile(np.abs(fast_log[valid]-ref_log[valid]),0.95)),
        "freqs_log10_hz": freqs.tolist(),
        "fast_log10_OmegaGW": fast_log.tolist(),
        "z10_log10_OmegaGW_at_fast_root": ref_log.tolist(),
        "z10_Ogw_at_fast_root": np.asarray(first["Ogw"]).tolist(),
        "z10_Oj_at_fast_root": np.asarray(first["Oj"]).tolist(),
        "valid_positive_reference_mask": valid.tolist(),
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:payload[k] for k in (
        "status","fast_root","closed_reference_root","root_relative_error",
        "reference_closure_residual","spectrum_max_abs_dex","spectrum_p95_abs_dex",
        "map_history")},indent=2))


if __name__ == "__main__":
    main()
