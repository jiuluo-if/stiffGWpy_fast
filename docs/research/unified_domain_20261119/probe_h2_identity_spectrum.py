"""One-sample diagnostic for the rejected H2-identity primitive candidate.

This runner changes no production file. It replaces only the Fast primitive
builder at runtime, then compares one fixed-point map and its exact frequency
grid against the unchanged Reference implementation. It is diagnostic only;
the candidate failed the current root and spectrum gates.
"""

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

exact_background = importlib.import_module("stiffgwpy_fast.exact_background")
fast_sgwb = importlib.import_module("stiffgwpy_fast.fast_sgwb")
run_reference = importlib.import_module("stiffgwpy_fast.reference").run_reference
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG
SGWB_iter_fast = fast_sgwb.SGWB_iter_fast
resolve_config = fast_sgwb.resolve_config


REFERENCE_Z_TAIL = 8.0
OUTPUT = Path(__file__).with_name("probe_h2_identity_wkb_adaptive_z8.json")
PARAMS = {
    "Omega_bh2": 1.0e5,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _h2_identity_builder(model, Nv, DN_eff, sigma_nodes=None, *, relative_scale=False):
    original = _h2_identity_builder.original
    old = original(model, Nv, DN_eff, sigma_nodes=sigma_nodes, relative_scale=True)
    nodes = np.asarray(Nv, dtype=np.float64)
    midpoints = 0.5 * (nodes[:-1] + nodes[1:])
    n0 = float(nodes[0])
    n_re = float(model.derived_param["N_inf"] - model.derived_param["N_re"])
    anchor = max(n0, n_re)
    h2_anchor = float(exact_background.H2_vec(np.array([anchor]), model, DN_eff)[0])
    post_nodes = nodes >= n_re
    post_mids = midpoints >= n_re
    h2_nodes = exact_background.H2_vec(nodes[post_nodes], model, DN_eff)
    h2_mids = exact_background.H2_vec(midpoints[post_mids], model, DN_eff)

    F_nodes = nodes - n0
    F_mids = midpoints - n0
    F_nodes[post_nodes] -= (np.log(h2_nodes) - math.log(h2_anchor)) / 3.0
    F_mids[post_mids] -= (np.log(h2_mids) - math.log(h2_anchor)) / 3.0
    phi = 1.5 * F_nodes - nodes + n0
    phi_mid = 1.5 * F_mids - midpoints + n0
    psi = 3.0 * F_nodes - 4.0 * nodes

    log_max = math.log(np.finfo(np.float64).max)
    log_min = math.log(np.nextafter(0.0, 1.0))
    relative = bool(
        np.any(psi > log_max)
        or np.any(psi < log_min)
        or np.any(-0.5 * psi > log_max)
        or np.any(-0.5 * psi < log_min)
    )
    if relative:
        s2 = np.ones_like(nodes)
        s2inv = np.ones_like(nodes)
    else:
        s2 = np.exp(psi)
        s2inv = np.exp(-0.5 * psi)
    model._fast_relative_background_required = relative
    phi_re = 0.5 * (n_re - n0) if n0 < n_re <= nodes[-1] else 0.0
    return phi, phi_mid, s2, s2inv, old[4], old[5], phi_re


_h2_identity_builder.original = exact_background.fast_phi_s2_split


def main() -> None:
    model = LCDM_SG(**PARAMS)
    original = exact_background.fast_phi_s2_split
    config = resolve_config("fast", threads=2)
    exact_background.fast_phi_s2_split = _h2_identity_builder
    try:
        start = time.perf_counter()
        SGWB_iter_fast(
            model,
            config=config,
            tol=1e-6,
        )
        fast_seconds = time.perf_counter() - start
    finally:
        exact_background.fast_phi_s2_split = original

    freqs = np.asarray(model.f, dtype=np.float64)
    total_dn_eff = float(model.cosmo_param["DN_eff"])
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    start = time.perf_counter()
    reference = run_reference(
        model,
        dn_eff=total_dn_eff,
        freq_subset=freqs,
        z_tail=REFERENCE_Z_TAIL,
        rtol=1e-9,
        self_consistent=False,
        workers=2,
    )
    reference_seconds = time.perf_counter() - start

    ref_log = np.asarray(reference["log10OmegaGW"], dtype=np.float64)
    positive = np.asarray(reference["Ogw"]) - np.asarray(reference["Oj"]) > 0.0
    valid = np.isfinite(fast_log) & np.isfinite(ref_log) & positive
    error = fast_log[valid] - ref_log[valid]
    source_paths = [
        ROOT / "stiffgwpy_fast" / "exact_background.py",
        ROOT / "stiffgwpy_fast" / "fast_sgwb.py",
        ROOT / "stiffgwpy_fast" / "reference.py",
    ]
    result = {
        "candidate_status": "exploratory_wkb_adaptive_h2_identity_candidate",
        "params": PARAMS,
        "fast_config": {
            "accuracy_mode": "fast",
            "h": 0.005,
            "z_tail": float(config.z_tail),
            "phase_max": 0.25,
            "threads": 2,
            "frequency_grid": "goal",
            "n_freq": int(freqs.size),
        },
        "fast": {
            "root": float(model.DN_gw[-1]),
            "total_DN_eff": total_dn_eff,
            "seconds": fast_seconds,
            "failure_reason": model.fast_failure_reason,
            "nonadiabatic_handoff_selected": int(
                getattr(model, "_fast_nonadiabatic_handoff_selected", 0)
            ),
        },
        "reference_fixed_root": {
            "DN_gw_map": float(reference["DN_gw"]),
            "fixed_point_residual": float(reference["DN_gw"] - total_dn_eff),
            "seconds": reference_seconds,
            "z_tail": REFERENCE_Z_TAIL,
            "rtol": 1e-9,
        },
        "spectrum_error_dex": {
            "valid_count": int(valid.sum()),
            "max_abs": float(np.max(np.abs(error))),
            "p95_abs": float(np.quantile(np.abs(error), 0.95)),
            "median_abs": float(np.median(np.abs(error))),
        },
        "raw": {
            "freqs_log10_hz": freqs.tolist(),
            "fast_log10_OmegaGW": fast_log.tolist(),
            "reference_log10_OmegaGW": ref_log.tolist(),
            "valid_positive_reference_mask": positive.tolist(),
            "signed_spectrum_error_dex": (fast_log - ref_log).tolist(),
        },
        "source_sha256": {str(path.relative_to(ROOT)): _sha256(path) for path in source_paths},
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("candidate_status", "fast", "reference_fixed_root", "spectrum_error_dex", "source_sha256")}, indent=2))
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
