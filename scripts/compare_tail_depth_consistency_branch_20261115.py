"""Compare Fast tails on one physical-amplitude cr>0 consistency case."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 1.0,
    "T_re": 2000.0,
    "DN_re": 20.0,
    "kappa10": 0.01,
    "DN_eff": 0.0,
}
ROOT_TOL = 1e-7
REFERENCE_RTOL = 1e-9
THREADS = 1
TAILS = (5.0, 8.0)
REFERENCE_TAILS = (8.0, 10.0)
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_tail_depth_consistency_branch_20261115.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reference_fixed(model: LCDM_SG, freqs: np.ndarray, root: float, z_tail: float) -> dict:
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        Ogw, Oj, Opgw, used_tail = REF.spectrum_reference(
            model,
            freqs,
            root,
            z_tail=z_tail,
            rtol=REFERENCE_RTOL,
            workers=THREADS,
            backend="cartesian",
        )
        base = np.where(used_tail, 3.0 * Opgw, Ogw - Oj)
        log10_spectrum = np.log(base) / math.log(10.0)
        g2, qerr, ierr = REF.integrate_spectrum(
            freqs, base, np.zeros_like(base)
        )
    runtime = time.perf_counter() - started
    omega_nu = gp.Omega_nh2 / float(model.derived_param["h"]) ** 2
    mapped_dn = gp.Neff0 * float(g2) / omega_nu
    return {
        "runtime_s": runtime,
        "mapped_DN_gw": float(mapped_dn),
        "log10_OmegaGW": np.asarray(log10_spectrum, dtype=np.float64).tolist(),
        "quadrature_error_relative": (
            None if qerr is None else float(qerr / g2)
        ),
        "interpolation_error_relative": (
            None if ierr is None else float(ierr / g2)
        ),
        "warnings": [str(item.message) for item in caught],
    }


def _fast_fixed(freqs: np.ndarray, root: float, z_tail: float) -> dict:
    model = LCDM_SG(**PARAMS)
    model.cosmo_param["DN_eff"] = root
    model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal, original_construct = FA.goal_oriented_freqs, FA.construct_fast_frequency_grid

    def locked_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = root
        return original_gen(*args, **kwargs)

    def fixed_grid(*args, **kwargs):
        return freqs.copy()

    model.gen_expansion = locked_gen
    FA.goal_oriented_freqs = fixed_grid
    FA.construct_fast_frequency_grid = fixed_grid
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=0.99,
            z_tail=z_tail,
            threads=THREADS,
            accuracy_mode="fast",
            freq_grid="goal",
            physical_guard=True,
        )
    finally:
        FA.goal_oriented_freqs, FA.construct_fast_frequency_grid = original_goal, original_construct
    runtime = time.perf_counter() - started
    if result is None or not np.array_equal(model.f, freqs):
        raise RuntimeError(f"Fast failed at z_tail={z_tail}: {model.fast_failure_reason}")
    return {
        "runtime_s": runtime,
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "map_residual_at_frozen_root": float(model.DN_gw[-1] - root),
        "log10_OmegaGW": np.asarray(model.log10OmegaGW, dtype=np.float64).tolist(),
    }


def main() -> int:
    root_model = LCDM_SG(**PARAMS)
    root_model.gen_expansion()
    root_started = time.perf_counter()
    root_result = root_model.SGWB_iter(
        engine="fast",
        fallback=False,
        tol=ROOT_TOL,
        z_tail=5.0,
        threads=THREADS,
        accuracy_mode="fast",
        freq_grid="goal",
        physical_guard=True,
    )
    root_runtime = time.perf_counter() - root_started
    if root_result is None:
        raise RuntimeError(f"Fast root failed: {root_model.fast_failure_reason}")
    converged = bool(getattr(root_model, "SGWB_converge", False))
    root_delta = float(root_model.dn_converged_delta)
    if not converged or root_delta > ROOT_TOL:
        raise RuntimeError(
            f"Fast root did not close: converged={converged}, delta={root_delta}"
        )
    root = float(root_model.cosmo_param["DN_eff"])
    freqs = np.asarray(root_model.f, dtype=np.float64)

    fast = {str(int(z)): _fast_fixed(freqs, root, z) for z in TAILS}
    reference_model = LCDM_SG(**PARAMS)
    reference_model.cosmo_param["DN_eff"] = root
    reference_model.gen_expansion()
    reference_at = float(reference_model.derived_param["A_t"])
    fast_at = float(root_model.derived_param["A_t"])
    if not math.isclose(reference_at, fast_at, rel_tol=0.0, abs_tol=0.0):
        raise AssertionError("Reference consistency branch changed physical A_t")
    reference = {
        str(int(z)): _reference_fixed(reference_model, freqs, root, z)
        for z in REFERENCE_TAILS
    }
    errors = {}
    for fast_tail in TAILS:
        for ref_tail in REFERENCE_TAILS:
            key = f"fast_z{int(fast_tail)}_vs_reference_z{int(ref_tail)}"
            delta = np.abs(
                np.asarray(fast[str(int(fast_tail))]["log10_OmegaGW"])
                - np.asarray(reference[str(int(ref_tail))]["log10_OmegaGW"])
            )
            errors[key] = {
                "max_abs_dex": float(np.max(delta)),
                "p95_abs_dex": float(np.quantile(delta, 0.95)),
            }
    output = {
        "schema_version": 1,
        "status": "complete_cr_positive_physical_amplitude_tail_screen",
        "parameters": PARAMS,
        "derived_nt": float(root_model.derived_param["nt"]),
        "physical_A_t": float(root_model.derived_param["A_t"]),
        "reference_physical_A_t": reference_at,
        "background_derived": {
            "N_inf": float(root_model.derived_param["N_inf"]),
            "N_re": float(root_model.derived_param["N_re"]),
            "N_re_abs": float(
                root_model.derived_param["N_inf"] - root_model.derived_param["N_re"]
            ),
            "input_DN_re": float(root_model.cosmo_param["DN_re"]),
        },
        "fast_root_DN_eff": root,
        "fast_root_converged": converged,
        "fast_root_delta": root_delta,
        "fast_root_runtime_s": root_runtime,
        "frequency_count": int(freqs.size),
        "frequency_grid_sha256_f64le": hashlib.sha256(
            np.asarray(freqs, dtype="<f8").tobytes()
        ).hexdigest(),
        "fast_fixed_root": fast,
        "reference_fixed_root_physical_amplitude": reference,
        "spectrum_error_summaries_dex": errors,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": THREADS,
        },
        "source_sha256": {
            "LCDM_stiff_Neff.py": _sha(ROOT / "stiffgwpy_fast" / "LCDM_stiff_Neff.py"),
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "diagnostic_script": _sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "root_DN_eff": root,
        "derived_nt": output["derived_nt"],
        "physical_A_t": output["physical_A_t"],
        "frequency_count": int(freqs.size),
        "fast_map_residuals": {
            z: row["map_residual_at_frozen_root"] for z, row in fast.items()
        },
        "spectrum_error_summaries_dex": errors,
        "reference_map_residuals": {
            z: row["mapped_DN_gw"] - root for z, row in reference.items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
