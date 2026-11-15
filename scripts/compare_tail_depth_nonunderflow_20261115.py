"""Compare z5/z8 Fast tails against direct Reference at a regular cr=0 point."""

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
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 20.0,
    "kappa10": 0.01,
    "DN_eff": 0.0,
}
FAST_TAILS = (5.0, 8.0)
REFERENCE_TAILS = (8.0, 10.0)
REFERENCE_RTOL = 1e-9
THREADS = 1
ROOT_TOL = 1e-7
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_tail_depth_nonunderflow_20261115.json"
ROOT_TOL = 1e-7


def _self_consistent_fast(
    params: dict,
) -> tuple[LCDM_SG, float, bool, bool, float, str | None]:
    model = LCDM_SG(**params)
    model.gen_expansion()
    started = time.perf_counter()
    result = model.SGWB_iter(
        engine="fast",
        fallback=False,
        tol=ROOT_TOL,
        z_tail=5.0,
        threads=THREADS,
        accuracy_mode="fast",
        freq_grid="goal",
        physical_guard=True,
    )
    runtime = time.perf_counter() - started
    return (
        model,
        runtime,
        result is not None,
        bool(getattr(model, "SGWB_converge", False)),
        (
            float(model.dn_converged_delta)
            if math.isfinite(float(getattr(model, "dn_converged_delta", math.inf)))
            else None
        ),
        model.fast_failure_reason,
    )


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fast_fixed(freqs: np.ndarray, root: float, z_tail: float) -> dict:
    model = LCDM_SG(**PARAMS)
    model.cosmo_param["DN_eff"] = root
    model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid

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
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    runtime = time.perf_counter() - started
    if result is None or not np.array_equal(model.f, freqs):
        raise RuntimeError(f"Fast failed at z_tail={z_tail}: {model.fast_failure_reason}")
    return {
        "runtime_s": runtime,
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "map_residual_at_frozen_root": float(model.DN_gw[-1] - root),
        "log10_OmegaGW": np.asarray(model.log10OmegaGW, dtype=np.float64).tolist(),
    }


def _reference_fixed(freqs: np.ndarray, root: float, z_tail: float, log_amp: float) -> dict:
    unit_params = dict(PARAMS)
    unit_params.update({"r": 1.0, "A_s": 1.0, "n_t": 0.0, "DN_eff": 0.0})
    model = LCDM_SG(**unit_params)
    model.gen_expansion()
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
        log_power = log_amp + PARAMS["n_t"] * (
            math.log(10.0) * freqs - math.log(gp.f_piv)
        )
        log_integrand = np.log(base) + log_power
        shift = float(np.max(log_integrand))
        normalized = np.exp(log_integrand - shift)
        g2_scaled, qerr_scaled, ierr_scaled = REF.integrate_spectrum(
            freqs, normalized, np.zeros_like(normalized)
        )
    runtime = time.perf_counter() - started
    log_spec = (np.log(base) + log_power) / math.log(10.0)
    omega_nu = gp.Omega_nh2 / float(model.derived_param["h"]) ** 2
    mapped_dn = gp.Neff0 * math.exp(shift + math.log(g2_scaled)) / omega_nu
    return {
        "runtime_s": runtime,
        "mapped_DN_gw": float(mapped_dn),
        "quadrature_error_relative": (
            None if qerr_scaled is None else float(qerr_scaled / g2_scaled)
        ),
        "interpolation_error_relative": (
            None if ierr_scaled is None else float(ierr_scaled / g2_scaled)
        ),
        "warnings": [str(item.message) for item in caught],
        "log10_OmegaGW": np.asarray(log_spec, dtype=np.float64).tolist(),
    }


def main() -> int:
    rejected_params = dict(PARAMS)
    rejected_params["n_t"] = 0.5
    (
        rejected_model,
        rejected_runtime,
        rejected_success,
        rejected_converged,
        rejected_delta,
        rejected_failure,
    ) = _self_consistent_fast(rejected_params)
    rejected_root_trial = {
        "status": "rejected_failure" if not rejected_success else (
            "accepted_root_preflight" if rejected_converged and rejected_delta <= ROOT_TOL
            else "rejected_nonconverged_root"
        ),
        "parameters": rejected_params,
        "root_tolerance": ROOT_TOL,
        "root_DN_eff": float(rejected_model.cosmo_param["DN_eff"]),
        "dn_converged_delta": rejected_delta,
        "SGWB_converge": rejected_converged,
        "fast_failure_reason": rejected_failure,
        "frequency_count": int(getattr(rejected_model, "f", np.empty(0)).size),
        "runtime_s": rejected_runtime,
    }

    root_model, root_runtime, root_success, root_converged, root_delta, root_failure = (
        _self_consistent_fast(PARAMS)
    )
    if not root_success or not root_converged or root_delta is None or root_delta > ROOT_TOL:
        raise RuntimeError(
            f"Fast root did not close: success={root_success}, converged={root_converged}, "
            f"delta={root_delta}, failure={root_failure}"
        )
    root = float(root_model.cosmo_param["DN_eff"])
    freqs = np.asarray(root_model.f, dtype=np.float64)
    log_amp = _log_tensor_amplitude(root_model.cosmo_param, root_model.derived_param)
    grid_hash = hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest()

    fast = {
        str(int(z)): _fast_fixed(freqs, root, z)
        for z in FAST_TAILS
    }
    reference = {
        str(int(z)): _reference_fixed(freqs, root, z, log_amp)
        for z in REFERENCE_TAILS
    }
    errors = {}
    for fast_tail in FAST_TAILS:
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
        "status": "complete_nonunderflow_fixed_root_tail_screen",
        "parameters": PARAMS,
        "fast_self_consistent_root": root,
        "fast_self_consistent_runtime_s": root_runtime,
        "dn_converged_delta": float(root_model.dn_converged_delta),
        "root_tolerance": ROOT_TOL,
        "root_converged": root_converged,
        "rejected_root_preflight_trial": rejected_root_trial,
        "fixed_root": root,
        "frequency_count": int(freqs.size),
        "frequency_grid_sha256_f64le": grid_hash,
        "fast_fixed_root": fast,
        "reference_fixed_root": reference,
        "spectrum_error_summaries_dex": errors,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": THREADS,
        },
        "source_sha256": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive.py": _sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "grid_hash_source": "SHA256(little-endian float64 frequency bytes)",
            "diagnostic_script": _sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "self_consistent_root": root,
        "dn_converged_delta": output["dn_converged_delta"],
        "frequency_count": int(freqs.size),
        "map_residuals": {
            z: row["mapped_DN_gw"] - root for z, row in fast.items()
        },
        "spectrum_error_summaries_dex": errors,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
