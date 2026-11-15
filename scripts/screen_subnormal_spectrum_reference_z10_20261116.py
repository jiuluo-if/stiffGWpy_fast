"""Compare the Fast underflow spectrum with direct Cartesian Reference z10."""

from __future__ import annotations

import hashlib
import json
import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude, construct_fast_frequency_grid  # noqa: E402
from stiffgwpy_fast.reference import integrate_spectrum, spectrum_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BASE = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "n_t": 0.0,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 1e-2,
}
R = 1e-315
OUTPUT = ROOT / "docs" / "research" / "results" / "subnormal_spectrum_reference_z10_20261116.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    physical = LCDM_SG(**{**BASE, "r": R})
    physical.gen_expansion()
    log_amp = _log_tensor_amplitude(physical.cosmo_param, physical.derived_param)
    grid_model = LCDM_SG(**{**BASE, "r": R})
    grid_model.gen_expansion()
    freqs = np.asarray(
        construct_fast_frequency_grid(grid_model, freq_res=2.0, log_amplitude=log_amp),
        dtype=np.float64,
    )

    reference_model = LCDM_SG(**{**BASE, "r": 1.0, "A_s": 1.0})
    reference_model.gen_expansion()
    ref_started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        Ogw, Oj, Opgw, used_tail = spectrum_reference(
            reference_model, freqs, 0.0, z_tail=10.0, rtol=1e-9,
            workers=1, backend="cartesian",
        )
        base_signal = np.where(used_tail, 3.0 * Opgw, Ogw - Oj)
        reference_log10 = np.log10(base_signal) + log_amp / math.log(10.0)
        reference_runtime = time.perf_counter() - ref_started

    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    try:
        FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
        FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
        fast_model = LCDM_SG(**{**BASE, "r": R})
        fast_started = time.perf_counter()
        result = fast_model.SGWB_iter(
            engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
            threads=1, accuracy_mode="fast", freq_grid="goal",
            physical_guard=True,
        )
        fast_runtime = time.perf_counter() - fast_started
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct

    fast_log10 = np.asarray(fast_model.log10OmegaGW, dtype=np.float64)
    if result is not fast_model or fast_log10.shape != reference_log10.shape:
        raise RuntimeError("Fast failed or returned a spectrum on a different frequency grid")
    delta = fast_log10 - reference_log10
    idx = int(np.argmax(np.abs(delta)))
    h = float(reference_model.derived_param["h"])
    log_map_shift = float(np.max(np.log(base_signal) + log_amp))
    normalized = np.exp(np.log(base_signal) + log_amp - log_map_shift)
    g2, qerr, ierr = integrate_spectrum(freqs, normalized, np.zeros_like(normalized))
    log_reference_map = (
        math.log(gp.Neff0 / (gp.Omega_nh2 / h**2))
        + log_map_shift + math.log(g2)
    )
    output = {
        "schema_version": 1,
        "status": "complete_same_grid_fast_reference_z10_spectrum_screen",
        "scope": "Pointwise Fast z5 versus direct unit-transfer Cartesian Reference z10 for cr=0,n_t=0,r=1e-315 on the same 457-node grid. Reference amplitude is restored in log10 space.",
        "parameters": {**BASE, "r": R},
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
        "fast": {
            "z_tail": 5.0,
            "tol": 1e-7,
            "threads": 1,
            "runtime_s": fast_runtime,
            "root_DN_eff": float(fast_model.cosmo_param["DN_eff"]),
            "mapped_DN_gw": float(fast_model.DN_gw[-1]),
            "failure_reason": fast_model.fast_failure_reason,
        },
        "reference": {
            "z_tail": 10.0,
            "rtol": 1e-9,
            "workers": 1,
            "runtime_s": reference_runtime,
            "log_DN_gw": log_reference_map,
            "DN_gw_float64": math.exp(log_reference_map),
            "used_tail_count": int(np.count_nonzero(used_tail)),
            "quadrature_error_relative": None if qerr is None else float(qerr / g2),
            "interpolation_error_relative": None if ierr is None else float(ierr / g2),
            "warnings": [str(item.message) for item in caught],
        },
        "spectrum_error": {
            "max_abs_dex": float(np.max(np.abs(delta))),
            "p95_abs_dex": float(np.percentile(np.abs(delta), 95.0)),
            "max_error_signed_dex": float(delta[idx]),
            "max_error_log10f": float(freqs[idx]),
            "fast_log10OmegaGW": fast_log10.tolist(),
            "reference_log10OmegaGW": reference_log10.tolist(),
            "passes_1e-3_dex": bool(np.max(np.abs(delta)) <= 1e-3),
        },
        "source_sha256": {
            "fast_sgwb.py": sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive.py": sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "reference.py": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "script": sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "frequency_count": freqs.size,
        "spectrum_error": output["spectrum_error"]["max_abs_dex"],
        "p95": output["spectrum_error"]["p95_abs_dex"],
        "worst_log10f": output["spectrum_error"]["max_error_log10f"],
        "fast_map": output["fast"]["mapped_DN_gw"],
        "reference_map": output["reference"]["DN_gw_float64"],
        "reference_runtime_s": reference_runtime,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
