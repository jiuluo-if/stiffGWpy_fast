"""Audit the Fast terminal moment blend on the frozen 469-node holdout."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 1e-316,
    "n_t": 17.38,
    "cr": 0.0,
    "T_re": 2e3,
    "DN_re": 10.0,
    "kappa10": 1e-2,
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _scale_by_log(values: np.ndarray, log_scale: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    output = np.zeros_like(values)
    nonzero = values != 0.0
    log_abs = np.log(np.abs(values[nonzero])) + log_scale[nonzero]
    signs = np.sign(values[nonzero])
    lo = math.log(np.nextafter(0.0, 1.0))
    hi = math.log(np.finfo(np.float64).max)
    scaled = np.zeros_like(log_abs)
    valid = (log_abs >= lo) & (log_abs <= hi)
    scaled[valid] = signs[valid] * np.exp(log_abs[valid])
    scaled[log_abs > hi] = signs[log_abs > hi] * np.inf
    output[nonzero] = scaled
    return output


def main() -> None:
    data_path = ROOT / "docs" / "underflow_phase_midpoint_reference_arrays_20261107.npz"
    data = np.load(data_path)
    frequencies = np.asarray(data["frequency_log10_hz"], dtype=np.float64)
    support = frequencies[::2]
    if frequencies.size != 469 or support.size != 235:
        raise RuntimeError("expected the frozen nested 235/469 frequency arrays")
    midpoints = 0.5 * (support[:-1] + support[1:])
    nested = np.unique(np.concatenate((support, midpoints)))[::-1]
    if not np.array_equal(nested, frequencies):
        raise RuntimeError("the stored 469-node grid is not the expected 235-node nesting")

    original_goal_builder = FA.goal_oriented_freqs
    FA.goal_oriented_freqs = lambda model, freq_res=1.0, **kwargs: support.copy()
    try:
        FS.apply_accuracy_mode("fast")
        FS.set_threads(2)
        model = LCDM_SG(**PARAMS)
        result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=1e-6,
            z_tail=5.0,
            freq_res=1.0,
            threads=2,
            accuracy_mode="fast",
            physical_guard=False,
            freq_grid="goal",
            eval_freqs=midpoints,
        )
    finally:
        FA.goal_oriented_freqs = original_goal_builder

    if result is not model or not model.SGWB_converge:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    if not np.array_equal(np.asarray(model.f, dtype=np.float64), frequencies):
        raise RuntimeError("Fast output does not match the frozen 469-node grid")

    log_power = math.log(PARAMS["r"]) + math.log(PARAMS["A_s"])
    log_power += PARAMS["n_t"] * (
        math.log(10.0) * frequencies - math.log(gp.f_piv)
    )
    reference_g_unit = np.asarray(data["reference_Ogw_unit"] - data["reference_Oj_unit"])
    if np.any(reference_g_unit <= 0.0):
        raise RuntimeError("reference G must be positive for log-spectrum comparison")
    reference_log_spectrum = np.log10(reference_g_unit) + log_power / math.log(10.0)
    spectrum_error = np.abs(np.asarray(model.log10OmegaGW) - reference_log_spectrum)

    reference_oj = _scale_by_log(data["reference_Oj_unit"], log_power)
    reference_op = _scale_by_log(data["reference_Opgw_unit"], log_power)
    fast_oj = np.asarray(model.Oj_today)
    fast_op = np.asarray(model.Opgw_today)
    oj_l1 = float(np.sum(np.abs(fast_oj - reference_oj)) / np.sum(np.abs(reference_oj)))
    op_l1 = float(np.sum(np.abs(fast_op - reference_op)) / np.sum(np.abs(reference_op)))
    fixed_root = float(data["fixed_DN_eff"][0])
    payload = {
        "schema_version": 1,
        "experiment": "underflow Fast terminal moment contract on frozen nested-frequency holdout",
        "fast_source_sha256": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "freq_source_sha256": _sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
        "input_npz_sha256": _sha(data_path),
        "params": PARAMS,
        "frequency_count": int(frequencies.size),
        "frequency_sha256": hashlib.sha256(frequencies.astype("<f8").tobytes()).hexdigest(),
        "fast_root": float(model.cosmo_param["DN_eff"]),
        "saved_reference_root": fixed_root,
        "root_delta_from_saved_reference": float(model.cosmo_param["DN_eff"] - fixed_root),
        "max_spectrum_error_dex_against_saved_deep_reference": float(spectrum_error.max()),
        "p95_spectrum_error_dex_against_saved_deep_reference": float(np.percentile(spectrum_error, 95)),
        "worst_frequency_log10_hz": float(frequencies[int(np.argmax(spectrum_error))]),
        "Oj_relative_L1_against_saved_deep_reference": oj_l1,
        "Opgw_relative_L1_against_saved_deep_reference": op_l1,
        "Opgw_finite": bool(np.isfinite(fast_op).all()),
        "w2_finite": bool(np.isfinite(model.w2).all()),
        "w2_final": float(model.w2[-1]),
        "scope_limit": "Reference moments are at the saved fixed root; w2 has no matching saved Reference field. This checks one interior point, not a continuous region.",
    }
    out = ROOT / "docs" / "underflow_phase_terminal_moment_contract_20261108.json"
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
