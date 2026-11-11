"""Check Fast's nonzero-underflow amplitude branch against scaled Reference."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
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


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _scale_by_log(values: np.ndarray, log_scale: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    result = np.zeros_like(values)
    nonzero = values != 0.0
    log_abs = np.log(np.abs(values[nonzero])) + log_scale[nonzero]
    signs = np.sign(values[nonzero])
    max_log = math.log(np.finfo(np.float64).max)
    min_log = math.log(np.nextafter(0.0, 1.0))
    scaled = np.zeros(log_abs.shape, dtype=np.float64)
    representable = (log_abs >= min_log) & (log_abs <= max_log)
    with np.errstate(over="ignore", under="ignore", invalid="raise"):
        scaled[representable] = signs[representable] * np.exp(log_abs[representable])
    scaled[log_abs > max_log] = signs[log_abs > max_log] * np.inf
    result[nonzero] = scaled
    return result


def _weighted_transfer_pchip_integral(
    freqs: np.ndarray, transfer: np.ndarray, log_at: float, n_t: float
) -> float:
    """Integrate a cubic PCHIP transfer multiplied by the exact blue power."""
    x = np.asarray(freqs, dtype=np.float64)
    y = np.asarray(transfer, dtype=np.float64)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    spline = PchipInterpolator(x, y)
    nodes, weights = np.polynomial.legendre.leggauss(16)
    total = 0.0
    for left, right in zip(x[:-1], x[1:]):
        half = 0.5 * (right - left)
        mid = 0.5 * (right + left)
        sample_x = mid + half * nodes
        sample_y = np.asarray(spline(sample_x), dtype=np.float64)
        log_power = log_at + n_t * (
            math.log(10.0) * sample_x - math.log(gp.f_piv)
        )
        positive = sample_y > 0.0
        if np.any(positive):
            log_integrand = np.log(sample_y[positive]) + log_power[positive]
            total += half * float(np.dot(weights[positive], np.exp(log_integrand)))
    return total


def run(
    z_tail: float, rtol: float, root_iterations: int,
    fast_z_tail: float, fast_freq_grid: str, overrides: dict | None = None,
) -> dict:
    params = dict(PARAMS)
    if overrides:
        params.update(overrides)
    candidate = LCDM_SG(**params)
    if candidate.derived_param["A_t"] != 0.0:
        raise AssertionError("test input no longer underflows in float64")
    started = time.perf_counter()
    candidate.SGWB_iter(
        engine="fast",
        fallback=False,
        tol=1e-6,
        z_tail=fast_z_tail,
        freq_res=1.0,
        threads=1,
        accuracy_mode="fast",
        physical_guard=False,
        freq_grid=fast_freq_grid,
    )
    fast_runtime = time.perf_counter() - started
    if not candidate.SGWB_converge:
        from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude

        raise RuntimeError(
            f"Fast failed: {candidate.fast_failure_reason}; "
            f"last_error={candidate.last_fast_error}; "
            f"A_t={candidate.derived_param['A_t']}; "
            f"log_A_t={_log_tensor_amplitude(candidate.cosmo_param, candidate.derived_param)}"
        )

    freqs = np.asarray(candidate.f, dtype=np.float64)
    root = float(candidate.cosmo_param["DN_eff"])
    log_at = math.log(params["r"]) + math.log(params["A_s"])
    log_pt = log_at + params["n_t"] * (
        math.log(10.0) * freqs - math.log(gp.f_piv)
    )

    # In the free-tilt branch the background and transfer do not depend on
    # A_t or n_t. Use an ordinary unit-amplitude, n_t=0 Reference model and
    # apply the exact primordial-power factor afterward. Reference source and
    # defaults remain untouched; this is a scaled oracle for the underflowing
    # target whose direct Reference path cannot represent A_t.
    reference_params = dict(params, r=1.0, A_s=1.0, n_t=0.0)
    oracle = LCDM_SG(**reference_params)
    started = time.perf_counter()
    ref_ogw, ref_oj, ref_opgw, used_tail = REF.spectrum_reference(
        oracle, freqs, DN_eff=root, z_tail=z_tail, rtol=rtol,
        workers=1, backend="cartesian",
    )
    reference_runtime = time.perf_counter() - started
    ref_integrand_unit = np.asarray(ref_ogw - ref_oj, dtype=np.float64)
    ref_log_omega = np.full(freqs.shape, -np.inf, dtype=np.float64)
    positive = ref_integrand_unit > 0.0
    ref_log_omega[positive] = (
        np.log10(ref_integrand_unit[positive]) + log_pt[positive] / math.log(10.0)
    )
    scaled_ogw = _scale_by_log(np.asarray(ref_ogw), log_pt)
    scaled_oj = _scale_by_log(np.asarray(ref_oj), log_pt)
    g2, quadrature_error, interpolation_error = REF.integrate_spectrum(
        freqs, scaled_ogw, scaled_oj
    )
    omega_nu = gp.Omega_nh2 / oracle.derived_param["h"] ** 2
    ref_dn_map = gp.Neff0 * g2 / omega_nu
    fast_dn_map = root - float(params["DN_eff"])
    fast_log_transfer = candidate.log10OmegaGW - log_pt / math.log(10.0)
    fast_transfer = np.power(10.0, fast_log_transfer)
    ref_transfer = np.power(10.0, ref_log_omega - log_pt / math.log(10.0))
    fast_weighted_pchip_g2 = math.log(10.0) * _weighted_transfer_pchip_integral(
        freqs, fast_transfer, log_at, params["n_t"]
    )
    ref_weighted_pchip_g2 = math.log(10.0) * _weighted_transfer_pchip_integral(
        freqs, ref_transfer, log_at, params["n_t"]
    )
    fast_weighted_pchip_map = gp.Neff0 * fast_weighted_pchip_g2 / omega_nu
    ref_weighted_pchip_map = gp.Neff0 * ref_weighted_pchip_g2 / omega_nu
    spectrum_mask = np.isfinite(ref_log_omega) & np.isfinite(candidate.log10OmegaGW)
    spectrum_delta = np.abs(
        candidate.log10OmegaGW[spectrum_mask] - ref_log_omega[spectrum_mask]
    )
    valid_indices = np.flatnonzero(spectrum_mask)
    worst_position = int(np.argmax(spectrum_delta))
    worst_index = int(valid_indices[worst_position])

    root_history = [{"input_DN_eff": root, "mapped_DN_gw": float(ref_dn_map)}]
    reference_root = root
    mapped = float(ref_dn_map)
    for _ in range(max(0, root_iterations - 1)):
        reference_root = mapped
        started = time.perf_counter()
        ref_ogw, ref_oj, _, _ = REF.spectrum_reference(
            oracle, freqs, DN_eff=reference_root, z_tail=z_tail, rtol=rtol,
            workers=1, backend="cartesian",
        )
        reference_runtime += time.perf_counter() - started
        scaled_ogw = _scale_by_log(np.asarray(ref_ogw), log_pt)
        scaled_oj = _scale_by_log(np.asarray(ref_oj), log_pt)
        g2, _, _ = REF.integrate_spectrum(freqs, scaled_ogw, scaled_oj)
        mapped = float(gp.Neff0 * g2 / omega_nu)
        root_history.append(
            {"input_DN_eff": reference_root, "mapped_DN_gw": mapped}
        )
    reference_root = mapped

    source_paths = [
        ROOT / "stiffgwpy_fast" / "fast_sgwb.py",
        ROOT / "stiffgwpy_fast" / "freq_adaptive.py",
        ROOT / "stiffgwpy_fast" / "exact_background.py",
        ROOT / "stiffgwpy_fast" / "reference.py",
        ROOT / "stiffgwpy_fast" / "stiff_SGWB.py",
        ROOT / "stiffgwpy_fast" / "LCDM_stiff_Neff.py",
        ROOT / "stiffgwpy_fast" / "global_param.py",
        Path(__file__),
    ]
    return {
        "schema_version": 1,
        "experiment": "Fast nonzero tensor amplitude underflow with unit-amplitude scaled Reference",
        "status": "complete",
        "reference_tail": float(z_tail),
        "fast_tail": float(fast_z_tail),
        "fast_frequency_grid": fast_freq_grid,
        "reference_rtol": float(rtol),
        "reference_root_iterations": int(root_iterations),
        "params": params,
        "float_A_t": float(candidate.derived_param["A_t"]),
        "log10_A_t_from_inputs": float(log_at / math.log(10.0)),
        "fast": {
            "converged": bool(candidate.SGWB_converge),
            "failure_reason": candidate.fast_failure_reason,
            "root_DN_eff": root,
            "DN_gw": fast_dn_map,
            "runtime_s": fast_runtime,
            "frequency_count": int(freqs.size),
            "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
            "log_spectrum_sha256": hashlib.sha256(
                np.asarray(candidate.log10OmegaGW, dtype="<f8").tobytes()
            ).hexdigest(),
            "linear_spectrum_sha256": hashlib.sha256(
                np.asarray(candidate.Ogw_today - candidate.Oj_today, dtype="<f8").tobytes()
            ).hexdigest(),
            "finite_log_spectrum_count": int(np.isfinite(candidate.log10OmegaGW).sum()),
            "linear_spectrum_zero_count": int(np.count_nonzero(candidate.Ogw_today - candidate.Oj_today == 0.0)),
        },
        "scaled_reference": {
            "backend": "cartesian unit A_t=1, n_t=0; exact log-power scaling by target inputs",
            "map_DN_gw_at_fast_root": float(ref_dn_map),
            "map_relative_residual_at_fast_root": float((ref_dn_map - fast_dn_map) / fast_dn_map),
            "weighted_transfer_pchip_map_at_fast_root": float(ref_weighted_pchip_map),
            "weighted_transfer_pchip_relative_root_residual": float(
                (ref_weighted_pchip_map - fast_dn_map) / fast_dn_map
            ),
            "fast_weighted_transfer_pchip_map": float(fast_weighted_pchip_map),
            "reference_weighted_transfer_pchip_map": float(ref_weighted_pchip_map),
            "root_iteration_history": root_history,
            "estimated_self_consistent_root_DN_gw": float(reference_root),
            "relative_root_difference_vs_fast": float(
                (reference_root - fast_dn_map) / fast_dn_map
            ),
            "runtime_s": reference_runtime,
            "spectrum_comparison_count": int(spectrum_mask.sum()),
            "spectrum_max_abs_delta_dex": float(np.max(spectrum_delta)),
            "spectrum_p95_abs_delta_dex": float(np.quantile(spectrum_delta, 0.95)),
            "quadrature_error": float(quadrature_error),
            "interpolation_error": float(interpolation_error),
            "used_tail_count": int(np.count_nonzero(used_tail)),
            "worst_spectrum_frequency_log10_hz": float(freqs[worst_index]),
            "worst_spectrum_delta_dex": float(spectrum_delta[worst_position]),
        },
        "frequency_rows": [
            {
                "log10_hz": float(freqs[i]),
                "log10_P_t": float(log_pt[i] / math.log(10.0)),
                "fast_log10_OmegaGW": float(candidate.log10OmegaGW[i]),
                "reference_log10_OmegaGW": float(ref_log_omega[i]),
                "delta_dex": float(candidate.log10OmegaGW[i] - ref_log_omega[i]),
                "reference_used_tail": bool(used_tail[i]),
            }
            for i in range(freqs.size)
        ],
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "source_sha256": {path.name: _sha256(path) for path in source_paths},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--z-tail", type=float, default=5.0)
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--root-iterations", type=int, default=1)
    parser.add_argument("--fast-z-tail", type=float, default=5.0)
    parser.add_argument("--fast-freq-grid", choices=("goal", "construct"), default="goal")
    parser.add_argument("--r", type=float)
    parser.add_argument("--a-s", type=float)
    parser.add_argument("--n-t", type=float)
    parser.add_argument("--t-re", type=float)
    parser.add_argument("--kappa10", type=float)
    parser.add_argument("--dn-re", type=float)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    overrides = {
        key: value for key, value in {
            "r": args.r,
            "A_s": args.a_s,
            "n_t": args.n_t,
            "T_re": args.t_re,
            "kappa10": args.kappa10,
            "DN_re": args.dn_re,
        }.items() if value is not None
    }
    result = run(
        args.z_tail, args.rtol, args.root_iterations,
        args.fast_z_tail, args.fast_freq_grid, overrides,
    )
    encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)
    if args.output is not None:
        args.output.write_text(encoded + "\n", encoding="utf-8")
    print(encoded)


if __name__ == "__main__":
    main()
