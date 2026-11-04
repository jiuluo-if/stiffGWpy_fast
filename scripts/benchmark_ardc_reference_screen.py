"""Screen the authors' ARDC package on the archived res4 Reference grid.

The external riccati package is supplied through --riccati-path and is never a
project dependency. The scalar tensor equation is propagated piecewise at the
reheating boundary, then compared at identical DN_eff and frequencies to a
saved Cartesian DOP853 Reference spectrum.
"""

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = dict(r=0.01, cr=0.0, n_t=2.0, T_re=1e4, kappa10=0.01,
              DN_re=10.0)


def _vectorized(fn):
    def evaluate(values):
        array = np.asarray(values)
        if array.ndim == 0:
            return fn(float(array))
        return np.asarray([fn(float(value)) for value in array.ravel()]).reshape(
            array.shape)
    return evaluate


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _solve_mode(riccati, model, bg, freq, dn_eff, z_tail, eps, epsh):
    n_start = REF._find_start_N(bg, freq, -12.0)

    def z_at(n):
        return (freq - REF.f_hor_abs(bg, float(n))) * REF.ln10

    z_today = z_at(bg.N_inf)
    used_tail = z_today >= z_tail
    n_end = (brentq(lambda n: z_at(n) - z_tail, n_start, bg.N_inf,
                    xtol=5e-15, rtol=1e-14)
             if used_tail else bg.N_inf)

    def omega_scalar(n):
        try:
            value = math.exp(z_at(n))
        except OverflowError as exc:
            raise FloatingPointError("ARDC local frequency overflowed") from exc
        if value == 0.0:
            raise FloatingPointError("ARDC local frequency underflowed to zero")
        return value

    def gamma_scalar(n):
        _, sigma = REF._H2_and_sigma(bg, float(n))
        return 1.5 - 0.75 * sigma

    omega = _vectorized(omega_scalar)
    gamma = _vectorized(gamma_scalar)
    boundaries = [n_start]
    if n_start < bg.N_re_abs < n_end:
        boundaries.append(bg.N_re_abs)
    boundaries.append(n_end)

    tensor = 1.0 + 0.0j
    derivative = 0.0 + 0.0j
    mode_steps = 0
    riccati_steps = 0
    chebyshev_steps = 0
    linear_solves = 0
    segment_rows = []
    start = time.perf_counter()
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        info = riccati.solversetup(omega, gamma, h0=min(0.1, right - left))
        xs, ys, dys, successes, _, step_types, _, _ = riccati.solve(
            info, left, right, tensor, derivative, eps=eps, epsh=epsh,
            hard_stop=True)
        if not xs or not ys or not dys:
            raise RuntimeError("ARDC returned no integration steps")
        if not all(int(value) == 1 for value in successes):
            raise RuntimeError("ARDC reported an unsuccessful internal step")
        endpoint_error = abs(float(xs[-1]) - right)
        if endpoint_error > 1e-8:
            raise RuntimeError(
                "ARDC endpoint miss %.17g exceeds its 1e-8 stopping tolerance"
                % endpoint_error)
        tensor, derivative = ys[-1], dys[-1]
        stats = info.output(step_types)
        types = [int(value) for value in step_types]
        mode_steps += max(0, len(xs) - 1)
        riccati_steps += sum(value == 1 for value in types)
        chebyshev_steps += sum(value == 0 for value in types)
        linear_solves += int(stats["linear solves"])
        segment_rows.append({"left": left, "right": right,
                             "endpoint_error": endpoint_error,
                             "stats": _jsonable(stats)})
    runtime = time.perf_counter() - start

    tensor_real = float(np.real(tensor))
    derivative_real = float(np.real(derivative))
    imaginary_ratio = max(abs(float(np.imag(tensor))) /
                          max(abs(tensor_real), 1e-300),
                          abs(float(np.imag(derivative))) /
                          max(abs(derivative_real), 1e-300))
    z_end = z_at(n_end)
    y_state = omega_scalar(n_end) * tensor_real
    x_state = derivative_real
    tensor_power = model.derived_param["A_t"] * (
        10.0 ** freq / gp.f_piv) ** model.derived_param["nt"]

    if used_tail:
        coefficient = math.sqrt(0.5 * (x_state * x_state + y_state * y_state))
        transfer = coefficient * math.exp(-z_end + n_end - bg.N_inf)
        x_today = transfer * 10.0 ** (
            freq - REF.f_hor_abs(bg, bg.N_inf))
        oj = -transfer * transfer / 3.0 * tensor_power
        opgw = x_today * x_today / 36.0 * tensor_power
        ogw = 3.0 * opgw + oj
    else:
        tensor_today = y_state / math.exp(z_end)
        oj = x_state * tensor_today / 3.0 * tensor_power
        ogw = (x_state * x_state + y_state * y_state) / 24.0 * tensor_power + oj
        opgw = (-5.0 * x_state * x_state + 7.0 * y_state * y_state) \
            / 72.0 * tensor_power

    omega_today = ogw - oj
    if not all(math.isfinite(float(value))
               for value in (tensor_real, derivative_real, omega_today,
                             oj, opgw, ogw, imaginary_ratio)):
        raise FloatingPointError("ARDC state or spectrum is non-finite")
    log_omega = math.log10(max(omega_today, 1e-40))
    return {
        "freq_log10_hz": float(freq),
        "used_tail": bool(used_tail),
        "event_N": float(n_end) if used_tail else None,
        "runtime_s": runtime,
        "steps": int(mode_steps),
        "riccati_steps": int(riccati_steps),
        "chebyshev_steps": int(chebyshev_steps),
        "linear_solves": int(linear_solves),
        "imaginary_state_ratio": imaginary_ratio,
        "log10OmegaGW": float(log_omega),
        "Ogw": float(ogw),
        "Oj": float(oj),
        "Opgw": float(opgw),
        "segments": segment_rows,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--riccati-path", type=Path, required=True,
                        help="isolated directory containing riccati, not a repo dependency")
    parser.add_argument("--reference", type=Path,
                        default=ROOT / "docs" / "res4_sigma_exact_reference_map_20261104.npz")
    parser.add_argument("--reference-map", type=Path,
                        default=ROOT / "docs" / "res4_sigma_exact_reference_map_20261104.json")
    parser.add_argument("--dn-eff", type=float, default=0.133732083004284)
    parser.add_argument("--z-tail", type=float, default=10.0)
    parser.add_argument("--eps", type=float, default=1e-9)
    parser.add_argument("--epsh", type=float, default=1e-9)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs" / "reference_ardc_res4_screen_20261104.json")
    parser.add_argument("--progress-output", type=Path, default=None,
                        help="flushed JSONL checkpoint, written after every mode")
    args = parser.parse_args(argv)
    sys.path.insert(0, str(args.riccati_path))
    import riccati  # pylint: disable=import-outside-toplevel

    with np.load(args.reference) as archive:
        freqs = np.asarray(archive["freqs"], dtype=float)
        reference_log = np.asarray(archive["reference_log"], dtype=float)
    if freqs.shape != reference_log.shape:
        parser.error("archived frequencies and spectrum must have equal shape")

    model = LCDM_SG(**PARAMS, DN_eff=args.dn_eff)
    bg = REF._Background(model, args.dn_eff)
    progress_path = (args.progress_output if args.progress_output is not None
                     else args.output.with_suffix(".progress.jsonl"))
    progress_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    with progress_path.open("w", encoding="utf-8") as progress:
        progress.write(json.dumps({
            "status": "started",
            "source_sha": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "package_version": riccati.__version__,
            "dn_eff": args.dn_eff,
            "z_tail": args.z_tail,
            "eps": args.eps,
            "epsh": args.epsh,
            "frequency_count": int(len(freqs)),
            "resources": telemetry(workers=1, threads=1),
        }) + "\n")
        progress.flush()
        for index, freq in enumerate(freqs):
            try:
                row = _solve_mode(riccati, model, bg, float(freq), args.dn_eff,
                                  args.z_tail, args.eps, args.epsh)
            except KeyboardInterrupt:
                progress.write(json.dumps({
                    "status": "interrupted",
                    "next_index": int(index),
                    "completed_modes": len(rows),
                }) + "\n")
                progress.flush()
                raise
            except Exception as exc:
                failure = {
                    "status": "failed",
                    "source_sha": subprocess.check_output(
                        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                    "failed_index": int(index),
                    "freq_log10_hz": float(freq),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "completed_modes": len(rows),
                    "resources": telemetry(workers=1, threads=1),
                }
                progress.write(json.dumps(failure) + "\n")
                progress.flush()
                args.output.write_text(json.dumps(failure, indent=2) + "\n",
                                       encoding="utf-8")
                raise
            row["index"] = int(index)
            row["abs_dex_error"] = abs(row["log10OmegaGW"] - reference_log[index])
            rows.append(row)
            progress.write(json.dumps({"status": "mode_completed", **row}) + "\n")
            progress.flush()

    g2, quadrature_error, interpolation_error = REF.integrate_spectrum(
        freqs, np.asarray([row["Ogw"] for row in rows]),
        np.asarray([row["Oj"] for row in rows]))
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    candidate_dn = gp.Neff0 * g2 / omega_nu
    with args.reference_map.open(encoding="utf-8") as handle:
        reference_map = json.load(handle)
    errors = np.asarray([row["abs_dex_error"] for row in rows])
    runtimes = np.asarray([row["runtime_s"] for row in rows])
    payload = {
        "schema_version": 1,
        "experiment": "author_ARDC_res4_full_grid_vs_cartesian_reference",
        "source_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "package": "riccati",
        "package_version": riccati.__version__,
        "parameters": {**PARAMS, "DN_eff": args.dn_eff},
        "z_tail": args.z_tail,
        "eps": args.eps,
        "epsh": args.epsh,
        "reference_rtol": 1e-11,
        "frequency_count": int(len(freqs)),
        "resources": telemetry(workers=1, threads=1),
        "candidate_DN_gw": float(candidate_dn),
        "reference_DN_gw": float(reference_map["reference_dn_map"]),
        "DN_gw_relative_error_vs_reference": float(
            abs(candidate_dn - reference_map["reference_dn_map"])
            / abs(reference_map["reference_dn_map"])),
        "candidate_map_residual_relative": float(
            (candidate_dn - args.dn_eff) / args.dn_eff),
        "quadrature_error": quadrature_error,
        "interpolation_error": interpolation_error,
        "spectrum_error_median_dex": float(np.median(errors)),
        "spectrum_error_p95_dex": float(np.percentile(errors, 95)),
        "spectrum_error_max_dex": float(np.max(errors)),
        "mode_runtime_sum_s": float(np.sum(runtimes)),
        "mode_runtime_p50_s": float(np.percentile(runtimes, 50)),
        "mode_runtime_p95_s": float(np.percentile(runtimes, 95)),
        "used_tail_count": int(sum(row["used_tail"] for row in rows)),
        "max_imaginary_state_ratio": float(max(
            row["imaginary_state_ratio"] for row in rows)),
        "rows": rows,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    with progress_path.open("a", encoding="utf-8") as progress:
        progress.write(json.dumps({"status": "complete",
                                   "completed_modes": len(rows),
                                   "candidate_DN_gw": float(candidate_dn),
                                   "spectrum_error_max_dex": payload[
                                       "spectrum_error_max_dex"]}) + "\n")
        progress.flush()
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
