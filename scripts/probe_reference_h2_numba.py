"""Prototype a Numba continuous-background H2/sigma kernel.

This leaves Reference production code untouched. It compares the compiled
kernel with the current Python implementation at all physical branch edges,
random in-range points, and representative scalar timings.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import statistics
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from numba import njit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "r": 0.19649188544593907,
    "cr": 0.0,
    "n_t": 1.0622501865308032,
    "T_re": 1761.314039059193,
    "kappa10": 0.011598483151630803,
    "DN_re": 10.0,
}
DN_EFF = 2.291411867022825
_FD_RHO_X, _FD_RHO_C = R._FD_RHO_CS.x, R._FD_RHO_CS.c
_FD_P_X, _FD_P_C = R._FD_P_CS.x, R._FD_P_CS.c
_THERMAL_RHO_X, _THERMAL_RHO_C = gp.spl_rho.x, gp.spl_rho.c
_THERMAL_RHOP_X, _THERMAL_RHOP_C = gp.spl_rhop.x, gp.spl_rhop.c


@njit(cache=False)
def _h2_sigma_compiled(
    N,
    Omh2,
    Osh2,
    Oerh2,
    Otrh2,
    Otreh2,
    OLh2,
    N_inf,
    N_re_abs,
    N_fin,
    N_max,
    nu_today,
    Omega_ph2,
    Omega_nh2,
    Omega_mnuh2,
):
    if N < N_re_abs:
        return 0.0, 1.0
    eN = math.exp(N_inf - N)
    e3N = eN * eN * eN
    nu = nu_today / eN
    if nu > 100.0:
        H2 = (
            Omh2 + Omega_mnuh2
            + (Omega_ph2 + 2.0 / 3.0 * Omega_nh2 + Oerh2) * eN
            + Osh2 * e3N + OLh2 / e3N
        )
        sigma = (
            Omh2 + Omega_mnuh2
            + 4.0 / 3.0 * (Omega_ph2 + 2.0 / 3.0 * Omega_nh2 + Oerh2) * eN
            + 2.0 * Osh2 * e3N
        ) / H2
    elif nu >= 0.1:
        lognu = math.log10(nu)
        rho_nu = R._cubic_spline_eval(_FD_RHO_X, _FD_RHO_C, lognu)
        p_nu = R._cubic_spline_eval(_FD_P_X, _FD_P_C, lognu)
        H2 = (
            Omh2
            + (Omega_ph2 + (2.0 / 3.0 + rho_nu / 3.0) * Omega_nh2 + Oerh2) * eN
            + Osh2 * e3N + OLh2 / e3N
        )
        sigma = (
            Omh2
            + 4.0 / 3.0 * (Omega_ph2 + 2.0 / 3.0 * Omega_nh2 + Oerh2) * eN
            + (rho_nu + p_nu) * Omega_nh2 / 3.0 * eN
            + 2.0 * Osh2 * e3N
        ) / H2
    elif N > N_inf - N_fin:
        H2 = Omh2 + Otrh2 * eN + Osh2 * e3N + OLh2 / e3N
        sigma = (Omh2 + 4.0 / 3.0 * Otrh2 * eN + 2.0 * Osh2 * e3N) / H2
    elif N >= N_inf - N_max:
        Nl = N_inf - N
        rho_i = R._cubic_spline_eval(_THERMAL_RHO_X, _THERMAL_RHO_C, Nl)
        rhop_i = R._cubic_spline_eval(_THERMAL_RHOP_X, _THERMAL_RHOP_C, Nl)
        H2 = Omh2 + (Omega_ph2 * rho_i + Oerh2) * eN + Osh2 * e3N + OLh2 / e3N
        sigma = (
            Omh2 + (Omega_ph2 * rhop_i + 4.0 / 3.0 * Oerh2) * eN
            + 2.0 * Osh2 * e3N
        ) / H2
    else:
        H2 = Omh2 + Otreh2 * eN + Osh2 * e3N + OLh2 / e3N
        sigma = (Omh2 + 4.0 / 3.0 * Otreh2 * eN + 2.0 * Osh2 * e3N) / H2
    return H2, sigma


def _args(bg):
    return (
        bg.Omh2, bg.Osh2, bg.Oerh2, bg.Otrh2, bg.Otreh2, bg.OLh2,
        bg.N_inf, bg.N_re_abs, bg.N_fin, bg.N_max, bg.nu_today,
        bg.Omega_ph2, bg.Omega_nh2, bg.Omega_mnuh2,
    )


def _compiled(bg, N):
    return _h2_sigma_compiled(float(N), *_args(bg))


def _git(*args):
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
    bg = R._Background(model, DN_EFF)
    R._warmup_numba_spline()

    boundary_points = [
        0.0,
        np.nextafter(bg.N_re_abs, -np.inf), bg.N_re_abs,
        np.nextafter(bg.N_re_abs, np.inf),
        bg.N_inf - bg.N_fin,
        np.nextafter(bg.N_inf - bg.N_fin, -np.inf),
        np.nextafter(bg.N_inf - bg.N_fin, np.inf),
        bg.N_inf - bg.N_max,
        np.nextafter(bg.N_inf - bg.N_max, -np.inf),
        np.nextafter(bg.N_inf - bg.N_max, np.inf),
        bg.N_inf,
    ]
    for nu_edge in (100.0, 0.1):
        n_edge = bg.N_inf + math.log(nu_edge / bg.nu_today)
        boundary_points.extend((np.nextafter(n_edge, -np.inf), n_edge, np.nextafter(n_edge, np.inf)))
    rng = np.random.default_rng(20261105)
    boundary_points.extend(rng.uniform(0.0, bg.N_inf, 2000))

    compile_started = time.perf_counter()
    _compiled(bg, bg.N_inf)
    compile_s = time.perf_counter() - compile_started
    worst_h2_rel = 0.0
    worst_sigma_abs = 0.0
    max_h2_ulp = 0
    compared = 0
    for N in boundary_points:
        expected = R._H2_and_sigma(bg, float(N))
        actual = _compiled(bg, float(N))
        if expected[0] != 0.0:
            worst_h2_rel = max(worst_h2_rel, abs(actual[0] - expected[0]) / abs(expected[0]))
            max_h2_ulp = max(max_h2_ulp, abs(int(np.float64(actual[0]).view(np.int64)) - int(np.float64(expected[0]).view(np.int64))))
        worst_sigma_abs = max(worst_sigma_abs, abs(actual[1] - expected[1]))
        compared += 1

    workload = np.asarray(boundary_points, dtype=float)
    # Warm each implementation, then alternate to limit drift.
    for N in workload[:16]:
        R._H2_and_sigma(bg, float(N))
        _compiled(bg, float(N))
    python_times, compiled_times = [], []
    for pair in range(8):
        order = ("python", "compiled") if pair % 2 == 0 else ("compiled", "python")
        for method in order:
            t0 = time.perf_counter()
            if method == "python":
                for N in workload:
                    R._H2_and_sigma(bg, float(N))
                python_times.append(time.perf_counter() - t0)
            else:
                for N in workload:
                    _compiled(bg, float(N))
                compiled_times.append(time.perf_counter() - t0)

    nu_100_N = bg.N_inf + math.log(200.0 / bg.nu_today)
    representative_times = {}
    for label, N in (("nu_gt_100", nu_100_N + 1.0), ("nu_near_100", nu_100_N), ("thermal", bg.N_inf - 10.0), ("late", 0.0)):
        t0 = time.perf_counter()
        for _ in range(20000):
            R._H2_and_sigma(bg, float(N))
        py_s = time.perf_counter() - t0
        t0 = time.perf_counter()
        for _ in range(20000):
            _compiled(bg, float(N))
        nb_s = time.perf_counter() - t0
        representative_times[label] = {"python_20k_s": py_s, "compiled_20k_s": nb_s, "speedup": py_s / nb_s}

    result = {
        "schema_version": 1,
        "experiment": "reference_h2_sigma_numba_prototype",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "reference_module_sha256": _hash(ROOT / "stiffgwpy_fast" / "reference.py"),
        "probe_script_sha256": _hash(Path(__file__)),
        "parameters": PARAMS,
        "DN_eff": DN_EFF,
        "background_N_inf": bg.N_inf,
        "background_N_re_abs": bg.N_re_abs,
        "point_count": compared,
        "random_seed": 20261105,
        "precision": {
            "worst_H2_relative_difference": worst_h2_rel,
            "worst_sigma_absolute_difference": worst_sigma_abs,
            "max_H2_ulp_difference": max_h2_ulp,
        },
        "compile_s_first_call": compile_s,
        "python_loop_s": python_times,
        "compiled_loop_s": compiled_times,
        "loop_median_speedup": statistics.median(python_times) / statistics.median(compiled_times),
        "representative_20k_call_timings": representative_times,
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "workers": 1,
        },
        "limitations": [
            "Prototype-only scalar background evaluation; no full DOP853 solve has been validated with this kernel.",
            "The compiled implementation uses the same equations and exact cubic coefficients, but still requires ODE-level accuracy testing before integration.",
        ],
    }
    out = ROOT / "docs" / "reference_h2_numba_probe_20261105.json"
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
