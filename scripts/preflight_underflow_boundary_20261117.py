"""English description: Screen amplitude-boundary backgrounds and exact Fast frequency grids. 中文说明：预筛振幅边界背景并记录 Fast 实际选用的频率网格。"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_preflight_20261117.json"
TILTS = (0.0, 1.6, 17.30, 17.38, 17.45)
TEMPERATURES = (10.0, 2000.0, 10000.0)
KAPPAS = (1e-3, 1e-2, 1.0)
A_S = 2.100549e-9
BASE_PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": A_S,
    "DN_re": 10.0,
    "cr": 0.0,
}
SOURCE_FILES = (
    "stiffgwpy_fast/fast_sgwb.py",
    "stiffgwpy_fast/freq_adaptive.py",
    "stiffgwpy_fast/exact_background.py",
    "stiffgwpy_fast/reference.py",
    "stiffgwpy_fast/stiff_SGWB.py",
    "stiffgwpy_fast/LCDM_stiff_Neff.py",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes(order="C")).hexdigest()


def amplitude_boundary_r(scalar_amplitude: float) -> tuple[float, float]:
    """Return adjacent r values whose rounded products straddle float64 tiny."""
    tiny = np.finfo(np.float64).tiny
    r_normal = float(tiny / scalar_amplitude)
    while r_normal * scalar_amplitude < tiny:
        r_normal = float(np.nextafter(r_normal, math.inf))
    while float(np.nextafter(r_normal, 0.0)) * scalar_amplitude >= tiny:
        r_normal = float(np.nextafter(r_normal, 0.0))
    return r_normal, float(np.nextafter(r_normal, 0.0))


def screen_case(tilt: float, temperature: float, kappa: float, label: str, r: float) -> dict[str, object]:
    params = {
        **BASE_PARAMS,
        "r": r,
        "n_t": tilt,
        "T_re": temperature,
        "kappa10": kappa,
    }
    model = LCDM_SG(**params)
    derived = model.derived_param
    log_amplitude = float(FA._log_tensor_amplitude(model.cosmo_param, derived))
    exact_zero = r == 0.0
    stored_product = float(derived["A_t"])
    product_status = (
        "exact_zero_limit" if exact_zero else
        "rounded_zero_nonzero_inputs" if stored_product == 0.0 else
        "normal" if stored_product >= np.finfo(np.float64).tiny else
        "subnormal"
    )
    output: dict[str, object] = {
        "case_id": f"nt{tilt:g}_T{temperature:g}_k{kappa:g}_{label}",
        "parameters": params,
        "amplitude_label": label,
        "A_t_stored": stored_product,
        "log_A_t": log_amplitude if math.isfinite(log_amplitude) else None,
        "product_status": product_status,
        "exact_zero_tensor_limit": exact_zero,
        "N_re": float(derived["N_re"]),
        "N_inf": None if derived["N_inf"] is None else float(derived["N_inf"]),
        "actual_n_t": float(derived["nt"]),
        "background_grid_sha256": None,
        "frequency_grid_requested": "goal",
        "frequency_grid_used": None,
        "frequency_count": None,
        "frequency_sha256_f64le": None,
        "frequency_min_log10_hz": None,
        "frequency_max_log10_hz": None,
        "log_power_requires_log_domain": None,
        "log_P_t_min": None,
        "log_P_t_max": None,
        "fast_background_failure": None,
        "preflight_error": None,
    }
    try:
        cfg = FS.resolve_config("fast")
        FS.gen_fast(model, cfg.h, kink_split=cfg.kink_split)
        if getattr(model, "_fast_background_failure", None) is not None:
            output["fast_background_failure"] = model._fast_background_failure
            output["preflight_error"] = "fast_background_failure"
            return output
        if cfg.kink_split:
            FS._correct_kink_background(model)
        output["background_grid_sha256"] = {
            "Nv": digest(np.asarray(model.Nv, dtype=np.float64)),
            "sigma": digest(np.asarray(model.sigma, dtype=np.float64)),
            "f_hor": digest(np.asarray(model.f_hor, dtype=np.float64)),
        }
        if not math.isfinite(log_amplitude) and not exact_zero:
            output["preflight_error"] = "nonzero_input_without_finite_log_amplitude"
            return output
        log_amplitude_mode = (
            stored_product == 0.0 and not exact_zero and math.isfinite(log_amplitude)
            and float(model.cosmo_param["cr"]) <= 0.0
        )
        goal = FA.goal_oriented_freqs(
            model,
            1.0,
            seed_n=64,
            max_points=120,
            eval_freqs=None,
            allow_log_amplitude=log_amplitude_mode,
        )
        requires_log_power = (
            math.isfinite(log_amplitude)
            and FA._requires_log_tensor_power(
                log_amplitude, float(derived["nt"]), goal, gp.f_piv
            )
        )
        used = goal
        grid_name = "goal"
        if math.isfinite(log_amplitude) and requires_log_power:
            used = FA.construct_fast_frequency_grid(
                model, 1.0, log_amplitude=log_amplitude
            ).copy()
            grid_name = "construct"
        used = np.asarray(used, dtype=np.float64)
        if used.size < 2 or not np.isfinite(used).all() or not np.all(np.diff(used) < 0.0):
            output["preflight_error"] = "invalid_frequency_grid"
            return output
        log_pt = None
        if math.isfinite(log_amplitude):
            log_pt = log_amplitude + float(derived["nt"]) * (
                math.log(10.0) * used - math.log(gp.f_piv)
            )
        output.update({
            "frequency_grid_used": grid_name,
            "frequency_count": int(used.size),
            "frequency_sha256_f64le": digest(used),
            "frequency_min_log10_hz": float(np.min(used)),
            "frequency_max_log10_hz": float(np.max(used)),
            "log_power_requires_log_domain": bool(requires_log_power),
            "log_P_t_min": None if log_pt is None else float(np.min(log_pt)),
            "log_P_t_max": None if log_pt is None else float(np.max(log_pt)),
        })
    except Exception as exc:
        output["preflight_error"] = f"{type(exc).__name__}: {exc}"
    return output


def main() -> int:
    r_normal, r_subnormal = amplitude_boundary_r(A_S)
    amplitudes = (
        ("normal_edge", r_normal),
        ("subnormal_edge", r_subnormal),
        ("subnormal_interior", 1e-300),
        ("rounded_zero_product", 1e-316),
        ("min_subnormal_r", float(np.nextafter(0.0, 1.0))),
        ("exact_zero", 0.0),
    )
    rows = [
        screen_case(tilt, temperature, kappa, label, r)
        for tilt in TILTS
        for temperature in TEMPERATURES
        for kappa in KAPPAS
        for label, r in amplitudes
    ]
    failures = [row for row in rows if row["preflight_error"] is not None]
    input_amplitude_checks = [
        {
            "label": label,
            "r": value,
            "A_t_rounded": value * A_S,
            "category": (
                "normal" if value * A_S >= np.finfo(np.float64).tiny else
                "subnormal" if value * A_S > 0.0 else
                "rounded_zero"
            ),
        }
        for label, value in amplitudes
    ]
    output = {
        "schema_version": 1,
        "status": "complete" if not failures else "preflight_failures_recorded",
        "scope": "cr=0 construction-only preflight; no transfer propagation, Reference solve, accuracy claim, or continuous-domain certification",
        "grid_contract": "Mirror current fast preset goal grid: gen_fast with resolved fast kink split, goal_oriented_freqs(64,120), and production log-amplitude-candidate fallback to construct whenever _requires_log_tensor_power says float64 power is unrepresentable.",
        "design": {
            "tilts": list(TILTS),
            "T_re": list(TEMPERATURES),
            "kappa10": list(KAPPAS),
            "amplitude_strata": input_amplitude_checks,
            "row_count": len(rows),
            "independent_holdout": "No holdout is claimed in this construction-only screen; any direct-solve design must reserve points before inspecting solver errors.",
        },
        "source_sha256": {rel: sha256(ROOT / rel) for rel in SOURCE_FILES},
        "runner_sha256": sha256(Path(__file__)),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "failure_count": len(failures),
        "failure_case_ids": [str(row["case_id"]) for row in failures],
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": output["status"],
        "case_count": len(rows),
        "failure_count": len(failures),
        "amplitude_strata": input_amplitude_checks,
        "grid_counts": sorted({(row["frequency_grid_used"], row["frequency_count"]) for row in rows}),
        "categories": sorted({(row["product_status"], row["frequency_grid_used"]) for row in rows}),
        "output": str(OUTPUT),
    }, ensure_ascii=False, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
