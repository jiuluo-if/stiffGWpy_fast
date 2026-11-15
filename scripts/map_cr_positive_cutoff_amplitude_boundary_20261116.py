"""Map the cr>0 physical cutoff boundary against float64 tensor amplitude."""

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

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.LCDM_stiff_Neff import _has_exact_zero_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BASE = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "n_t": 0.0,
    "cr": 1.0,
    "DN_re": 20.0,
    "DN_eff": 0.0,
}
OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_cutoff_amplitude_boundary_20261116.json"
R_MIN = float(np.nextafter(0.0, 1.0))
AT_NORMAL_R = float(np.finfo(np.float64).tiny / BASE["A_s"])


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def delta_n_from_log_product(model: LCDM_SG, derived: dict[str, object], r: float) -> float:
    return (
        (gp.N_10 - float(derived["N_re"])) * 4.0 / 3.0
        + math.log(gp.M_Pl) * 4.0 / 3.0
        + (
            math.log(45.0 / 2.0)
            + math.log(r)
            + math.log(float(model.cosmo_param["A_s"]))
            - math.log(model.rhorad_re + model.rhostiff_re)
        ) / 3.0
    )


def evaluate(t_re: float, kappa: float, r: float) -> dict[str, object]:
    params = {**BASE, "T_re": t_re, "kappa10": kappa, "r": r}
    model = LCDM_SG(**params)
    derived = model.derived_param
    delta_n = delta_n_from_log_product(model, derived, r)
    at = float(derived["A_t"])
    return {
        "r": r,
        "log10_r": math.log10(r),
        "stored_A_t": at,
        "A_t_class": "rounded_zero" if at == 0.0 else "normal" if at >= np.finfo(float).tiny else "subnormal",
        "log_r_As": math.log(r) + math.log(float(params["A_s"])),
        "Delta_N_from_log_product": float(delta_n),
        "N_inf": None if derived["N_inf"] is None else float(derived["N_inf"]),
        "N_re": float(derived["N_re"]),
        "V_inf": derived.get("V_inf"),
        "exact_zero_tensor_limit": bool(_has_exact_zero_tensor_amplitude(model.cosmo_param)),
    }


def main() -> int:
    rows = []
    for t_re in (10.0, 2000.0, 10000.0):
        for kappa in (1e-3, 1e-2, 1.0):
            reference = LCDM_SG(**{**BASE, "T_re": t_re, "kappa10": kappa, "r": R_MIN})
            d = reference.derived_param
            log_r_cutoff = (
                -4.0 * (gp.N_10 - float(d["N_re"]))
                - 4.0 * math.log(gp.M_Pl)
                - math.log(45.0 / 2.0)
                - math.log(BASE["A_s"])
                + math.log(reference.rhorad_re + reference.rhostiff_re)
            )
            r_cutoff = math.exp(log_r_cutoff)
            below = evaluate(t_re, kappa, r_cutoff * (1.0 - 1e-6))
            above = evaluate(t_re, kappa, r_cutoff * (1.0 + 1e-6))
            subnormal = evaluate(t_re, kappa, AT_NORMAL_R * (1.0 - 1e-6))
            normal_boundary = evaluate(t_re, kappa, AT_NORMAL_R)
            product_zero = evaluate(t_re, kappa, 1e-315)
            exact_zero_model = LCDM_SG(**{**BASE, "T_re": t_re, "kappa10": kappa, "r": 0.0})
            exact_zero = bool(_has_exact_zero_tensor_amplitude(exact_zero_model.cosmo_param))
            rows.append({
                "T_re": t_re,
                "kappa10": kappa,
                "cutoff_root_log10_r": log_r_cutoff / math.log(10.0),
                "cutoff_root_r": r_cutoff,
                "delta_N_below_cutoff_root": below,
                "delta_N_above_cutoff_root": above,
                "immediate_subnormal_A_t_below_normal_boundary": subnormal,
                "A_t_normal_subnormal_boundary": normal_boundary,
                "rounded_zero_product_control_r_1e-315": product_zero,
                "cutoff_root_decades_above_A_t_normal_boundary": math.log10(r_cutoff / AT_NORMAL_R),
                "exact_r_zero_is_separate": exact_zero,
            })
    output = {
        "schema_version": 1,
        "status": "complete_cr_positive_cutoff_vs_A_t_representability_boundary",
        "scope": "For cr>0 at current low/mid/high T_re and kappa10, analytically solve Delta_N=0 for log(r) from the exact consistency equations, then verify direct derived_param signs on both sides and compare with the stored A_t normal/subnormal and rounded-zero limits. No invalid-cutoff point is propagated.",
        "parameters_except_T_re_kappa10_r": BASE,
        "A_t_normal_r_boundary": AT_NORMAL_R,
        "r_min_positive_subnormal": R_MIN,
        "cases": rows,
        "environment": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_sha256": {
            name: sha(ROOT / path)
            for name, path in {
                "LCDM_stiff_Neff.py": "stiffgwpy_fast/LCDM_stiff_Neff.py",
                "stiff_SGWB.py": "stiffgwpy_fast/stiff_SGWB.py",
                "global_param.py": "stiffgwpy_fast/global_param.py",
                "scan_consistency_underflow_boundary.py": "scripts/scan_consistency_underflow_boundary_20261115.py",
                "script": "scripts/map_cr_positive_cutoff_amplitude_boundary_20261116.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "cases": [{
            "T_re": row["T_re"], "kappa10": row["kappa10"],
            "log10_r_cutoff": row["cutoff_root_log10_r"],
            "cutoff_to_A_t_normal_decades": row["cutoff_root_decades_above_A_t_normal_boundary"],
            "below_delta_N": row["delta_N_below_cutoff_root"]["Delta_N_from_log_product"],
            "above_delta_N": row["delta_N_above_cutoff_root"]["Delta_N_from_log_product"],
            "subnormal_delta_N": row["immediate_subnormal_A_t_below_normal_boundary"]["Delta_N_from_log_product"],
            "subnormal_N_inf": row["immediate_subnormal_A_t_below_normal_boundary"]["N_inf"],
        } for row in rows],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
