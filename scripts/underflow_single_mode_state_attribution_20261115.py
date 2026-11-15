"""Compare Fast and Cartesian Reference states at the same start and tail node."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.exact_background import exact_phi_s2_split  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
HANDOFF_INPUT = ROOT / "docs" / "research" / "results" / "underflow_handoff_pair_diagnostic_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_single_mode_state_attribution_20261115.json"
MODE_INDEX = 217
TAIL_ZS = (5.0, 8.0)
PHASE_MAX = 0.25
VARIANTS = ((0.01, False), (0.01, True), (0.005, False), (0.000625, True))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fast_handoff_state(
    params: dict, freq: float, dn_eff: float, h: float, sigma_exact: bool, tail_z: float
) -> dict:
    model = LCDM_SG(**params)
    model.cosmo_param["DN_eff"] = dn_eff
    model.gen_expansion()
    FS.gen_fast(model, h)
    kink_index, kink_fraction = FS._correct_kink_background(model)
    sigma, f_hor, j0s, z0s, _ = FS.prep_frequency_only(
        model, model.Nv, np.array([freq])
    )
    exact_primitive = None
    h_arr = None
    if sigma_exact:
        exact_primitive = exact_phi_s2_split(
            model, model.Nv, dn_eff, sigma_nodes=sigma
        )
    else:
        from stiffgwpy_fast.exact_background import fast_phi_s2_split

        exact_primitive = fast_phi_s2_split(
            model, model.Nv, dn_eff, sigma_nodes=sigma
        )
    Phi, Phi_mid, S2, S2inv, kink_index, kink_fraction, phi_re = exact_primitive

    j0 = int(j0s[0])
    z0 = float(z0s[0])
    xh = 0.0
    yh = math.exp(z0) * float(S2inv[j0])
    k = j0
    z = z0
    last_xh, last_yh, last_z = xh, yh, z
    kend = j0
    while k < len(model.Nv) - 1 and z < tail_z:
        step = float(h_arr[k]) if h_arr is not None else h
        z_mid = z0 + float(Phi_mid[k]) - float(Phi[j0])
        if k == kink_index and 0.0 < kink_fraction < 1.0:
            z_break = z0 + float(phi_re) - float(Phi[j0])
            z_end = z0 + float(Phi[k + 1]) - float(Phi[j0])
            xh, yh = FS._phase_segment(
                xh, yh, z, z_break, step * kink_fraction, PHASE_MAX
            )
            xh, yh = FS._phase_segment(
                xh, yh, z_break, z_end, step * (1.0 - kink_fraction), PHASE_MAX
            )
        else:
            z_end = 2.0 * z_mid - z
            xh, yh = FS._phase_segment(xh, yh, z, z_end, step, PHASE_MAX)
        k += 1
        z = z0 + float(Phi[k]) - float(Phi[j0])
        if z < tail_z:
            last_xh, last_yh, last_z = xh, yh, z
            kend = k

    if z >= tail_z:
        kend = max(j0, k - 1)
    scale = math.sqrt(float(S2[kend]))
    return {
        "h": h,
        "sigma_exact": sigma_exact,
        "kink_split": True,
        "primitive_method": "exact_phi_s2_split" if sigma_exact else "fast_phi_s2_split",
        "tail_z": tail_z,
        "time_node_count": int(len(model.Nv)),
        "frequency_start_index": j0,
        "frequency_start_N": float(model.Nv[j0]),
        "frequency_start_z": z0,
        "grid_f_hor_start": float(f_hor[j0]),
        "last_node_before_tail_index": int(kend),
        "last_node_before_tail_N": float(model.Nv[kend]),
        "last_node_before_tail_z": float(last_z),
        "fast_state_scaled": [float(last_xh), float(last_yh)],
        "fast_state_physical": [float(scale * last_xh), float(scale * last_yh)],
        "fast_S2_at_handoff": float(S2[kend]),
        "grid_sigma_at_handoff": float(sigma[kend]),
        "loop_exit_z": float(z),
    }


def _reference_state_at(params: dict, freq: float, dn_eff: float, N_start: float, z_start: float, N_end: float) -> dict:
    unit_params = dict(params)
    unit_params.update({"r": 1.0, "A_s": 1.0, "n_t": 0.0})
    model = LCDM_SG(**unit_params)
    model.cosmo_param["DN_eff"] = dn_eff
    model.gen_expansion()
    background = REF._Background(model, dn_eff)
    started = time.perf_counter()
    result = solve_ivp(
        REF._tensor_orig,
        (N_start, N_end),
        (z_start, 0.0, math.exp(z_start)),
        method="DOP853",
        rtol=1e-9,
        atol=(1e-13, 1e-13, 1e-13),
        args=(background,),
        dense_output=False,
    )
    elapsed = time.perf_counter() - started
    if not result.success:
        raise RuntimeError(f"Reference state integration failed: {result.message}")
    z, x, y = map(float, result.y[:, -1])
    return {
        "z": z,
        "state_physical": [x, y],
        "nfev": int(result.nfev),
        "runtime_s": elapsed,
    }


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    handoff_data = json.loads(HANDOFF_INPUT.read_text(encoding="utf-8"))
    params = data["parameters"]
    dn_eff = float(data["background"]["common_DN_eff"])
    freqs = np.asarray(data["frequency_grids"]["235"]["log10_hz"], dtype=np.float64)
    freq = float(freqs[MODE_INDEX])
    variants = []
    for h, sigma_exact in VARIANTS:
        for tail_z in TAIL_ZS:
            fast = _fast_handoff_state(params, freq, dn_eff, h, sigma_exact, tail_z)
            reference_state = _reference_state_at(
                params,
                freq,
                dn_eff,
                fast["frequency_start_N"],
                fast["frequency_start_z"],
                fast["last_node_before_tail_N"],
            )
            fast_xy = np.asarray(fast["fast_state_physical"], dtype=np.float64)
            ref_xy = np.asarray(reference_state["state_physical"], dtype=np.float64)
            rel_state = float(np.linalg.norm(fast_xy - ref_xy) / np.linalg.norm(ref_xy))
            rel_norm2 = float(np.dot(fast_xy, fast_xy) / np.dot(ref_xy, ref_xy) - 1.0)
            variants.append({
                **fast,
                "reference_state": reference_state,
                "z_difference": fast["last_node_before_tail_z"] - reference_state["z"],
                "relative_state_vector_difference": rel_state,
                "relative_quadratic_norm_difference": rel_norm2,
            })

    handoff_index = int(handoff_data["worst_node"]["index"])
    handoff_frequency = float(handoff_data["worst_node"]["log10_f_Hz"])
    production_state = next(
        row for row in variants
        if row["h"] == 0.005 and not row["sigma_exact"] and row["tail_z"] == 8.0
    )
    production_model = LCDM_SG(**params)
    production_model.cosmo_param["DN_eff"] = dn_eff
    production_model.gen_expansion()
    FS.gen_fast(production_model, 0.005)
    FS._correct_kink_background(production_model)
    f_hor_today = float(production_model.f_hor[-1])
    n_inf = float(production_model.derived_param["N_inf"])
    x_fast, y_fast = map(float, production_state["fast_state_physical"])
    z_handoff = float(production_state["last_node_before_tail_z"])
    n_handoff = float(production_state["last_node_before_tail_N"])
    inverse_omega = math.exp(-z_handoff)
    gamma = float(FS._tail_match_gamma(production_state["grid_sigma_at_handoff"]))
    coefficient_squared = 0.5 * (
        x_fast * x_fast
        + y_fast * y_fast * (1.0 + gamma * gamma * inverse_omega * inverse_omega)
        + 2.0 * gamma * x_fast * y_fast * inverse_omega
    )
    transfer_squared = coefficient_squared * math.exp(
        2.0 * (n_handoff - z_handoff - n_inf)
    )
    today_frequency_factor_squared = 10.0 ** (
        2.0 * (handoff_frequency - f_hor_today)
    )
    log_tensor_power = _log_tensor_amplitude(params, {"A_t": 0.0}) + float(
        params["n_t"]
    ) * (math.log(10.0) * handoff_frequency - math.log(REF.gp.f_piv))
    predicted_opgw = (
        transfer_squared * today_frequency_factor_squared / 36.0
        * math.exp(log_tensor_power)
    )
    actual_opgw = float(
        handoff_data["fast_z8"]["Opgw_today"][handoff_index]
    )
    reference_opgw = float(
        handoff_data["worst_node"]["Reference_z8_unit_terminal_moments"][
            "scaled_Opgw"
        ]
    )
    z_today = (handoff_frequency - f_hor_today) * math.log(10.0)
    phase_weight = float(
        FS._phase_resolved_tail_weights(np.array([z_today]), 8.0, enabled=True)[0]
    )
    terminal_scale_attribution = {
        "frequency_index": handoff_index,
        "frequency_log10_hz": handoff_frequency,
        "z_today": z_today,
        "phase_tail_weight": phase_weight,
        "reference_z_tail": 8.0,
        "fast_handoff_N": n_handoff,
        "fast_handoff_z": z_handoff,
        "fast_handoff_eps": float(
            handoff_data["fast_z8"]["handoff_eps"][handoff_index]
        ),
        "fast_gamma": gamma,
        "fast_coefficient_squared": coefficient_squared,
        "fast_f_hor_today": f_hor_today,
        "fast_N_inf": n_inf,
        "log_tensor_power": log_tensor_power,
        "predicted_fast_Opgw_today": predicted_opgw,
        "actual_fast_Opgw_today": actual_opgw,
        "reference_z8_scaled_Opgw_today": reference_opgw,
        "predicted_to_actual_fast_ratio": predicted_opgw / actual_opgw,
        "actual_fast_to_reference_ratio": actual_opgw / reference_opgw,
        "reference_state_relative_vector_difference": production_state[
            "relative_state_vector_difference"
        ],
        "reference_state_relative_quadratic_norm_difference": production_state[
            "relative_quadratic_norm_difference"
        ],
    }

    output = {
        "schema_version": 1,
        "status": "complete_fixed_mode_pre_handoff_state_attribution",
        "parameters": params,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": data["frequency_grids"]["235"]["sha256_f64le"],
        "frequency_index": MODE_INDEX,
        "frequency_log10_hz": freq,
        "frozen_DN_eff": dn_eff,
        "tail_z_targets": list(TAIL_ZS),
        "phase_max": PHASE_MAX,
        "reference_state_rtol": 1e-9,
        "variants": variants,
        "terminal_scale_attribution": terminal_scale_attribution,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "exact_background.py": _sha(ROOT / "stiffgwpy_fast" / "exact_background.py"),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "input_json": _sha(INPUT),
            "handoff_input_json": _sha(HANDOFF_INPUT),
            "diagnostic_script": _sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "variants": [
            {
                "h": row["h"],
                "sigma_exact": row["sigma_exact"],
                "tail_z": row["tail_z"],
                "N_start": row["frequency_start_N"],
                "z_start": row["frequency_start_z"],
                "N_handoff": row["last_node_before_tail_N"],
                "z_handoff": row["last_node_before_tail_z"],
                "z_difference": row["z_difference"],
                "relative_state_vector_difference": row["relative_state_vector_difference"],
                "relative_quadratic_norm_difference": row["relative_quadratic_norm_difference"],
            }
            for row in variants
        ],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
