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
from stiffgwpy_fast.exact_background import exact_phi_s2  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_single_mode_state_attribution_20261115.json"
MODE_INDEX = 215
TAIL_Z = 5.0
VARIANTS = ((0.01, False), (0.01, True), (0.000625, True))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fast_handoff_state(params: dict, freq: float, dn_eff: float, h: float, sigma_exact: bool) -> dict:
    model = LCDM_SG(**params)
    model.cosmo_param["DN_eff"] = dn_eff
    FS.gen_fast(model, h)
    prepared = FS.prep_fast(model, model.Nv, np.array([freq]), h, variable_grid=False)
    sigma, f_hor, Phi, Phi_mid, _, S2, S2inv, j0s, z0s, _ = prepared
    h_arr = None
    if sigma_exact:
        Phi, Phi_mid, S2, S2inv = exact_phi_s2(model, model.Nv, dn_eff, h)
        h_arr = np.diff(model.Nv)

    j0 = int(j0s[0])
    z0 = float(z0s[0])
    xh = 0.0
    yh = math.exp(z0) * float(S2inv[j0])
    k = j0
    z = z0
    last_xh, last_yh, last_z = xh, yh, z
    kend = j0
    while k < len(model.Nv) - 1 and z < TAIL_Z:
        step = float(h_arr[k]) if h_arr is not None else h
        z_mid = z0 + float(Phi_mid[k]) - float(Phi[j0])
        z_end = 2.0 * z_mid - z
        xh, yh = FS._phase_segment(xh, yh, z, z_end, step, 0.0)
        k += 1
        z = z0 + float(Phi[k]) - float(Phi[j0])
        if z < TAIL_Z:
            last_xh, last_yh, last_z = xh, yh, z
            kend = k

    if z >= TAIL_Z:
        kend = max(j0, k - 1)
    scale = math.sqrt(float(S2[kend]))
    return {
        "h": h,
        "sigma_exact": sigma_exact,
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
    params = data["parameters"]
    dn_eff = float(data["background"]["common_DN_eff"])
    freqs = np.asarray(data["frequency_grids"]["235"]["log10_hz"], dtype=np.float64)
    freq = float(freqs[MODE_INDEX])
    variants = []
    for h, sigma_exact in VARIANTS:
        fast = _fast_handoff_state(params, freq, dn_eff, h, sigma_exact)
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

    output = {
        "schema_version": 1,
        "status": "complete_fixed_mode_pre_handoff_state_attribution",
        "parameters": params,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": data["frequency_grids"]["235"]["sha256_f64le"],
        "frequency_index": MODE_INDEX,
        "frequency_log10_hz": freq,
        "frozen_DN_eff": dn_eff,
        "tail_z": TAIL_Z,
        "reference_state_rtol": 1e-9,
        "variants": variants,
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
