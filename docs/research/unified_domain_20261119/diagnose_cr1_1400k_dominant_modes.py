"""Compare one Fast trajectory with the independent tensor ODE at shared nodes."""

# The script is run by path and adds the repository root before importing the
# local checkout rather than whichever package version is installed.
# ruff: noqa: E402, I001

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from numba import njit
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast import reference as REF
from stiffgwpy_fast.exact_background import fast_phi_s2_split
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


RESEARCH = Path(__file__).resolve().parent
_phase_segment = FS._phase_segment
PARAMS = {
    "Omega_bh2": 7.0e5,
    "Omega_ch2": 7.0e5,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 1.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
FREQUENCIES = (-15.40875682964938, -15.040578320076837, -14.672399810504294, -14.251624372050038)


@njit(cache=True)
def propagate_mode(Nv, Phi, Phi_mid, S2inv, j0, z0, col_step, h, phase_max,
                   kink_index, kink_fraction, phi_re, tail_z, curved_tolerance):
    count = (len(Nv) - 1 - j0) // col_step + 2
    indices = np.empty(count, dtype=np.int64)
    zs = np.empty(count)
    xs = np.empty(count)
    ys = np.empty(count)
    size = 0
    xh = 0.0
    yh = math.exp(z0) * S2inv[j0]
    log_gain = 0.0
    P0 = Phi[j0] - Nv[j0]
    k = j0
    z = z0
    indices[size] = k
    zs[size] = z
    xs[size] = xh
    ys[size] = yh
    size += 1
    while k < len(Nv) - 1 and z < tail_z:
        step = Nv[k + 1] - Nv[k]
        if k == kink_index and 0.0 < kink_fraction < 1.0:
            z_break = z0 + phi_re - Phi[j0]
            z_end = z0 + Phi[k + 1] - Phi[j0]
            xh, yh = _phase_segment(
                xh, yh, z, z_break, step * kink_fraction, phase_max
            )
            xh, yh = _phase_segment(
                xh, yh, z_break, z_end, step * (1.0 - kink_fraction), phase_max
            )
        else:
            z_mid_step = z0 + Phi_mid[k] - Phi[j0]
            z_end = 2.0 * z_mid_step - z
            z_exact_end = z0 + Phi[k + 1] - Phi[j0]
            if abs(z_exact_end - z_end) > curved_tolerance:
                # Quadratic interpolation through exact left/mid/right z values;
                # evaluate two composite-midpoint transfer panels.
                zq1 = 0.375 * z + 0.75 * z_mid_step - 0.125 * z_exact_end
                zq3 = -0.125 * z + 0.75 * z_mid_step + 0.375 * z_exact_end
                xh, yh = _phase_segment(
                    xh, yh, zq1, zq1, 0.5 * step, phase_max
                )
                xh, yh = _phase_segment(
                    xh, yh, zq3, zq3, 0.5 * step, phase_max
                )
            else:
                xh, yh = _phase_segment(xh, yh, z, z_end, step, phase_max)
        size_state = max(abs(xh), abs(yh))
        if size_state > 2.0**450 or 0.0 < size_state < 2.0**-450:
            xh /= size_state
            yh /= size_state
            log_gain += math.log(size_state)
        k += 1
        z = z0 + Phi[k] - Phi[j0]
        if k % col_step == 0 or k == len(Nv) - 1:
            scale = math.exp((Phi[k] - Nv[k] - P0) + log_gain)
            indices[size] = k
            zs[size] = z
            xs[size] = scale * xh
            ys[size] = scale * yh
            size += 1
    return indices[:size], zs[:size], xs[:size], ys[:size]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    model = LCDM_SG(**PARAMS)
    if FS.SGWB_iter_fast(
        model,
        config=FS.resolve_config("fast", threads=2),
                tol=1e-6,
    ) is None:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    dn_eff = float(model.DN_gw[-1])
    Nv = np.asarray(model.Nv, dtype=float)
    freqs = np.asarray(FREQUENCIES, dtype=float)
    sigma, f_hor, j0s, z0s, _ = FS.prep_frequency_only(model, Nv, freqs)
    Phi, Phi_mid, _, S2inv, kink_index, kink_fraction, phi_re = fast_phi_s2_split(
        model, Nv, dn_eff, sigma_nodes=model.sigma, relative_scale=True
    )
    bg = REF._Background(model, dn_eff)
    rows = []
    for mode, freq in enumerate(freqs):
        j0 = int(j0s[mode])
        z0 = float(z0s[mode])
        variants = {}
        for name, threshold in (("base", math.inf), ("curved_1e-5", 1e-5), ("curved_1e-4", 1e-4)):
            indices_v, z_v, x_v, y_v = propagate_mode(
                Nv, Phi, Phi_mid, S2inv, j0, z0, 8, 0.005, 0.25,
                int(kink_index), float(kink_fraction), float(phi_re), 5.0,
                threshold,
            )
            before_tail = z_v <= 5.0
            variants[name] = (
                indices_v[before_tail], z_v[before_tail],
                x_v[before_tail], y_v[before_tail],
            )
        indices, z_fast, x_fast, y_fast = variants["base"]
        sample_N = Nv[indices]
        event_n = float(bg.N_re_abs)
        if sample_N[0] < event_n < sample_N[-1]:
            left_samples = sample_N[sample_N <= event_n]
            right_samples = sample_N[sample_N > event_n]
            left = solve_ivp(
                lambda n, state: REF._tensor_orig(n, state, bg),
                (float(sample_N[0]), event_n),
                (z0, 0.0, math.exp(z0)),
                method="DOP853",
                rtol=2e-11,
                atol=(1e-12, 1e-22, 1e-22),
                t_eval=left_samples,
            )
            right = solve_ivp(
                lambda n, state: REF._tensor_orig(n, state, bg),
                (event_n, float(sample_N[-1])),
                left.y[:, -1],
                method="DOP853",
                rtol=2e-11,
                atol=(1e-12, 1e-22, 1e-22),
                t_eval=right_samples,
            )
            if not left.success or not right.success:
                raise RuntimeError(left.message if not left.success else right.message)
            reference_values = np.concatenate((left.y, right.y), axis=1)
        else:
            reference = solve_ivp(
                lambda n, state: REF._tensor_orig(n, state, bg),
                (float(sample_N[0]), float(sample_N[-1])),
                (z0, 0.0, math.exp(z0)),
                method="DOP853",
                rtol=2e-11,
                atol=(1e-12, 1e-22, 1e-22),
                t_eval=sample_N,
            )
            if not reference.success:
                raise RuntimeError(reference.message)
            reference_values = reference.y
        z_ref, x_ref, y_ref = reference_values
        state_scale = np.maximum(np.maximum(np.abs(x_ref), np.abs(y_ref)), 1e-250)
        for name, (indices_v, z_v, x_v, y_v) in variants.items():
            state_error = np.maximum(np.abs(x_v - x_ref), np.abs(y_v - y_ref)) / state_scale
            rows.append({
                "variant": name,
                "frequency": float(freq),
                "j0": j0,
                "N_start": float(Nv[j0]),
                "z_start": z0,
                "n_shared_nodes": int(indices_v.size),
                "max_abs_z_difference": float(np.max(np.abs(z0 + Phi[indices_v] - Phi[j0] - z_ref))),
                "max_relative_state_difference": float(np.max(state_error)),
                "p95_relative_state_difference": float(np.percentile(state_error, 95)),
                "sample_N": Nv[indices_v].tolist(),
                "fast_z": z_v.tolist(),
                "reference_z": z_ref.tolist(),
                "fast_x": x_v.tolist(),
                "reference_x": x_ref.tolist(),
                "fast_y": y_v.tolist(),
                "reference_y": y_ref.tolist(),
            })
    record = {
        "parameters": PARAMS,
        "DN_eff": dn_eff,
        "comparison_stops_before_z_tail_5": True,
        "source_head": __import__("subprocess").check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256": {
            "fast_sgwb.py": sha256(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "exact_background.py": sha256(ROOT / "stiffgwpy_fast" / "exact_background.py"),
            "reference.py": sha256(ROOT / "stiffgwpy_fast" / "reference.py"),
        },
        "modes": rows,
    }
    (RESEARCH / "cr1_1400k_dominant_mode_trajectories.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "DN_eff": dn_eff,
        "modes": [{k: v for k, v in row.items() if not isinstance(v, list)} for row in rows],
    }, indent=2))


if __name__ == "__main__":
    main()
