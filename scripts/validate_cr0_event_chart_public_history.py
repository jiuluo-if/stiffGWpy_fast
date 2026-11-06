"""Compare event-chart public history with full-grid Fast and two References.

This is an isolated validation harness. It does not modify the prototype,
production Fast source, or Reference implementation.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, telemetry  # noqa: E402

apply_environment()

import numpy as np  # noqa: E402
from numba import njit, prange  # noqa: E402

from scripts import prototype_cr0_event_chart_kernel as prototype  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.exact_background import (  # noqa: E402
    build_transition_grid,
    exact_phi_s2_split,
)
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "endpoint_admissible_wkb_tail_20261105.json"
OUTPUT = ROOT / "docs" / "cr0_event_chart_public_history_20261106.json"
OUTPUT_REF5 = ROOT / "docs" / "cr0_event_chart_public_history_reference_z5_20261106.json"
OUTPUT_REF8 = ROOT / "docs" / "cr0_event_chart_public_history_reference_z8_20261106.json"
PRIOR_PROFILE = ROOT / "docs" / "cr0_event_chart_deferred_profile_20261106.json"
Z_TAIL = prototype.Z_TAIL
H = prototype.H
COL_STEP = prototype.COL_STEP
PHASE_MAX = prototype.PHASE_MAX


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _state_value(value):
    scalar = float(value)
    return scalar if math.isfinite(scalar) else None


def public_history_indices(Nv: np.ndarray, col_step: int = COL_STEP) -> np.ndarray:
    """Return the production coarse history layout, including terminal node."""
    grid = np.asarray(Nv, dtype=np.float64)
    if grid.ndim != 1 or grid.size < 2 or col_step < 1:
        raise ValueError("Nv must be a 1-D grid with at least two nodes and positive col_step")
    return np.unique(np.append(np.arange(0, grid.size, col_step), grid.size - 1)).astype(np.int64)


def chart_history_indices(
    public_Nv: np.ndarray,
    idx_out: np.ndarray,
    chart_N_abs: np.ndarray,
    N_star: float,
    *,
    atol: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray]:
    """Map post-event public columns to chart nodes; pre-event slots map to -1."""
    nodes = np.asarray(public_Nv, dtype=np.float64)[np.asarray(idx_out, dtype=np.int64)]
    chart_nodes = np.asarray(chart_N_abs, dtype=np.float64)
    mapping = np.full(nodes.size, -1, dtype=np.int64)
    for slot, value in enumerate(nodes):
        if value < N_star:
            continue
        pos = int(np.searchsorted(chart_nodes, value))
        candidates = (pos - 1, pos)
        best = min((q for q in candidates if 0 <= q < chart_nodes.size),
                   key=lambda q: abs(float(chart_nodes[q]) - float(value)))
        if abs(float(chart_nodes[best]) - float(value)) > atol:
            raise ValueError(f"public post-event node {value:.17g} is absent from chart")
        mapping[slot] = best
    return nodes, mapping


@njit(inline="always")
def _assemble_tail_history(
    Ogw, Oj, Opgw, mode, sample_chart_indices, sample_tau, anchor_tau,
    last_z, xh, yh, s2k, sigma, fp_freq, P_t, ev_minus, fp_minus,
    stop_slot,
):
    e_z = math.exp(-last_z)
    gamma = FS._tail_match_gamma(sigma)
    amp2 = (
        xh * xh + yh * yh * (1.0 + gamma * gamma * e_z * e_z)
        + 2.0 * gamma * xh * yh * e_z
    )
    coeff = math.sqrt(0.5 * s2k * amp2)
    e_nz = math.exp(anchor_tau - last_z)
    for slot in range(stop_slot):
        chart_index = sample_chart_indices[slot]
        if chart_index >= 0 and sample_tau[slot] >= anchor_tau:
            FS.assemble_tail(
                Ogw, Oj, Opgw, mode, slot, chart_index, coeff, e_nz,
                fp_freq, P_t, ev_minus, fp_minus,
            )


@njit(parallel=True, cache=True)
def _candidate_history_kernel(
    Nv, Phi_grid, Phi_mid, S2, start_indices, start_z, z_offsets,
    start_tau, direct_tail, start_xh, start_yh, P_t, ev_minus,
    fp_minus, fp_freq, h_arr, Sv, phase_max, z_tail, deferred_mask,
    sample_chart_indices, sample_tau, Ogw, Oj, Opgw,
    anchor_z, anchor_x, anchor_y, terminal_z, terminal_x, terminal_y,
):
    terminal = len(Nv) - 1
    last_slot = len(sample_chart_indices) - 1
    for mode in prange(len(start_indices)):
        if direct_tail[mode]:
            z_anchor = start_z[mode]
            tau_anchor = start_tau[mode]
            e_z = math.exp(-z_anchor)
            s2k = math.exp(-tau_anchor)
            gamma = FS._tail_match_gamma(1.0)
            amp2 = (
                start_xh[mode] * start_xh[mode]
                + start_yh[mode] * start_yh[mode] * (1.0 + gamma * gamma * e_z * e_z)
                + 2.0 * gamma * start_xh[mode] * start_yh[mode] * e_z
            )
            coeff = math.sqrt(0.5 * s2k * amp2)
            e_nz = math.exp(tau_anchor - z_anchor)
            anchor_z[mode] = z_anchor
            anchor_x[mode] = start_xh[mode]
            anchor_y[mode] = start_yh[mode]
            terminal_z[mode] = z_anchor + Nv[terminal] - tau_anchor
            terminal_x[mode] = 0.0
            terminal_y[mode] = 0.0
            for slot in range(len(sample_chart_indices)):
                chart_index = sample_chart_indices[slot]
                if chart_index >= 0 and sample_tau[slot] >= tau_anchor:
                    FS.assemble_tail(
                        Ogw, Oj, Opgw, mode, slot, chart_index, coeff, e_nz,
                        fp_freq[mode], P_t[mode], ev_minus, fp_minus,
                    )
            continue

        start = start_indices[mode]
        k = start
        zz = start_z[mode]
        zoff = z_offsets[mode]
        # Post-event starts have x=0; pre-event prefixes arrive at the event
        # with the analytically transferred state.
        xh = start_xh[mode]
        yh = start_yh[mode]
        lxh, lyh, last_z = xh, yh, zz
        for slot in range(len(sample_chart_indices)):
            if sample_chart_indices[slot] == k:
                FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[k], xh, yh, zz, P_t[mode])

        deferred = deferred_mask[mode]
        phase_max_mode = phase_max
        if deferred and 0.25 <= phase_max <= FS._WKB_TAIL_DEFERRED_PHASE_MAX:
            phase_max_mode = FS._WKB_TAIL_DEFERRED_PHASE_MAX
        z_tail_mode = math.inf if deferred else z_tail
        tail_saved = False
        anchor_z[mode] = math.nan
        anchor_x[mode] = math.nan
        anchor_y[mode] = math.nan

        if deferred and zz >= z_tail:
            anchor_z[mode], anchor_x[mode], anchor_y[mode] = last_z, lxh, lyh
            _assemble_tail_history(
                Ogw, Oj, Opgw, mode, sample_chart_indices, sample_tau,
                Nv[k], last_z, lxh, lyh, S2[k], Sv[k], fp_freq[mode],
                P_t[mode], ev_minus, fp_minus, last_slot,
            )
            tail_saved = True

        while k < terminal and zz < z_tail_mode:
            h_step = h_arr[k]
            z_node = zz
            z_mid = zoff + Phi_mid[k]
            z_end = 2.0 * z_mid - z_node
            xh, yh = FS._phase_segment(xh, yh, z_node, z_end, h_step, phase_max_mode)
            k += 1
            zz = zoff + Phi_grid[k]

            if deferred and not tail_saved and zz >= z_tail:
                anchor_z[mode], anchor_x[mode], anchor_y[mode] = last_z, lxh, lyh
                _assemble_tail_history(
                    Ogw, Oj, Opgw, mode, sample_chart_indices, sample_tau,
                    Nv[k - 1], last_z, lxh, lyh, S2[k - 1], Sv[k - 1],
                    fp_freq[mode], P_t[mode], ev_minus, fp_minus, last_slot,
                )
                tail_saved = True

            for slot in range(len(sample_chart_indices)):
                if sample_chart_indices[slot] == k and not (
                    deferred and tail_saved and k < terminal
                ):
                    FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[k], xh, yh, zz, P_t[mode])
            if zz < z_tail:
                lxh, lyh, last_z = xh, yh, zz

        terminal_z[mode] = zz
        terminal_x[mode] = xh
        terminal_y[mode] = yh
        if deferred:
            if not tail_saved and zz >= z_tail:
                anchor_z[mode], anchor_x[mode], anchor_y[mode] = last_z, lxh, lyh
                _assemble_tail_history(
                    Ogw, Oj, Opgw, mode, sample_chart_indices, sample_tau,
                    Nv[k - 1], last_z, lxh, lyh, S2[k - 1], Sv[k - 1],
                    fp_freq[mode], P_t[mode], ev_minus, fp_minus, last_slot,
                )
            for slot in range(len(sample_chart_indices)):
                if sample_chart_indices[slot] == terminal:
                    FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[terminal], xh, yh, zz, P_t[mode])
        elif zz < z_tail:
            for slot in range(len(sample_chart_indices)):
                if sample_chart_indices[slot] == terminal:
                    FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[terminal], xh, yh, zz, P_t[mode])
        else:
            kend = k - 1
            if kend < start:
                kend = start
            anchor_z[mode], anchor_x[mode], anchor_y[mode] = last_z, lxh, lyh
            _assemble_tail_history(
                Ogw, Oj, Opgw, mode, sample_chart_indices, sample_tau,
                Nv[kend], last_z, lxh, lyh, S2[kend], Sv[kend],
                fp_freq[mode], P_t[mode], ev_minus, fp_minus,
                len(sample_chart_indices),
            )


def _history_aggregate(freqs, j0, idx_out, outputs):
    Xf = np.flip(np.asarray(freqs, dtype=np.float64))
    Wmat = np.zeros((freqs.size, freqs.size), dtype=np.float64)
    FS.build_Wmat(freqs.size, Xf, np.diff(Xf), Wmat)
    j_hi = np.searchsorted(j0, idx_out, side="right") - 1
    g2c = np.zeros(idx_out.size, dtype=np.float64)
    w2c = np.zeros(idx_out.size, dtype=np.float64)
    FS.int_SGWB_W(
        freqs.size, idx_out.size, j_hi.astype(np.int64), Wmat,
        outputs[0], outputs[1], outputs[2], g2c, w2c, math.log(10.0),
    )
    g2c[-1] = FS.integrate_frequency_quadrature(
        freqs, outputs[0][:, -1] - outputs[1][:, -1], "pchip"
    ) * math.log(10.0)
    return j_hi.astype(np.int64), g2c, w2c


def _public_history_arrays(model, freqs, defer_endpoint):
    Nv, sigma, f_hor, _ = build_transition_grid(model, H)
    Phi, Phi_mid, S2, S2inv, kink_index, kink_fraction, phi_re = exact_phi_s2_split(
        model, Nv, model.cosmo_param["DN_eff"], sigma_nodes=sigma
    )
    _, _, j0, z0, fp_minus = FS.prep_frequency_only(model, Nv, freqs)
    endpoint_z = z0 + Phi[-1] - Phi[j0]
    endpoint_defect = np.fromiter(
        (FS._endpoint_wkb_defect(float(sigma[-1]), float(value)) for value in endpoint_z),
        dtype=np.float64,
        count=endpoint_z.size,
    )
    idx_out = public_history_indices(Nv)
    outputs = [np.zeros((freqs.size, idx_out.size), dtype=np.float64) for _ in range(3)]
    ev_minus = np.exp(-Nv)
    fp_freq = np.power(10.0, freqs)
    P_t = model.derived_param["A_t"] * np.power(
        (10.0 ** freqs) / gp.f_piv, model.derived_param["nt"]
    )
    FS.solve_kernel(
        Nv, Phi, Phi_mid, S2, S2inv, j0, z0, P_t, ev_minus, fp_minus,
        fp_freq, int(defer_endpoint), idx_out.size, COL_STEP, H, Z_TAIL,
        *outputs, np.diff(Nv), sigma, PHASE_MAX, None,
        kink_index, kink_fraction, phi_re,
    )
    j_hi, g2c, w2c = _history_aggregate(freqs, j0, idx_out, outputs)
    nv = Nv.size
    g2 = np.empty(nv, dtype=np.float64)
    w2 = np.empty(nv, dtype=np.float64)
    FS.pchip_fine(idx_out.astype(np.float64), g2c, nv, g2)
    FS.pchip_fine(idx_out.astype(np.float64), w2c, nv, w2)
    Omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    DN_gw = gp.Neff0 * g2 * np.exp(
        2.0 * (f_hor - f_hor[-1]) * gp.ln10 + 2.0 * (Nv - Nv[-1])
    ) / Omega_nu
    return {
        "model": model,
        "Nv": Nv,
        "sigma": sigma,
        "f_hor": f_hor,
        "j0": j0,
        "deferred_mask": endpoint_defect > FS._WKB_TAIL_ENDPOINT_EPS_LIMIT,
        "idx_out": idx_out,
        "j_hi": j_hi,
        "history": outputs,
        "g2c": g2c,
        "w2c": w2c,
        "g2": g2,
        "w2": w2,
        "DN_gw": DN_gw,
    }


def _candidate_public_history(model, freqs, defer_endpoint):
    chart = prototype._post_event_chart(model)
    start_data = prototype._make_candidates(model, freqs, chart)
    start_indices, start_z, z_offsets, start_tau, direct_tail, start_xh, start_yh, pre_mask, N0 = start_data

    public_Nv, idx_out = _public_grid(model)
    sample_abs_N, sample_chart_indices = chart_history_indices(
        public_Nv, idx_out, chart["N_abs"], chart["N_star"]
    )
    sample_tau = sample_abs_N - chart["N_star"]
    pt = model.derived_param["A_t"] * np.power(
        (10.0 ** freqs) / gp.f_piv, model.derived_param["nt"]
    )
    fp_freq = np.power(10.0, freqs)
    ev_minus = np.exp(-chart["Nv"])
    fp_minus = np.exp(-chart["f_hor"] * gp.ln10)
    endpoint_z = start_z + chart["Phi"][-1] - chart["Phi"][start_indices]
    endpoint_defect = np.fromiter(
        (FS._endpoint_wkb_defect(float(chart["sigma"][-1]), float(v)) for v in endpoint_z),
        dtype=np.float64,
        count=endpoint_z.size,
    )
    deferred_mask = bool(defer_endpoint) & (endpoint_defect > FS._WKB_TAIL_ENDPOINT_EPS_LIMIT)
    if np.any(deferred_mask & direct_tail):
        raise NotImplementedError("deferred continuation from a pre-event direct-tail anchor is unsupported")

    outputs = [np.zeros((freqs.size, idx_out.size), dtype=np.float64) for _ in range(3)]
    anchor_z = np.full(freqs.size, np.nan)
    anchor_x = np.full(freqs.size, np.nan)
    anchor_y = np.full(freqs.size, np.nan)
    terminal_z = np.full(freqs.size, np.nan)
    terminal_x = np.full(freqs.size, np.nan)
    terminal_y = np.full(freqs.size, np.nan)
    _candidate_history_kernel(
        chart["Nv"], chart["Phi"], chart["Phi_mid"], chart["S2"],
        start_indices, start_z, z_offsets, start_tau, direct_tail,
        start_xh, start_yh, pt, ev_minus, fp_minus, fp_freq,
        np.diff(chart["Nv"]), chart["sigma"], PHASE_MAX, Z_TAIL,
        deferred_mask, sample_chart_indices, sample_tau, *outputs,
        anchor_z, anchor_x, anchor_y, terminal_z, terminal_x, terminal_y,
    )

    # Reconstruct the analytic sigma=1 prefix at the exact public coarse nodes.
    bessel_info = []
    for mode in np.flatnonzero(pre_mask):
        n0 = float(N0[mode])
        z0 = float((freqs[mode] - chart["f_re"]) * math.log(10.0) + 0.5 * (n0 - chart["N_star"]))
        c0, c1 = FS._bessel_sigma_one_coefficients(2.0 * math.exp(z0), 1.0)
        last_N = float(chart["N_star"])
        if direct_tail[mode]:
            last_N = float(chart["N_star"] + start_tau[mode])
        scale = math.exp(-0.5 * chart["N_star"])
        for slot, N_abs in enumerate(sample_abs_N):
            if n0 <= N_abs <= last_N:
                xh, yh = FS._bessel_sigma_one_fast_state(float(N_abs), n0, z0, c0, c1)
                zz = z0 + 0.5 * (float(N_abs) - n0)
                FS.assemble_main(
                    *outputs, int(mode), slot,
                    math.exp(-(float(N_abs) - chart["N_star"])),
                    scale * xh, scale * yh, zz, pt[mode],
                )
        bessel_info.append({
            "mode": int(mode), "N0": n0, "z0": z0,
            "coeff0": float(c0), "coeff1": float(c1),
            "last_N": last_N,
        })

    N_start = np.where(pre_mask, N0, chart["N_abs"][start_indices])
    j0_candidate = np.searchsorted(public_Nv, N_start).astype(np.int64)
    j_hi, g2c, w2c = _history_aggregate(freqs, j0_candidate, idx_out, outputs)
    nv = public_Nv.size
    g2 = np.empty(nv, dtype=np.float64)
    w2 = np.empty(nv, dtype=np.float64)
    FS.pchip_fine(idx_out.astype(np.float64), g2c, nv, g2)
    FS.pchip_fine(idx_out.astype(np.float64), w2c, nv, w2)
    # Public f_hor/Nv histories are restored on the old grid for this output.
    _, _, public_f_hor, _ = build_transition_grid(model, H)
    Omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    DN_gw = gp.Neff0 * g2 * np.exp(
        2.0 * (public_f_hor - public_f_hor[-1]) * gp.ln10
        + 2.0 * (public_Nv - public_Nv[-1])
    ) / Omega_nu
    return {
        "Nv": public_Nv,
        "f_hor": public_f_hor,
        "j0": j0_candidate,
        "idx_out": idx_out,
        "j_hi": j_hi,
        "history": outputs,
        "g2c": g2c,
        "w2c": w2c,
        "g2": g2,
        "w2": w2,
        "DN_gw": DN_gw,
        "deferred_mask": deferred_mask,
        "endpoint_defect": endpoint_defect,
        "anchor_states": {
            "z": anchor_z,
            "x": anchor_x,
            "y": anchor_y,
            "terminal_z": terminal_z,
            "terminal_x": terminal_x,
            "terminal_y": terminal_y,
        },
        "bessel_prefix": bessel_info,
        "chart": chart,
        "start_data": start_data,
    }


def _public_grid(model):
    d = model.derived_param
    N_inf = float(d["N_inf"])
    N_star = float(d["N_inf"] - d["N_re"])
    n_uniform = math.floor(N_inf / H) + 1
    Nv = np.unique(np.sort(np.concatenate((np.arange(0, n_uniform) * H, [N_star, N_inf]))))
    return Nv, public_history_indices(Nv)


def _difference_summary(baseline, candidate, freqs, idx_out):
    names = ("Ogw", "Oj", "Opgw")
    result = {}
    for channel, base, cand in zip(names, baseline["history"], candidate["history"]):
        delta = cand - base
        scale = max(float(np.max(np.abs(base))), float(np.max(np.abs(cand))), np.finfo(float).tiny)
        rel = np.abs(delta) / np.maximum(np.abs(base), scale * 1e-12)
        worst = np.unravel_index(int(np.argmax(rel)), rel.shape)
        mode, slot = map(int, worst)
        result[channel] = {
            "max_abs": float(np.max(np.abs(delta))),
            "max_relative_with_global_1e-12_floor": float(rel[worst]),
            "nonzero_difference_count": int(np.count_nonzero(delta)),
            "worst": {
                "mode": mode,
                "frequency": float(freqs[mode]),
                "frequency_offset_from_f_re": float(freqs[mode] - candidate["chart"]["f_re"]),
                "coarse_slot": slot,
                "original_idx_out": int(idx_out[slot]),
                "N": float(baseline["Nv"][idx_out[slot]]),
                "baseline": float(base[mode, slot]),
                "candidate": float(cand[mode, slot]),
                "delta": float(delta[mode, slot]),
            },
        }
    for name in ("g2c", "w2c", "g2", "w2", "DN_gw"):
        base = baseline[name]
        cand = candidate[name]
        delta = cand - base
        scale = max(float(np.max(np.abs(base))), float(np.max(np.abs(cand))), np.finfo(float).tiny)
        rel = np.abs(delta) / np.maximum(np.abs(base), scale * 1e-12)
        slot = int(np.argmax(rel))
        result[name] = {
            "max_abs": float(np.max(np.abs(delta))),
            "max_relative_with_global_1e-12_floor": float(rel[slot]),
            "worst_index": slot,
            "baseline": float(base[slot]),
            "candidate": float(cand[slot]),
            "terminal_relative": float(abs(cand[-1] - base[-1]) / max(abs(base[-1]), 1e-300)),
        }
    return result


def _run_reference(params, dn_re, dn_eff, freqs, z_tail):
    model = LCDM_SG(**params, DN_re=dn_re)
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = R.run_reference(
            model,
            dn_eff=dn_eff,
            freq_subset=freqs,
            self_consistent=False,
            z_tail=z_tail,
            rtol=1e-9,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
    elapsed = time.perf_counter() - started
    if not np.array_equal(np.asarray(result["freqs"], dtype=np.float64), freqs):
        raise RuntimeError(f"Reference z_tail={z_tail:g} changed the saved frequency grid")
    omega = np.asarray(result["Ogw"], dtype=np.float64) - np.asarray(result["Oj"], dtype=np.float64)
    return omega, [str(item.message) for item in caught], elapsed


def _load_profile():
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    params = dict(source["parameters"])
    dn_re = float(params.pop("DN_re"))
    dn_eff = float(source["baseline"]["root"])
    params["DN_eff"] = dn_eff
    freqs = np.asarray(source["baseline"]["frequencies_log10_hz"], dtype=np.float64)
    digest = hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest()
    if freqs.size != 183 or not np.all(np.diff(freqs) < 0.0):
        raise RuntimeError("expected the saved descending 183-frequency final-root support")
    if digest != source["baseline"]["frequency_sha256"]:
        raise RuntimeError("saved final-root frequency checksum mismatch")
    prior = json.loads(PRIOR_PROFILE.read_text(encoding="utf-8"))
    if _sha(Path(R.__file__)) != prior["source_hashes"]["reference"]:
        raise RuntimeError("live Reference source differs from the stored branch-profile Reference source")
    return source, params, dn_re, dn_eff, freqs


def run_validation() -> dict:
    source, params, dn_re, dn_eff, freqs = _load_profile()
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    baseline_model = LCDM_SG(**params, DN_re=dn_re)
    baseline = _public_history_arrays(baseline_model, freqs, defer_endpoint=True)
    candidate_model = LCDM_SG(**params, DN_re=dn_re)
    candidate = _candidate_public_history(candidate_model, freqs, defer_endpoint=True)

    if not np.array_equal(candidate["idx_out"], baseline["idx_out"]):
        raise RuntimeError("candidate public idx_out differs from Fast")
    if not np.array_equal(candidate["j0"], baseline["j0"]):
        raise RuntimeError("candidate public j0 differs from Fast")
    if not np.array_equal(candidate["j_hi"], baseline["j_hi"]):
        raise RuntimeError("candidate public j_hi differs from Fast")
    if not np.array_equal(candidate["deferred_mask"], baseline["deferred_mask"]):
        raise RuntimeError("candidate endpoint-deferred classification differs from full Fast")
    if candidate["history"][0].shape != baseline["history"][0].shape:
        raise RuntimeError("candidate public channel history shape differs from Fast")

    # Preserve the exact all-zero prefix before each public mode start.
    zero_prefix_failures = []
    for mode, j0 in enumerate(baseline["j0"]):
        before = np.flatnonzero(baseline["idx_out"] < j0)
        if before.size:
            for channel, values in zip(("Ogw", "Oj", "Opgw"), candidate["history"]):
                if np.any(values[mode, before] != 0.0):
                    zero_prefix_failures.append({"mode": int(mode), "channel": channel})

    base_spec = baseline["history"][0][:, -1] - baseline["history"][1][:, -1]
    cand_spec = candidate["history"][0][:, -1] - candidate["history"][1][:, -1]
    positive = (base_spec > 0.0) & (cand_spec > 0.0)
    delta_dex = np.full(freqs.size, np.nan)
    delta_dex[positive] = np.abs(np.log10(cand_spec[positive]) - np.log10(base_spec[positive]))
    omega_nu = gp.Omega_nh2 / candidate_model.derived_param["h"] ** 2
    baseline_dn_endpoint = float(baseline["DN_gw"][-1])
    candidate_dn_endpoint = float(candidate["DN_gw"][-1])
    root_relative = abs(candidate_dn_endpoint - baseline_dn_endpoint) / max(abs(baseline_dn_endpoint), 1e-300)
    self_consistency_relative = abs(candidate_dn_endpoint - dn_eff) / max(abs(dn_eff), 1e-300)

    ref5, warnings5, runtime5 = _run_reference(params, dn_re, dn_eff, freqs, 5.0)
    ref8, warnings8, runtime8 = _run_reference(params, dn_re, dn_eff, freqs, 8.0)
    ref_summaries = {}
    for ztail, reference in ((5.0, ref5), (8.0, ref8)):
        ref_positive = (reference > 0.0) & (cand_spec > 0.0)
        err = np.abs(np.log10(cand_spec[ref_positive]) - np.log10(reference[ref_positive]))
        ref_summaries[str(int(ztail))] = {
            "positive_mode_count": int(np.count_nonzero(ref_positive)),
            "max_abs_dex": float(np.max(err)),
            "p95_abs_dex": float(np.percentile(err, 95)),
            "warnings": warnings5 if ztail == 5.0 else warnings8,
            "runtime_s": float(runtime5 if ztail == 5.0 else runtime8),
            "map_relative_difference": float(
                abs(FS.integrate_frequency_quadrature(freqs, cand_spec, "pchip")
                    - FS.integrate_frequency_quadrature(freqs, reference, "pchip"))
                / max(abs(FS.integrate_frequency_quadrature(freqs, reference, "pchip")), 1e-300)
            ),
        }
    diff = _difference_summary(baseline, candidate, freqs, baseline["idx_out"])
    record = {
        "schema_version": 1,
        "status": "complete",
        "experiment": "cr0 event-chart public history reconstruction on the branch-valid final-root profile",
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        "input_artifact": str(INPUT.relative_to(ROOT)),
        "parameters": {**params, "DN_re": dn_re},
        "final_root_dn_eff": dn_eff,
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "history_layout": {
            "idx_out_count": int(baseline["idx_out"].size),
            "idx_out_sha256": hashlib.sha256(baseline["idx_out"].astype("<i8").tobytes()).hexdigest(),
            "idx_out_exact_match": bool(np.array_equal(candidate["idx_out"], baseline["idx_out"])),
            "j0_exact_match": bool(np.array_equal(candidate["j0"], baseline["j0"])),
            "j_hi_exact_match": bool(np.array_equal(candidate["j_hi"], baseline["j_hi"])),
            "history_shape": list(baseline["history"][0].shape),
            "zero_prefix_failures": zero_prefix_failures,
            "all_channel_values_finite": bool(all(np.isfinite(a).all() for a in candidate["history"])),
            "deferred_mode_count": int(np.count_nonzero(candidate["deferred_mask"])),
            "deferred_crossing_count": int(np.count_nonzero(
                candidate["deferred_mask"] & (candidate["anchor_states"]["terminal_z"] >= Z_TAIL)
            )),
            "deferred_no_crossing_count": int(np.count_nonzero(
                candidate["deferred_mask"] & (candidate["anchor_states"]["terminal_z"] < Z_TAIL)
            )),
            "deferred_mode_indices": np.flatnonzero(candidate["deferred_mask"]).tolist(),
            "deferred_direct_tail_count": int(np.count_nonzero(
                candidate["deferred_mask"] & candidate["start_data"][4]
            )),
        },
        "spectrum_and_map_gates": {
            "max_abs_spectrum_delta_dex_vs_full_fast": float(np.nanmax(delta_dex)),
            "spectrum_gate_1e-3_dex_pass": bool(np.nanmax(delta_dex) <= 1e-3),
            "candidate_dn_gw_endpoint": candidate_dn_endpoint,
            "baseline_dn_gw_endpoint": baseline_dn_endpoint,
            "candidate_vs_baseline_map_relative_difference": root_relative,
            "root_relative_gate_2e-4_pass": bool(root_relative <= 2e-4),
            "fixed_point_residual_relative_at_saved_baseline_dn_eff": self_consistency_relative,
            "fixed_point_residual_definition": "abs(candidate DN_gw endpoint map - saved baseline DN_eff root) / abs(saved baseline DN_eff root); the candidate outer root solver was not run",
            "fixed_point_residual_gate_2e-4_pass": bool(self_consistency_relative <= 2e-4),
            "candidate_self_consistent_root_solve_executed": False,
            "candidate_root_error_vs_independent_reference": None,
            "Omega_nu": float(omega_nu),
        },
        "reference_gates": {
            "tail_conventions_separate": True,
            "z_tail_5_frozen_default": ref_summaries["5"],
            "z_tail_8_deep_reference": ref_summaries["8"],
        },
        "_reference_payloads": {
            "z5": {"frequencies": freqs, "omega": ref5, "warnings": warnings5, "runtime_s": runtime5},
            "z8": {"frequencies": freqs, "omega": ref8, "warnings": warnings8, "runtime_s": runtime8},
        },
        "history_differences": diff,
        "candidate_states": {
            "missing_tail_anchor_semantics": "null means this mode did not cross z_tail and has no WKB anchor",
            "tail_anchor_z_range": [
                float(np.nanmin(candidate["anchor_states"]["z"])),
                float(np.nanmax(candidate["anchor_states"]["z"])),
            ],
            "tail_anchor_state_examples": [
                {
                    "mode": int(i),
                    "frequency": float(freqs[i]),
                    "z": _state_value(candidate["anchor_states"]["z"][i]),
                    "x": _state_value(candidate["anchor_states"]["x"][i]),
                    "y": _state_value(candidate["anchor_states"]["y"][i]),
                    "terminal_z": float(candidate["anchor_states"]["terminal_z"][i]),
                    "terminal_x": float(candidate["anchor_states"]["terminal_x"][i]),
                    "terminal_y": float(candidate["anchor_states"]["terminal_y"][i]),
                }
                for i in np.flatnonzero(candidate["deferred_mask"])[::max(1, int(np.count_nonzero(candidate["deferred_mask"])) // 6)]
            ],
        },
        "limitations": [
            "The baseline uses production full assembly and endpoint-deferred exact continuation; candidate public history is reconstructed in this independent harness.",
            "Internal baseline x/y states are not exposed by solve_kernel; state comparisons report candidate tail anchors/terminal states and compare the three resulting public channels column-by-column.",
            "The input is one branch-valid fixed background at a saved final-root DN_eff, not a self-consistent root solve over parameters.",
            "z_tail=5 frozen-default and z_tail=8 deep Reference results are separate conventions.",
        ],
        "resources": telemetry(workers=1, threads=2),
        "source_hashes": {
            "script": _sha(Path(__file__)),
            "prototype": _sha(Path(prototype.__file__)),
            "fast_sgwb": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "exact_background": _sha(ROOT / "stiffgwpy_fast" / "exact_background.py"),
            "reference": _sha(Path(R.__file__)),
        },
    }
    return record


def main() -> int:
    started = time.perf_counter()
    record = run_validation()
    record["elapsed_s"] = time.perf_counter() - started
    payloads = record.pop("_reference_payloads")
    prior, params, dn_re, dn_eff, freqs = _load_profile()
    ref_artifacts = {}
    for label, z_tail, path in (("z5", 5.0, OUTPUT_REF5), ("z8", 8.0, OUTPUT_REF8)):
        payload = payloads[label]
        ref_record = {
            "schema_version": 1,
            "status": "complete",
            "experiment": "independent Prüfer Reference spectrum on final-root frequency support",
            "tail_convention": f"z_tail={z_tail:g}",
            "input_artifact": str(INPUT.relative_to(ROOT)),
            "parameters": {**params, "DN_re": dn_re},
            "dn_eff": dn_eff,
            "frequency_count": int(freqs.size),
            "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
            "frequencies_log10_hz": payload["frequencies"].tolist(),
            "omega_gw": payload["omega"].tolist(),
            "runtime_s": payload["runtime_s"],
            "warnings": payload["warnings"],
            "reference_source_sha256": _sha(Path(R.__file__)),
        }
        path.write_text(
            json.dumps(ref_record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        ref_artifacts[label] = {
            "path": str(path.relative_to(ROOT)),
            "sha256": _sha(path),
        }
    record["reference_gates"]["separate_artifacts"] = ref_artifacts
    OUTPUT.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "frequency_count": record["frequency_count"],
        "history_layout": record["history_layout"],
        "spectrum_and_map_gates": record["spectrum_and_map_gates"],
        "reference_gates": record["reference_gates"],
        "artifact": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
