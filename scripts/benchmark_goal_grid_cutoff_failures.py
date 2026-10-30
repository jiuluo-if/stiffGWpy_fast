"""Reproduce the legacy goal-grid cutoff runaway on four stiff cases.

The legacy grid logic is reconstructed locally as ``legacy_goal_grid``; the
production grid and solver remain untouched. The comparison isolates the
feature-reserve guard added to :mod:`stiffgwpy_fast.freq_adaptive`.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import psutil

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    2: dict(r=0.5302241362263932, n_t=1.0622501865308032,
            kappa10=132.3355544029978, T_re=0.5870367189712415),
    5: dict(r=0.04220983994484846, n_t=1.2684771843868714,
            kappa10=0.12474949452626993, T_re=8.191629936517634),
    11: dict(r=0.000641847083065821, n_t=1.4400831295705445,
             kappa10=1.1693981508537732, T_re=3.101174976979448),
    18: dict(r=0.0017368769573485178, n_t=1.3996377339218413,
             kappa10=194.7347634468269, T_re=9.068561481578515),
}
COMMON = dict(DN_re=10.0, Omega_bh2=0.0223828, Omega_ch2=0.1201075,
              H0=67.32117, A_s=2.100549e-9, cr=0.0, DN_eff=0.0)


def legacy_goal_grid(m, freq_res=1.0, seed_n=64, max_points=120,
                     eval_freqs=None, feature_half_width=1.0, feature_n=17):
    """Pre-fix goal grid: always create the feature reserve, even if disjoint."""
    seed_n = max(2, int(seed_n))
    max_points = max(seed_n, int(max_points))
    feature_n = max(3, int(feature_n))
    base, fmin, fmax = FA.grid_independent_freqs(m, freq_res)
    if base.size <= seed_n:
        seed = base
    else:
        idx = np.rint(np.linspace(0, base.size - 1, seed_n)).astype(int)
        seed = base[np.unique(idx)]
    f_re = float(getattr(m, "f_re", 0.5 * (fmin + fmax)))
    lo = max(float(fmin), f_re - float(feature_half_width))
    hi = min(float(fmax), f_re + float(feature_half_width))
    feature = np.linspace(hi, lo, feature_n)
    if eval_freqs is None:
        eval = np.empty(0, dtype=float)
    else:
        eval = np.asarray(eval_freqs, dtype=float).reshape(-1)
        eval = eval[np.isfinite(eval)]
        eval = eval[(eval >= fmin) & (eval <= fmax)]
    mandatory = np.unique(np.concatenate((feature, eval)))
    if mandatory.size >= max_points:
        selected = np.unique(np.concatenate((eval, feature)))
        if selected.size > max_points and eval.size < max_points:
            selected = np.unique(np.concatenate((eval, feature[:max_points - eval.size])))
        return np.sort(selected)[::-1]
    optional = seed[(seed < lo) | (seed > hi)]
    room = max_points - mandatory.size
    if optional.size > room:
        idx = np.rint(np.linspace(0, optional.size - 1, room)).astype(int)
        optional = optional[np.unique(idx)]
    return np.sort(np.unique(np.concatenate((mandatory, optional))))[::-1]


def run_case(sample_id: int, use_legacy_grid: bool) -> dict:
    params = {**COMMON, **CASES[sample_id]}
    model = LCDM_SG(**params)
    grids, kernel = [], []
    original_grid = FA.goal_oriented_freqs
    original_kernel = FS.solve_kernel

    def grid_builder(*args, **kwargs):
        builder = legacy_goal_grid if use_legacy_grid else original_grid
        freqs = np.asarray(builder(*args, **kwargs), dtype=float)
        d = args[0].derived_param if args else model.derived_param
        if not args or not hasattr(args[0], "derived_param"):
            d = model.derived_param
        nt = float(d["nt"])
        at = float(d["A_t"])
        cutoff = (-np.log10(at) / nt + np.log10(FS.gp.f_piv)) if nt > 0 else np.inf
        grids.append({
            "count": int(freqs.size),
            "max_log10f": float(np.max(freqs)),
            "amplitude_fmax": float(cutoff),
            "nodes_above_cutoff": int(np.count_nonzero(freqs > cutoff)),
            "max_overshoot_dex": float(max(0.0, np.max(freqs) - cutoff)),
        })
        return freqs

    def kernel_probe(*args, **kwargs):
        before_s2 = np.asarray(args[3])
        bad_background = not np.isfinite(before_s2).all()
        result = original_kernel(*args, **kwargs)
        ogw, oj, opgw = (np.asarray(args[i]) for i in (16, 17, 18))
        kernel.append({
            "dn_eff_input": float(model.cosmo_param["DN_eff"]),
            "background_s2_nonfinite": int(before_s2.size - np.count_nonzero(np.isfinite(before_s2))),
            "output_nonfinite": {
                name: int(arr.size - np.count_nonzero(np.isfinite(arr)))
                for name, arr in (("Ogw", ogw), ("Oj", oj), ("Opgw", opgw))
            },
            "first_nonfinite_layer": "background_S2" if bad_background else None,
        })
        return result

    FA.goal_oriented_freqs = grid_builder
    FS.solve_kernel = kernel_probe
    try:
        config = FS.resolve_config("fast", threads=2, z_tail=5.0, freq_grid="goal")
        started = time.perf_counter()
        result = FS.SGWB_iter_fast(model, config=config, kink_split=True,
                                   freq_res=1.0, tol=1e-6)
        elapsed = time.perf_counter() - started
    finally:
        FA.goal_oriented_freqs = original_grid
        FS.solve_kernel = original_kernel

    status = "failure" if result is not model else "finite"
    reason = getattr(model, "fast_failure_reason", None)
    if result is model and not np.isfinite(np.asarray(model.log10OmegaGW)).all():
        status = "nonfinite_output"
    first_bad = next((item for item in kernel if item["first_nonfinite_layer"]), None)
    return {
        "sample_id": sample_id,
        "grid_path": "legacy_unconditional_reserve" if use_legacy_grid else "current_cutoff_guard",
        "status": status,
        "failure_reason": reason,
        "first_nonfinite_layer": first_bad["first_nonfinite_layer"] if first_bad else None,
        "first_nonfinite_dn_eff": first_bad["dn_eff_input"] if first_bad else None,
        "dn_gw": float(model.DN_gw[-1]) if result is model else None,
        "frequency_grid_builds": grids,
        "kernel_calls": kernel,
        "runtime_s": elapsed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs/stiff_goal_grid_cutoff_replay_20261030.json")
    args = parser.parse_args()
    process = psutil.Process()
    affinity = process.cpu_affinity()
    if len(affinity) > 2:
        process.cpu_affinity(affinity[4:6])
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    records = []
    for sample_id in CASES:
        pair = {
            "sample_id": sample_id,
            "legacy": run_case(sample_id, True),
            "current": run_case(sample_id, False),
        }
        records.append(pair)
        print(json.dumps(pair, ensure_ascii=False), flush=True)
    payload = {
        "legacy_grid_commit": "9f1a70a9885fea8f03db9e43d948a2a7daed5f19",
        "cutoff_fix_commit": "d363f9b2a497d15628f947234ec6fe7790c4c710",
        "current_code_commit": __import__("subprocess").check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_file": FS.__file__,
        "resources": {"threads": 2, "cpu_affinity": process.cpu_affinity()},
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
