import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    "default": dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    "lowT": dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    "highT": dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    "stiff": dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
}


def run(case, legacy, freqs):
    model = LCDM_SG(**CASES[case])
    original = FS.solve_kernel
    if legacy:

        def old_step(*args, **kwargs):
            args = list(args)
            args[19] = None
            if len(args) > 26:
                args[26] = 0.0
            return original(*args, **kwargs)

        FS.solve_kernel = old_step
    try:
        start = time.perf_counter()
        result = FS.SGWB_iter_fast(model, kink_split=True, eval_freqs=freqs)
        seconds = time.perf_counter() - start
    finally:
        FS.solve_kernel = original
    if result is not model:
        raise RuntimeError(f"{case}: {model.fast_failure_reason}")
    order = np.argsort(np.asarray(model.f, dtype=float))
    ff = np.asarray(model.f)[order]
    ix = np.searchsorted(ff, freqs)
    if np.any(ix >= ff.size) or not np.array_equal(ff[ix], freqs):
        raise RuntimeError("oracle frequency alignment failed")
    return model, np.asarray(model.log10OmegaGW, dtype=float)[order][ix], seconds


def main():
    parser = argparse.ArgumentParser(
        description="Compare legacy and actual terminal-step kink propagation against Oracle C."
    )
    parser.add_argument("--output", type=Path, default=ROOT / "docs/kink_terminal_step_fix_oracle_c_20261030.json")
    args = parser.parse_args()
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    summary = {
        "source_file": FS.__file__,
        "source_commit": __import__("subprocess")
        .check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True)
        .strip(),
        "cases": {},
    }
    for case in CASES:
        oracle = json.loads((ROOT / f"docs/oracle_c_wkb_{case}.json").read_text(encoding="utf-8"))
        freqs = np.asarray(oracle["frequencies"], dtype=float)
        ref = []
        use = []
        for row in oracle["rows"]:
            selected = row["wkb_today"] if row.get("wkb_today") is not None else row["deep_today"]
            value = float(selected["Ogw_today"] - selected["Oj_today"])
            ref.append(value)
            use.append("wkb" if row.get("wkb_today") is not None else "deep")
        ref = np.asarray(ref, dtype=float)
        positive = ref > 0.0
        if not np.any(positive):
            raise RuntimeError(f"{case}: no positive oracle power bins")
        old, old_log, old_s = run(case, True, freqs)
        new, new_log, new_s = run(case, False, freqs)
        ref_log = np.log10(ref[positive])
        old_err = np.abs(old_log[positive] - ref_log)
        new_err = np.abs(new_log[positive] - ref_log)
        positive_freqs = freqs[positive]
        valid_ref = ref[positive]
        old_idx = int(np.argmax(old_err))
        new_idx = int(np.argmax(new_err))
        summary["cases"][case] = {
            "oracle_source_commit": oracle["generated_commit"],
            "legacy_dn_gw": float(old.DN_gw[-1]),
            "corrected_dn_gw": float(new.DN_gw[-1]),
            "oracle_dn_wkb": float(oracle["summary"]["DN_wkb_corrected"]),
            "legacy_dn_rel": abs(float(old.DN_gw[-1]) - float(oracle["summary"]["DN_wkb_corrected"]))
            / max(abs(float(oracle["summary"]["DN_wkb_corrected"])), 1e-300),
            "corrected_dn_rel": abs(float(new.DN_gw[-1]) - float(oracle["summary"]["DN_wkb_corrected"]))
            / max(abs(float(oracle["summary"]["DN_wkb_corrected"])), 1e-300),
            "legacy_spectrum_median_p95_max_dex": np.percentile(old_err, [50, 95, 100]).tolist(),
            "corrected_spectrum_median_p95_max_dex": np.percentile(new_err, [50, 95, 100]).tolist(),
            "legacy_modes_over_1e3": int(np.count_nonzero(old_err > 1e-3)),
            "corrected_modes_over_1e3": int(np.count_nonzero(new_err > 1e-3)),
            "compared_positive_modes": int(np.count_nonzero(positive)),
            "excluded_nonpositive_reference_modes": int(np.count_nonzero(~positive)),
            "legacy_worst_mode": {
                "f": float(positive_freqs[old_idx]),
                "reference_omega": float(valid_ref[old_idx]),
                "fast_log10omega": float(old_log[positive][old_idx]),
                "used_reference_tail": use[old_idx],
            },
            "corrected_worst_mode": {
                "f": float(positive_freqs[new_idx]),
                "reference_omega": float(valid_ref[new_idx]),
                "fast_log10omega": float(new_log[positive][new_idx]),
                "used_reference_tail": use[new_idx],
            },
            "legacy_runtime_s": old_s,
            "corrected_runtime_s": new_s,
            "reference_tail_sources": {name: use.count(name) for name in set(use)},
        }
        print(case, json.dumps(summary["cases"][case], separators=(",", ":")), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
