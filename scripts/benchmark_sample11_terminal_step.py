"""Reproduce the terminal-step correction on the retained random sample11 case."""
from __future__ import annotations

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

PARAMS = {
    "r": 0.000641847083065821,
    "n_t": 1.4400831295705445,
    "kappa10": 1.1693981508537732,
    "T_re": 3.101174976979448,
    "DN_re": 10.0,
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "cr": 0.0,
    "DN_eff": 0.0,
}


def run(z_tail: float, legacy_step: bool, frequencies: np.ndarray) -> dict:
    model = LCDM_SG(**PARAMS)
    config = FS.resolve_config(
        "fast", freq_grid="grid_independent", h=0.005,
        z_tail=z_tail, threads=2,
    )
    original = FS.solve_kernel
    if legacy_step:
        def use_legacy_terminal_step(*args, **kwargs):
            call_args = list(args)
            call_args[19] = None
            if len(call_args) > 26:
                call_args[26] = 0.0
            return original(*call_args, **kwargs)

        FS.solve_kernel = use_legacy_terminal_step
    try:
        started = time.perf_counter()
        result = FS.SGWB_iter_fast(
            model,
            config=config,
            freq_res=1.0,
            tol=1e-6,
            kink_split=True,
            eval_freqs=frequencies,
        )
        elapsed = time.perf_counter() - started
    finally:
        FS.solve_kernel = original
    if result is not model:
        raise RuntimeError(model.fast_failure_reason or "Fast did not converge")

    order = np.argsort(np.asarray(model.f, dtype=float))
    ordered_freqs = np.asarray(model.f, dtype=float)[order]
    indexes = np.searchsorted(ordered_freqs, frequencies)
    if np.any(indexes >= ordered_freqs.size) or not np.array_equal(
        ordered_freqs[indexes], frequencies
    ):
        raise RuntimeError("Fast omitted a required Reference frequency")
    return {
        "runtime_s": elapsed,
        "dn_eff": float(model.cosmo_param["DN_eff"]),
        "dn_gw": float(model.DN_gw[-1]),
        "dn_residual": float(model.dn_converged_delta),
        "log_spectrum": np.asarray(model.log10OmegaGW, dtype=float)[order][indexes],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "docs/sample11_terminal_step_20261030.json",
    )
    args = parser.parse_args()
    reference_path = ROOT / "docs/sample11_reference_tail10_spectrum_20261030.json"
    reference_meta_path = ROOT / "docs/sample11_reference_tail10_20261030.json"
    second_map_path = ROOT / "docs/sample11_reference_second_map_20261030.json"
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    frequencies = np.asarray(reference["frequencies"], dtype=float)
    reference_log = np.asarray(reference["reference_log"], dtype=float)
    reference_meta = json.loads(reference_meta_path.read_text(encoding="utf-8"))
    reference_root = json.loads(second_map_path.read_text(encoding="utf-8"))
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)

    results = []
    for z_tail in (5.0, 8.0):
        legacy = run(z_tail, True, frequencies)
        corrected = run(z_tail, False, frequencies)
        legacy_error = np.abs(legacy["log_spectrum"] - reference_log)
        corrected_error = np.abs(corrected["log_spectrum"] - reference_log)
        results.append(
            {
                "z_tail": z_tail,
                "legacy": {
                    key: value for key, value in legacy.items()
                    if key != "log_spectrum"
                },
                "corrected": {
                    key: value for key, value in corrected.items()
                    if key != "log_spectrum"
                },
                "legacy_spectrum_median_p95_max_dex": np.percentile(
                    legacy_error, [50, 95, 100]
                ).tolist(),
                "corrected_spectrum_median_p95_max_dex": np.percentile(
                    corrected_error, [50, 95, 100]
                ).tolist(),
                "legacy_modes_over_5e-4": int(np.count_nonzero(legacy_error > 5e-4)),
                "corrected_modes_over_5e-4": int(np.count_nonzero(corrected_error > 5e-4)),
                "corrected_dn_root_rel_vs_reference": abs(
                    corrected["dn_gw"] - reference_root["reference_dn_map"]
                ) / abs(reference_root["reference_dn_map"]),
            }
        )

    payload = {
        "source_commit": __import__("subprocess").check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_file": FS.__file__,
        "case": "random_solver_stiff_sample11",
        "params": PARAMS,
        "frequency_count": int(frequencies.size),
        "reference_background_dn": reference_meta["dn_eff"],
        "reference_fixed_point_dn": reference_root["reference_dn_map"],
        "reference_map_residual": reference_root["map_relative_residual"],
        "reference_spectrum_source": str(reference_path.relative_to(ROOT)),
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
