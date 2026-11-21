"""Paired warm Fast benchmark for adaptive goal/global seeds and mode reuse."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "benchmark_adaptive_support_reuse.json"
PARAMS = {
    "Omega_bh2": 13_890_136.71875,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": -0.022705078125,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
REPEATS = 7


def solve_case(label: str) -> dict[str, object]:
    model = LCDM_SG(**PARAMS)
    goal_builder = FA.goal_oriented_freqs
    kernel = FS.solve_kernel
    subset_solver = FS._solve_kernel_mode_subset
    mode_rows: list[int] = []

    if label == "global_seed_reuse":
        FA.goal_oriented_freqs = lambda m, freq_res=1.0, **kwargs: FA.grid_independent_freqs(
            m, freq_res, allow_log_amplitude=kwargs.get("allow_log_amplitude", False)
        )[0]
    if label == "goal_seed_recompute_all":
        FS._solve_kernel_mode_subset = lambda args, indices, **kwargs: FS.solve_kernel(
            *args, **kwargs
        )

    def count_modes(*args, **kwargs):
        mode_rows.append(int(len(args[5])))
        return kernel(*args, **kwargs)

    FS.solve_kernel = count_modes
    config = FS.resolve_config("fast", threads=2, freq_grid="adaptive")
    started = time.perf_counter()
    result = FS.SGWB_iter_fast(
        model,
        config=config,
        tol=1.0e-6,
        physical_guard=True,
        freq_grid_target=1.0e-6,
        freq_grid_max_points=1500,
    )
    elapsed = time.perf_counter() - started

    FS.solve_kernel = kernel
    FS._solve_kernel_mode_subset = subset_solver
    FA.goal_oriented_freqs = goal_builder
    if result is not model:
        raise RuntimeError((label, model.fast_failure_reason))
    return {
        "seconds": elapsed,
        "frequency_count": int(model.f.size),
        "propagated_mode_rows": int(sum(mode_rows)),
        "kernel_calls_mode_rows": mode_rows,
        "root": float(model.DN_gw[-1]),
    }


def main() -> None:
    labels = (
        "goal_seed_reuse",
        "goal_seed_recompute_all",
        "global_seed_reuse",
    )
    samples = {label: [] for label in labels}
    for repeat in range(REPEATS):
        order = labels[repeat % len(labels) :] + labels[: repeat % len(labels)]
        for label in order:
            row = solve_case(label)
            row["repeat"] = repeat
            samples[label].append(row)
    payload = {
        "status": "complete",
        "hardware": {"threads": 2},
        "case_params": PARAMS,
        "adaptive_settings": {"freq_grid_target": 1.0e-6, "freq_grid_max_points": 1500},
        "source_sha256": {
            name: hashlib.sha256(
                (ROOT / "stiffgwpy_fast" / name).read_bytes()
            ).hexdigest()
            for name in ("fast_sgwb.py", "reference.py", "exact_background.py", "_wkb_transfer.py")
        },
        "samples": samples,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
