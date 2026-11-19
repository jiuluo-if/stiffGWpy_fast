"""Paired Fast A/B: exact matter transfer enabled versus disabled."""

# The scripts are run by path and add the repository root before importing the
# local checkout rather than whichever package version is installed.
# ruff: noqa: E402, I001

import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


RESEARCH = Path(__file__).resolve().parent
PARAMS = {
    "cr": 0.0,
    "T_re": 2.0e3,
    "kappa10": 1.0e-2,
    "n_t": 1.0,
    "r": 0.01,
    "A_s": 1.0e-15,
    "DN_re": 740.003,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def solve(exact: bool) -> tuple[dict, float]:
    original = FS.solve_kernel
    eligibility = []

    def wrapped(*args, **kwargs):
        call_args = list(args)
        nv = len(call_args[0])
        j0s, z0s, sigma = call_args[5], call_args[6], call_args[20]
        kink_index, kink_fraction, phi_re = call_args[23:26]
        eligible = 0
        if (
            kink_index >= 0
            and kink_index < nv - 1
            and math.isfinite(kink_fraction)
            and 0.0 <= kink_fraction < 1.0
        ):
            for j0, z0 in zip(j0s, z0s):
                j0 = int(j0)
                if j0 >= kink_index or sigma[j0] != 1.0 or sigma[kink_index] != 1.0:
                    continue
                d_n = call_args[0][kink_index] - call_args[0][j0]
                d_phi = call_args[1][kink_index] - call_args[1][j0]
                z_re = z0 + phi_re - call_args[1][j0]
                if (
                    abs(d_phi - 0.5 * d_n)
                    <= 128.0 * np.finfo(np.float64).eps * max(1.0, abs(d_n), abs(d_phi))
                    and z_re < call_args[15]
                ):
                    eligible += 1
        eligibility.append({
            "count": eligible,
            "kink_index": int(kink_index),
            "kink_fraction": float(kink_fraction),
            "sigma_start_min": float(np.min(sigma[j0s.astype(np.int64)])),
            "sigma_kink": float(sigma[kink_index]) if 0 <= kink_index < nv else None,
            "min_z_event": float(np.min(z0s + phi_re - call_args[1][j0s.astype(np.int64)])),
            "max_z_event": float(np.max(z0s + phi_re - call_args[1][j0s.astype(np.int64)])),
        })
        if not exact:
            call_args[24] = float("nan")  # retain event metadata; bypass only MD shortcut
        return original(*call_args, **kwargs)

    FS.solve_kernel = wrapped
    try:
        model = LCDM_SG(**PARAMS)
        start = time.perf_counter()
        ok = FS.SGWB_iter_fast(model, config=FS.resolve_config("fast", threads=2), tol=1e-6)
        elapsed = time.perf_counter() - start
        if not ok:
            raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
        return {
            "root": float(model.DN_gw[-1]),
            "frequency": np.asarray(model.f, dtype=float).tolist(),
            "logspec": np.asarray(model.log10OmegaGW, dtype=float).tolist(),
            "eligibility": eligibility,
        }, elapsed
    finally:
        FS.solve_kernel = original


def main() -> None:
    solve(True)
    solve(False)
    rows = []
    for _ in range(4):
        for label, exact in (("exact", True), ("baseline", False)):
            result, elapsed = solve(exact)
            rows.append({"variant": label, "seconds": elapsed, **result})
    exact_rows = [row for row in rows if row["variant"] == "exact"]
    base_rows = [row for row in rows if row["variant"] == "baseline"]
    frequency_grid_differences = [
        float(np.max(np.abs(np.asarray(a["frequency"]) - np.asarray(b["frequency"]))))
        if len(a["frequency"]) == len(b["frequency"])
        else float("inf")
        for a, b in zip(exact_rows, base_rows)
    ]
    spectrum_differences = []
    for a, b in zip(exact_rows, base_rows):
        af = np.asarray(a["frequency"])
        bf = np.asarray(b["frequency"])
        if len(af) != len(bf):
            spectrum_differences.append(float("inf"))
            continue
        b_interp = np.interp(af, bf[::-1], np.asarray(b["logspec"])[::-1])
        spectrum_differences.append(float(np.max(np.abs(np.asarray(a["logspec"]) - b_interp))))
    record = {
        "parameters": PARAMS,
        "source_head": __import__("subprocess").check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256": sha256(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "pairs": rows,
        "max_frequency_grid_difference": max(frequency_grid_differences),
        "max_abs_root_difference": max(
            abs(a["root"] - b["root"]) for a, b in zip(exact_rows, base_rows)
        ),
        "max_spectrum_difference_dex_after_frequency_interpolation": max(spectrum_differences),
    }
    output = RESEARCH / "exact_matter_transfer_branch_ab.json"
    output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in record.items() if k != "pairs"}, indent=2))
    print("timings", [(r["variant"], r["seconds"]) for r in rows])


if __name__ == "__main__":
    main()
