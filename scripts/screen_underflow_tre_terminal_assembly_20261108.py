"""Screen the underflow Fast branch across a T_re line for zero transfers."""

from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BASE_PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 4.985023e-317,
    "n_t": 17.375,
    "cr": 0.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}


def main() -> None:
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    tre_values = np.logspace(1.0, math.log10(2000.0), 41)
    rows = []
    for tre in tre_values:
        model = LCDM_SG(**(BASE_PARAMS | {"T_re": float(tre)}))
        started = time.perf_counter()
        result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=1e-6,
            z_tail=5.0,
            freq_res=1.0,
            threads=2,
            accuracy_mode="fast",
            physical_guard=False,
            freq_grid="construct",
        )
        rows.append(
            {
                "T_re_GeV": float(tre),
                "converged": result is model and bool(model.SGWB_converge),
                "failure_reason": model.fast_failure_reason,
                "root_DN_gw": (
                    float(model.DN_gw[-1])
                    if result is model and model.SGWB_converge
                    else None
                ),
                "frequency_count": int(len(model.f)),
                "runtime_s": time.perf_counter() - started,
            }
        )
    source = ROOT / "stiffgwpy_fast" / "fast_sgwb.py"
    payload = {
        "schema_version": 1,
        "experiment": "underflow Fast T_re line screen after terminal coarse-sample fix",
        "fast_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "fixed_params": BASE_PARAMS,
        "T_re_grid_GeV": "41-point log-spaced [10, 2000]",
        "threads": 2,
        "points": rows,
        "scope_limit": "Fast convergence/failure screen only; no Reference accuracy claim.",
    }
    out = ROOT / "docs" / "underflow_tre_terminal_assembly_scan_20261108.json"
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
