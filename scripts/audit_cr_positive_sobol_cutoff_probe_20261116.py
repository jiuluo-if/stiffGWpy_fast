"""Trace the expected warning from the invalid-amplitude cutoff probe."""

from __future__ import annotations

import hashlib
import io
import json
import math
import platform
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_cutoff_holdout_20261116.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_cutoff_probe_audit_20261116.json"
SAFE_PROBE_R = 1e-2


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(params: dict[str, float]) -> tuple[LCDM_SG, dict[str, object], str]:
    stream = io.StringIO()
    with redirect_stdout(stream), redirect_stderr(stream):
        model = LCDM_SG(**params)
        derived = model.derived_param
    return model, derived, stream.getvalue()


def main() -> int:
    preflight = json.loads(INPUT.read_text(encoding="utf-8"))
    if preflight["source_sha256"]["stiff_SGWB.py"] != sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"):
        raise RuntimeError("Sobol preflight source differs from current stiff_SGWB.py")
    base = {key: float(value) for key, value in preflight["base_parameters"].items()}
    r_min = float(np.nextafter(0.0, 1.0))
    rows = []
    for sample in preflight["samples"]:
        sobol_index = int(sample["sobol_index"])
        t_re = float(sample["parameters"]["T_re"])
        kappa = float(sample["parameters"]["kappa10"])
        common = {**base, "T_re": t_re, "kappa10": kappa}

        min_model, min_derived, min_stdout = build({**common, "r": r_min})
        safe_model, safe_derived, safe_stdout = build({**common, "r": SAFE_PROBE_R})
        log_r_cut = (
            -4.0 * (gp.N_10 - float(safe_derived["N_re"]))
            - 4.0 * math.log(gp.M_Pl)
            - math.log(45.0 / 2.0)
            - math.log(float(common["A_s"]))
            + math.log(safe_model.rhorad_re + safe_model.rhostiff_re)
        )
        r_cut = math.exp(log_r_cut)
        logged_cutoff = float(sample["cutoff"]["r_cutoff"])
        background_equal = {
            "N_re_bitwise_equal": bool(float(min_derived["N_re"]) == float(safe_derived["N_re"])),
            "rhorad_re_bitwise_equal": bool(min_model.rhorad_re == safe_model.rhorad_re),
            "rhostiff_re_bitwise_equal": bool(min_model.rhostiff_re == safe_model.rhostiff_re),
        }
        if not all(background_equal.values()):
            raise RuntimeError(f"Sobol point {sobol_index}: cutoff background changed with r")
        if not math.isclose(r_cut, logged_cutoff, rel_tol=2e-14, abs_tol=0.0):
            raise RuntimeError(f"Sobol point {sobol_index}: safe-probe cutoff differs from preflight cutoff")
        safe_delta = float(safe_derived["N_inf"] - safe_derived["N_re"])
        if safe_delta <= 0.0:
            raise RuntimeError(f"Sobol point {sobol_index}: safe probe is not cutoff-valid")

        rows.append({
            "sobol_index": sobol_index,
            "T_re": t_re,
            "kappa10": kappa,
            "preflight_sample_status": sample["status"],
            "r_min_probe": {
                "r": r_min,
                "stored_A_t": float(min_derived["A_t"]),
                "N_inf": min_derived["N_inf"],
                "captured_stdout": min_stdout.splitlines(),
                "fast_or_reference_propagated": False,
            },
            "safe_probe": {
                "r": SAFE_PROBE_R,
                "N_inf": float(safe_derived["N_inf"]),
                "N_re": float(safe_derived["N_re"]),
                "Delta_N": safe_delta,
                "captured_stdout": safe_stdout.splitlines(),
                "analytic_r_cut": r_cut,
                "preflight_r_cut": logged_cutoff,
                "relative_cutoff_difference": r_cut / logged_cutoff - 1.0,
                "background_r_invariance": background_equal,
            },
        })

    output = {
        "schema_version": 1,
        "status": "complete_sobol_cutoff_probe_warning_audit",
        "scope": "Reconstruct the preflight's r=minsub cutoff-only probe and compare its stdout warning with a safe r=.01 background probe. Verify that N_re and radiation/stiff densities are r-independent and that the exact analytic cutoff matches every Sobol record. No invalid point is propagated by Fast or Reference.",
        "sobol_preflight_input": str(INPUT.relative_to(ROOT)).replace("\\", "/"),
        "r_min": r_min,
        "safe_probe_r": SAFE_PROBE_R,
        "samples": rows,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "source_sha256": {
            "stiff_SGWB.py": sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "LCDM_stiff_Neff.py": sha(ROOT / "stiffgwpy_fast" / "LCDM_stiff_Neff.py"),
            "global_param.py": sha(ROOT / "stiffgwpy_fast" / "global_param.py"),
            "sobol_preflight_script": sha(ROOT / "scripts" / "screen_cr_positive_sobol_cutoff_holdout_20261116.py"),
            "sobol_preflight_input": sha(INPUT),
            "audit_script": sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "sample_count": len(rows),
        "r_min_warning_samples": sum(bool(row["r_min_probe"]["captured_stdout"]) for row in rows),
        "safe_probe_warning_samples": sum(bool(row["safe_probe"]["captured_stdout"]) for row in rows),
        "cutoff_formula_mismatches": sum(abs(row["safe_probe"]["relative_cutoff_difference"]) > 2e-14 for row in rows),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
