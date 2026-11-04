"""Measure one fresh-process Reference mode with lazy spline compilation."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = dict(r=0.01, cr=0.0, n_t=2.0, T_re=1e4, kappa10=0.01,
              DN_re=10.0)
DN_EFF = 0.133732083004284


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freq", type=float, required=True)
    parser.add_argument("--rtol", type=float, default=1e-11)
    parser.add_argument("--z-tail", type=float, default=10.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
    signatures_before = len(REF._cubic_spline_eval.signatures)
    source_path = ROOT / "stiffgwpy_fast" / "reference.py"
    start = time.perf_counter()
    result = REF.solve_reference_mode(
        model, args.freq, DN_EFF, z_tail=args.z_tail, rtol=args.rtol,
        backend="cartesian")
    runtime = time.perf_counter() - start
    payload = {
        "schema_version": 1,
        "experiment": "reference_spline_single_mode_cold_dispatch",
        "source_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "source_script_sha256": hashlib.sha256(
            Path(__file__).read_bytes()).hexdigest(),
        "reference_module_sha256": hashlib.sha256(
            source_path.read_bytes()).hexdigest(),
        "parameters": {**PARAMS, "DN_eff": DN_EFF},
        "freq_log10_hz": args.freq,
        "z_tail": args.z_tail,
        "rtol": args.rtol,
        "atol": [1e-12, 1e-22, 1e-22],
        "runtime_s": runtime,
        "jit_signatures_before": signatures_before,
        "jit_signatures_after": len(REF._cubic_spline_eval.signatures),
        "nfev": int(result["n_steps"]),
        "used_tail": bool(result["used_tail"]),
        "log10_omega": math.log10(
            result["Ogw_today"] - result["Oj_today"]),
        "resources": telemetry(workers=1, threads=1),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
