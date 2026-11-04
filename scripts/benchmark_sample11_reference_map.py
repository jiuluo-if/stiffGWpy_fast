"""Run one independent tail-10 Reference map for the retained sample11 point.

This checks the Fast fixed point at an explicit DN_eff; it does not iterate
the Reference solver to its own root.
中文：在指定 DN_eff 上独立计算 sample11 的 tail-10 Reference 映射；本脚本不迭代 Reference 自洽根。
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import psutil

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = dict(
    r=0.000641847083065821,
    n_t=1.4400831295705445,
    kappa10=1.1693981508537732,
    T_re=3.101174976979448,
    DN_re=10.0,
    Omega_bh2=0.0223828,
    Omega_ch2=0.1201075,
    H0=67.32117,
    A_s=2.100549e-9,
    cr=0.0,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dn-eff", type=float, default=4.857519441699757)
    parser.add_argument("--z-tail", type=float, default=10.0)
    parser.add_argument("--rtol", type=float, default=1e-11)
    parser.add_argument(
        "--frequencies", type=Path,
        default=ROOT / "docs" / "sample11_reference_tail10_20261030.npz",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "docs" / "sample11_h02_reference_map_20261104.json",
    )
    parser.add_argument("--affinity", default=None,
                        help="optional comma-separated CPU IDs")
    args = parser.parse_args()
    process = psutil.Process()
    available = list(process.cpu_affinity())
    if args.affinity:
        selected = [int(value) for value in args.affinity.split(",")]
        if not selected or any(cpu not in available for cpu in selected):
            parser.error(f"affinity must be a nonempty subset of {available}")
        process.cpu_affinity(selected)
    else:
        selected = list(process.cpu_affinity())

    saved = np.load(args.frequencies)
    frequencies = np.asarray(saved["freqs"], dtype=float)
    model = LCDM_SG(**PARAMS, DN_eff=float(args.dn_eff))
    started = time.perf_counter()
    reference = run_reference(
        model, dn_eff=float(args.dn_eff), freq_subset=frequencies,
        z_tail=float(args.z_tail), rtol=float(args.rtol),
        self_consistent=False, physical_guard=False,
    )
    record = {
        "case": "sample11_h02_fast_root_reference_map",
        "source_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "cpu_affinity": selected,
        "seconds": time.perf_counter() - started,
        "dn_eff": float(args.dn_eff),
        "reference_dn_map": float(reference["DN_gw"]),
        "map_relative_residual": float(
            (reference["DN_gw"] - args.dn_eff) / args.dn_eff
        ),
        "z_tail": float(args.z_tail),
        "rtol": float(args.rtol),
        "frequency_count": int(len(reference["freqs"])),
        "quadrature_error": reference.get("quadrature_error"),
        "interpolation_error": reference.get("interpolation_error"),
    }
    args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
