"""Re-evaluate the independent fixed-background Reference map for res4.

中文：在指定 ``DN_eff`` 上用独立 Reference 求同一频率网格的固定背景 map，并记录尾部/插值/求积误差。
该积分较慢；参数和网格通过命令行显式绑定。
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
# 中文：先加入仓库根路径，再导入当前工作树内的 Reference 与模型。
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = dict(r=0.01, cr=0.0, n_t=2.0, T_re=1e4, kappa10=0.01,
              DN_re=10.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dn-eff", type=float, default=0.13373681817755848)
    parser.add_argument("--z-tail", type=float, default=10.0)
    parser.add_argument("--rtol", type=float, default=1e-11)
    parser.add_argument("--frequencies", type=Path,
                        default=ROOT / "docs" / "res4_reference_secant_root_tail10_20261104.npz")
    parser.add_argument("--output-prefix", type=Path,
                        default=ROOT / "docs" / "res4_reference_root_recheck_20261104")
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
    start = time.perf_counter()
    reference = run_reference(
        model, dn_eff=float(args.dn_eff), freq_subset=frequencies,
        z_tail=float(args.z_tail), rtol=float(args.rtol),
        self_consistent=False, physical_guard=False,
    )
    elapsed = time.perf_counter() - start
    record = {
        "case": "blue_res4_fixed_background_reference_map",
        "source_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "reference_source": str((ROOT / "stiffgwpy_fast" / "reference.py").resolve()),
        "cpu_affinity": selected,
        "seconds": elapsed,
        "dn_eff": float(args.dn_eff),
        "reference_dn_map": float(reference["DN_gw"]),
        "map_relative_residual": float(
            (reference["DN_gw"] - args.dn_eff) / args.dn_eff),
        "z_tail": float(args.z_tail),
        "rtol": float(args.rtol),
        "frequency_count": int(len(reference["freqs"])),
        "quadrature_error": reference.get("quadrature_error"),
        "interpolation_error": reference.get("interpolation_error"),
    }
    prefix = args.output_prefix
    prefix.with_suffix(".json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    np.savez_compressed(
        prefix.with_suffix(".npz"),
        freqs=np.asarray(reference["freqs"], dtype=float),
        reference_log=np.asarray(reference["log10OmegaGW"], dtype=float),
    )
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
