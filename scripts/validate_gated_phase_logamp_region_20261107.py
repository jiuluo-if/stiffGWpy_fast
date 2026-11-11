"""Validate a sampled two-parameter underflow region: 验证二维下溢采样区域。"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _array_sha(values) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def _compare_point(reference_dir: Path, r_name: str, nt_name: str) -> dict:
    stem = f"fast_log_amplitude_grid_{r_name}_{nt_name}_seed20261107"
    z5_path = reference_dir / f"{stem}_z5.json"
    z8_path = reference_dir / f"{stem}_z8.json"
    z5 = json.loads(z5_path.read_text(encoding="utf-8"))
    z8 = json.loads(z8_path.read_text(encoding="utf-8"))
    params = z8["params"]
    rows5, rows8 = z5["frequency_rows"], z8["frequency_rows"]
    if len(rows5) != len(rows8):
        raise RuntimeError(f"z5/z8 frequency-count mismatch for {stem}")
    freq5 = np.asarray([row["log10_hz"] for row in rows5], dtype=np.float64)
    freq8 = np.asarray([row["log10_hz"] for row in rows8], dtype=np.float64)
    if not np.array_equal(freq5, freq8):
        raise RuntimeError(f"z5/z8 frequencies differ for {stem}")
    ref_sha5 = z5["fast"]["frequency_sha256"]
    ref_sha8 = z8["fast"]["frequency_sha256"]
    if ref_sha5 != ref_sha8 or _array_sha(freq8) != ref_sha8:
        raise RuntimeError(f"reference frequency SHA mismatch for {stem}")

    model = LCDM_SG(**params)
    started = time.perf_counter()
    model.SGWB_iter(
        engine="fast",
        fallback=False,
        tol=1e-6,
        z_tail=5.0,
        freq_res=1.0,
        threads=2,
        accuracy_mode="fast",
        physical_guard=False,
    )
    elapsed = time.perf_counter() - started
    candidate_freq_sha = _array_sha(model.f)
    if candidate_freq_sha != ref_sha8:
        raise RuntimeError(f"candidate/reference frequencies differ for {stem}")
    candidate = np.asarray(model.log10OmegaGW, dtype=np.float64)
    refs = {
        "z5": np.asarray([row["reference_log10_OmegaGW"] for row in rows5]),
        "z8": np.asarray([row["reference_log10_OmegaGW"] for row in rows8]),
    }
    errors = {tail: np.abs(candidate - values) for tail, values in refs.items()}
    roots = {
        "z5": float(z5["scaled_reference"]["estimated_self_consistent_root_DN_gw"]),
        "z8": float(z8["scaled_reference"]["estimated_self_consistent_root_DN_gw"]),
    }
    root = float(model.DN_gw[-1])
    return {
        "point": f"{r_name}_{nt_name}",
        "params": params,
        "frequency_count": int(freq8.size),
        "frequency_sha256": candidate_freq_sha,
        "runtime_s": elapsed,
        "fast_root_DN_gw": root,
        "reference_root_z5": roots["z5"],
        "reference_root_z8": roots["z8"],
        "root_relative_error_z5": (root - roots["z5"]) / roots["z5"],
        "root_relative_error_z8": (root - roots["z8"]) / roots["z8"],
        "max_spectrum_error_dex_z5": float(errors["z5"].max()),
        "p95_spectrum_error_dex_z5": float(np.percentile(errors["z5"], 95)),
        "max_spectrum_error_dex_z8": float(errors["z8"].max()),
        "p95_spectrum_error_dex_z8": float(np.percentile(errors["z8"], 95)),
        "worst_z8_log10_hz": float(freq8[int(np.argmax(errors["z8"]))]),
        "reference_artifacts": {
            "z5_path": str(z5_path),
            "z5_sha256": _sha(z5_path),
            "z8_path": str(z8_path),
            "z8_sha256": _sha(z8_path),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    rows = [
        _compare_point(args.reference_dir, r_name, nt_name)
        for r_name in ("rlo", "rmid", "rhi")
        for nt_name in ("ntlo", "ntmid", "nthi")
    ]
    payload = {
        "schema_version": 1,
        "experiment": "gated phase-tail correction on 3x3 underflow rectangle",
        "git_head": "a5c7afdce50074ac16a1555c0bc8a204ddae6688",
        "fast_sgwb_sha256": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "freq_adaptive_sha256": _sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
        "reference_source_sha256": {
            name: _sha(ROOT / "stiffgwpy_fast" / name)
            for name in (
                "reference.py",
                "stiff_SGWB.py",
                "LCDM_stiff_Neff.py",
                "exact_background.py",
                "global_param.py",
            )
        },
        "reference_dir": str(args.reference_dir.resolve()),
        "threads": 2,
        "reference_tail_conventions": [5.0, 8.0],
        "points": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
