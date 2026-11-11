"""Re-run Fast against frozen z5/z8 References on a T_re x r x n_t grid."""

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

R_NAMES = ("rlo", "rmid", "rhi")
NT_NAMES = ("ntlo", "ntmid", "nthi")


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_point(tre: int, r_name: str, nt_name: str) -> dict:
    stem = f"fast_momentblend_3d_tre{tre}_{r_name}_{nt_name}"
    refs = {
        tail: ROOT / "docs" / f"{stem}_{tail}_20261108.json"
        for tail in ("z5", "z8")
    }
    data = {tail: json.loads(path.read_text(encoding="utf-8")) for tail, path in refs.items()}
    params = data["z8"]["params"]
    if any(item["params"] != params for item in data.values()):
        raise RuntimeError(f"z5/z8 parameter mismatch: {stem}")
    rows = {tail: item["frequency_rows"] for tail, item in data.items()}
    frequencies = np.asarray([row["log10_hz"] for row in rows["z8"]], dtype=np.float64)
    if not np.array_equal(
        frequencies,
        np.asarray([row["log10_hz"] for row in rows["z5"]], dtype=np.float64),
    ):
        raise RuntimeError(f"z5/z8 frequency mismatch: {stem}")
    digest = hashlib.sha256(frequencies.astype("<f8").tobytes()).hexdigest()

    model = LCDM_SG(**params)
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
    elapsed = time.perf_counter() - started
    if result is not model or not model.SGWB_converge:
        raise RuntimeError(f"Fast failed for {stem}: {model.fast_failure_reason}")
    if not np.array_equal(np.asarray(model.f, dtype=np.float64), frequencies):
        raise RuntimeError(f"Fast/Reference frequency-grid mismatch: {stem}")
    if hashlib.sha256(np.asarray(model.f, dtype="<f8").tobytes()).hexdigest() != digest:
        raise RuntimeError(f"frequency hash mismatch: {stem}")

    root = float(model.DN_gw[-1])
    spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
    errors = {
        tail: np.abs(
            spectrum
            - np.asarray([row["reference_log10_OmegaGW"] for row in rows[tail]])
        )
        for tail in ("z5", "z8")
    }
    ref_roots = {
        tail: float(data[tail]["scaled_reference"]["estimated_self_consistent_root_DN_gw"])
        for tail in ("z5", "z8")
    }
    return {
        "point": f"tre{tre}_{r_name}_{nt_name}",
        "params": params,
        "frequency_count": int(frequencies.size),
        "frequency_sha256": digest,
        "runtime_s": elapsed,
        "fast_root": root,
        "reference_root_z5": ref_roots["z5"],
        "reference_root_z8": ref_roots["z8"],
        "root_relative_error_z5": (root - ref_roots["z5"]) / ref_roots["z5"],
        "root_relative_error_z8": (root - ref_roots["z8"]) / ref_roots["z8"],
        "max_spectrum_error_dex_z5": float(errors["z5"].max()),
        "max_spectrum_error_dex_z8": float(errors["z8"].max()),
        "p95_spectrum_error_dex_z5": float(np.percentile(errors["z5"], 95)),
        "p95_spectrum_error_dex_z8": float(np.percentile(errors["z8"], 95)),
        "reference_files": {
            tail: {"path": str(path), "sha256": _sha(path)} for tail, path in refs.items()
        },
        "reference_root_method_limit": "single Reference map iteration at Fast root; not a converged Reference root certificate",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tre", type=int, choices=(10, 10000), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    points = [
        _run_point(args.tre, r_name, nt_name)
        for r_name in R_NAMES
        for nt_name in NT_NAMES
    ]
    payload = {
        "schema_version": 1,
        "experiment": "momentblend Fast verification on frozen 3D underflow reference surfaces",
        "tre_GeV": args.tre,
        "fast_source_sha256": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "frequency_source_sha256": _sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
        "threads": 2,
        "points": points,
        "scope_limit": "Reference roots use one map update; finite grid samples do not certify every point in the continuous box.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
