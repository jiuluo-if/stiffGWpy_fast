"""Check the exact matter-era path when the solved spectrum spans reheating."""

# The script is run by path and adds the repository root before importing the
# local checkout rather than whichever package version is installed.
# ruff: noqa: E402, I001

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config
from stiffgwpy_fast.reference import run_reference
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


RESEARCH = Path(__file__).resolve().parent
PARAMS = {
    "cr": 0.0,
    "T_re": 1.0e-3,
    "kappa10": 1.0e-3,
    "n_t": 1.0,
    "r": 0.01,
    "A_s": 6.7e-160,
    "DN_re": 740.0,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    model = LCDM_SG(**PARAMS)
    started = time.perf_counter()
    ok = SGWB_iter_fast(model, config=resolve_config("fast", threads=2), tol=1e-6)
    elapsed = time.perf_counter() - started
    if not ok:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")

    frequencies = np.asarray(model.f, dtype=float)
    fast_spectrum = np.asarray(model.log10OmegaGW, dtype=float)
    reference = run_reference(
        LCDM_SG(**PARAMS),
        dn_eff=float(model.DN_gw[-1]),
        freq_subset=frequencies,
        z_tail=5.0,
        workers=1,
        rtol=1e-11,
    )
    ref_spectrum = np.log10(np.asarray(reference["Ogw"] - reference["Oj"], dtype=float))
    fre = float(model.f_re)
    record = {
        "parameters": PARAMS,
        "fast_status": "success",
        "n_frequency": int(frequencies.size),
        "frequency_reheating_log10_hz": fre,
        "frequency_min_log10_hz": float(frequencies.min()),
        "frequency_max_log10_hz": float(frequencies.max()),
        "reheating_inside_frequency_support": bool(frequencies.min() < fre < frequencies.max()),
        "fast_root": float(model.DN_gw[-1]),
        "reference_map_at_fast_root": float(reference["DN_gw"]),
        "map_relative_difference": abs(float(model.DN_gw[-1]) - float(reference["DN_gw"]))
        / abs(float(reference["DN_gw"])),
        "max_spectrum_error_dex_same_root": float(np.max(np.abs(fast_spectrum - ref_spectrum))),
        "p95_spectrum_error_dex_same_root": float(np.percentile(np.abs(fast_spectrum - ref_spectrum), 95)),
        "fast_seconds": elapsed,
        "frequency": frequencies.tolist(),
        "fast_logspec": fast_spectrum.tolist(),
        "reference_logspec_same_root": ref_spectrum.tolist(),
        "git_head": __import__("subprocess").check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256": {
            "fast_sgwb.py": sha256(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference.py": sha256(ROOT / "stiffgwpy_fast" / "reference.py"),
        },
    }
    (RESEARCH / "exact_matter_transfer_reheating_probe.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in record.items() if k not in {
        "frequency", "fast_logspec", "reference_logspec_same_root", "source_sha256",
    }}, indent=2))


if __name__ == "__main__":
    main()
