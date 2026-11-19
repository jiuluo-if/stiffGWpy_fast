"""Recheck the rejected feature-band point after exact matter-era transfer."""

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
    "T_re": 10.0,
    "kappa10": 1e-6,
    "n_t": 1.0,
    "r": 0.01,
    "A_s": 1.9680104830424577e-7,
    "DN_re": 740.0,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    baseline = json.loads(
        (RESEARCH / "per_mode_error_feature_band.json").read_text(encoding="utf-8")
    )
    model = LCDM_SG(**PARAMS)
    start = time.perf_counter()
    ok = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2),
        freq_grid="grid_independent",
        freq_res=4,
        tol=1e-6,
    )
    elapsed = time.perf_counter() - start
    if not ok:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")

    fast_freq = np.asarray(model.f, dtype=float)
    fast_spectrum = np.asarray(model.log10OmegaGW, dtype=float)
    ref = run_reference(
        LCDM_SG(**PARAMS),
        dn_eff=float(model.DN_gw[-1]),
        freq_subset=fast_freq,
        z_tail=5.0,
        workers=1,
        rtol=1e-11,
    )
    ref_logspec = np.log10(np.asarray(ref["Ogw"] - ref["Oj"], dtype=float))
    if len(fast_freq) != len(baseline["freq"]):
        raise AssertionError("Fast mode count differs from the archived comparison grid")
    old_freq = np.asarray(baseline["freq"], dtype=float)
    if not np.allclose(fast_freq, old_freq, rtol=0.0, atol=2e-12):
        raise AssertionError("Fast frequency support differs from the archived comparison grid")

    record = {
        "parameters": PARAMS,
        "fast_status": "success",
        "fast_root": float(model.DN_gw[-1]),
        "fast_seconds": elapsed,
        "reference_map_at_fast_root": float(ref["DN_gw"]),
        "root_relative_error_vs_closed_reference": abs(
            float(model.DN_gw[-1]) - 0.1614848239010858
        )
        / 0.1614848239010858,
        "max_spectrum_error_dex_same_root": float(
            np.max(np.abs(fast_spectrum - ref_logspec))
        ),
        "p95_spectrum_error_dex_same_root": float(
            np.percentile(np.abs(fast_spectrum - ref_logspec), 95)
        ),
        "max_old_spectrum_error_dex_same_root": float(
            np.max(np.abs(np.asarray(baseline["fast_logspec"]) - ref_logspec))
        ),
        "frequency": fast_freq.tolist(),
        "fast_logspec": fast_spectrum.tolist(),
        "reference_logspec_same_root": ref_logspec.tolist(),
        "git_head": __import__("subprocess").check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256": {
            "fast_sgwb.py": sha256(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference.py": sha256(ROOT / "stiffgwpy_fast" / "reference.py"),
        },
    }
    output = RESEARCH / "exact_matter_transfer_feature_probe.json"
    output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in record.items() if k not in {
        "frequency", "fast_logspec", "reference_logspec_same_root",
        "source_sha256",
    }}, indent=2))


if __name__ == "__main__":
    main()
