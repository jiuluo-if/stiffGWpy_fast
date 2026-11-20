from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

OUTPUT = Path(__file__).resolve().parent / "diagnose_cr0_density_14m_reference_step_underflow.json"
PARAMS = {
    "Omega_bh2": 14.0e6, "Omega_ch2": 0.1201075, "H0": 67.32117,
    "DN_eff": 0.0, "A_s": 2.100549e-9, "r": 0.01, "n_t": -0.5,
    "cr": 0.0, "T_re": 2000.0, "DN_re": 10.0, "kappa10": 0.01,
}


def main() -> None:
    model = LCDM_SG(**PARAMS)
    result = SGWB_iter_fast(model, config=resolve_config("fast", threads=2), tol=1e-6, physical_guard=True)
    if result is not model:
        raise RuntimeError(f"Fast unexpectedly failed: {model.fast_failure_reason}")
    freqs = model.f.copy()
    dn_eff = float(model.cosmo_param["DN_eff"])
    cases = []
    for backend in ("cartesian", "prufer"):
        row = {"backend": backend, "frequency_count": int(freqs.size), "dn_eff": dn_eff, "fast_root": float(model.DN_gw[-1])}
        try:
            out = run_reference(
                model, dn_eff=dn_eff, freq_subset=freqs, z_tail=8.0, rtol=1.0e-9,
                self_consistent=False, workers=1, backend=backend,
            )
            row.update({"status": "success", "finite_modes": int(len(out["Ogw"]))})
        except Exception as exc:
            row.update({"status": "reference_integration_failure", "exception_type": type(exc).__name__, "message": str(exc)})
        cases.append(row)
    payload = {
        "status": "reference_backends_fail_at_14m_low_tilt_boundary" if all(c["status"] != "success" for c in cases) else "mixed_reference_boundary_outcome",
        "params": PARAMS,
        "fast_success": True,
        "fast_root": float(model.DN_gw[-1]),
        "frequency_count": int(freqs.size),
        "source_sha256": {
            name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
            for name in ("fast_sgwb.py", "reference.py", "exact_background.py")
        },
        "reference_code_modified": False,
        "reference_settings": {"z_tail": 8.0, "rtol": 1.0e-9, "self_consistent": False},
        "backends": cases,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
