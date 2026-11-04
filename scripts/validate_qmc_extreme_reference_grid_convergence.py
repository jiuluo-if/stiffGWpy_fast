"""Check full-grid Reference convergence for QMC Sobol point 146."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from scripts.pressure_fast_domain import sobol_cases  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.freq_adaptive import grid_independent_freqs  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

SEED = 20261105
POWER = 8
SAMPLE_INDEX = 146
DEFAULT_FREQ_RESOLUTIONS = (2.0, 4.0)
Z_TAIL = 5.0
RTOL = 1e-9
OUTPUT = ROOT / "docs" / "qmc_sobol146_reference_grid_convergence_20261105.json"
RES1_PATH = ROOT / "docs" / "qmc_sobol146_reference_res1_prufer_20261105.json"


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _frequency_sha(freqs: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest()


def _source_hashes() -> dict[str, str]:
    return {
        "fast_module": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "reference_module": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "frequency_module": _sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
        "pressure_script": _sha(ROOT / "scripts" / "pressure_fast_domain.py"),
        "script": _sha(Path(__file__)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--freq-res", type=float, action="append",
        help="Run only these frequency resolutions and merge prior recorded maps.",
    )
    args = parser.parse_args()
    resolutions = tuple(args.freq_res or DEFAULT_FREQ_RESOLUTIONS)
    item = sobol_cases(POWER, SEED)[SAMPLE_INDEX]
    params = item["params"]
    res1 = json.loads(RES1_PATH.read_text(encoding="utf-8"))
    dn_eff = float(res1["DN_eff_fixed"])
    model = LCDM_SG(**params)
    model.cosmo_param["DN_eff"] = dn_eff
    res1_freqs = np.asarray(grid_independent_freqs(model, 1.0)[0], dtype=float)
    res1_logomega = np.asarray(res1["reference_log10OmegaGW"], dtype=float)
    if res1_freqs.size != res1_logomega.size:
        raise ValueError("res=1 frequency coordinates do not match archived spectrum")
    bg = R._Background(model, dn_eff)
    warm_started = time.perf_counter()
    R._warmup_numba_reference_background(bg)
    warm_s = time.perf_counter() - warm_started

    maps: dict[str, dict[str, object]] = {}
    previous_warnings: list[object] = []
    if OUTPUT.exists():
        prior = json.loads(OUTPUT.read_text(encoding="utf-8"))
        maps.update(prior.get("maps", {}))
        previous_warnings = list(
            prior.get("pre_scale_revision_warnings", prior.get("warnings", []))
        )
    if "res1" not in maps or "script_sha256" not in maps["res1"]:
        maps["res1"] = {
            "frequency_res": 1.0,
            "frequency_count": int(res1["reference_frequency_count"]),
            "frequency_sha256": res1["reference_frequency_sha256"],
            "log10freqs": res1_freqs.tolist(),
            "DN_gw": float(res1["Reference_res1_DN_gw"]),
            "log10OmegaGW": res1_logomega.tolist(),
            "runtime_s": float(res1["reference_runtime_s"]),
            "quadrature_error": float(res1["reference_quadrature_error"]),
            "interpolation_error": float(res1["reference_interpolation_error"]),
            "runtime_warnings": [],
            "source_sha256": {
                "fast_module": res1["fast_module_sha256"],
                "reference_module": res1["reference_module_sha256"],
                "frequency_module": _sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
                "pressure_script": res1["pressure_script_sha256"],
                "script": res1["script_sha256"],
            },
        }
    current_warnings: list[dict[str, str | int]] = []

    for freq_res in resolutions:
        started = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = R.run_reference(
                model,
                dn_eff=dn_eff,
                freq_res=freq_res,
                z_tail=Z_TAIL,
                rtol=RTOL,
                self_consistent=False,
                workers=1,
                physical_guard=False,
                backend="prufer",
            )
        elapsed = time.perf_counter() - started
        run_warnings = [
            {
                "message": str(w.message),
                "filename": w.filename,
                "lineno": w.lineno,
                "category": w.category.__name__,
            }
            for w in caught
        ]
        current_warnings.extend(run_warnings)
        freqs = np.asarray(result["freqs"], dtype=float)
        omega = np.asarray(result["Ogw"], dtype=float) - np.asarray(
            result["Oj"], dtype=float
        )
        if not np.all(np.isfinite(omega) & (omega > 0.0)):
            raise FloatingPointError(
                f"Reference res={freq_res:g} has nonpositive/nonfinite spectrum"
            )
        if not np.all(np.diff(freqs) < 0.0):
            raise ValueError(f"Reference res={freq_res:g} frequency order changed")
        maps[f"res{int(freq_res)}"] = {
            "frequency_res": freq_res,
            "frequency_count": int(freqs.size),
            "frequency_sha256": _frequency_sha(freqs),
            "frequency_range": [float(freqs[-1]), float(freqs[0])],
            "log10freqs": freqs.tolist(),
            "DN_gw": float(result["DN_gw"]),
            "log10OmegaGW": np.log10(omega).tolist(),
            "runtime_s": elapsed,
            "quadrature_error": float(result["quadrature_error"]),
            "interpolation_error": float(result["interpolation_error"]),
            "n_iter": result["n_iter"],
            "runtime_warnings": run_warnings,
            "source_sha256": _source_hashes(),
        }

    base_dn = float(maps["res1"]["DN_gw"])
    for key, entry in maps.items():
        entry["relative_DN_gw_delta_from_res1"] = abs(
            float(entry["DN_gw"]) - base_dn
        ) / max(abs(float(entry["DN_gw"])), abs(base_dn), np.finfo(float).tiny)
        fine_f = np.asarray(entry["log10freqs"], dtype=float)
        fine_spectrum = np.asarray(entry["log10OmegaGW"], dtype=float)
        interp = np.interp(res1_freqs[::-1], fine_f[::-1], fine_spectrum[::-1])
        entry["res1_nodes_max_abs_dex_from_this_grid"] = float(
            np.max(np.abs(res1_logomega[::-1] - interp))
        )

    record = {
        "schema_version": 1,
        "experiment": "sobol_20261105_0146_reference_grid_convergence",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "qmc_seed": SEED,
        "qmc_power": POWER,
        "sample_index": SAMPLE_INDEX,
        "label": item["label"],
        "parameters": params,
        "DN_eff_fixed": dn_eff,
        "reference_backend": "prufer",
        "reference_rtol": RTOL,
        "reference_z_tail": Z_TAIL,
        "reference_workers": 1,
        "warmup_s": warm_s,
        "maps": maps,
        "warnings": current_warnings,
        "pre_scale_revision_warnings": previous_warnings,
        "source_sha256": {
            **_source_hashes(),
            "res1_artifact": _sha(RES1_PATH),
        },
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "limitations": [
            "Fixed-background maps at one extreme QMC point; this is not a self-consistent root or a domain certification.",
            "Frequency-grid convergence is assessed via res=1/2/4/8 integrated maps; ODE tolerance and tail depth have separate archived checks.",
            "The coarse-node spectral metric compares frequency-grid interpolation, not two ODE solutions evaluated at identical frequencies.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    summary = {
        key: {
            "count": value["frequency_count"],
            "DN_gw": value["DN_gw"],
            "relative_DN_gw_delta_from_res1": value["relative_DN_gw_delta_from_res1"],
            "runtime_s": value["runtime_s"],
        }
        for key, value in maps.items()
    }
    print(json.dumps({"maps": summary, "warmup_s": warm_s, "warnings": current_warnings, "output": str(OUTPUT)}, indent=2))


if __name__ == "__main__":
    main()
