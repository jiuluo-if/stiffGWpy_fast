"""Continue the direct Reference tail-depth check at the 457-node grid."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import screen_subnormal_map_grid_convergence_20261115 as screen  # noqa: E402

from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude, construct_fast_frequency_grid  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

OUTPUT = ROOT / "docs" / "research" / "results" / "subnormal_map_tail_z12_20261115.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    physical = LCDM_SG(**{**screen.BASE, "r": screen.R})
    physical.gen_expansion()
    log_amp = _log_tensor_amplitude(physical.cosmo_param, physical.derived_param)
    grid_model = LCDM_SG(**{**screen.BASE, "r": screen.R})
    grid_model.gen_expansion()
    freqs = np.asarray(
        construct_fast_frequency_grid(grid_model, freq_res=2.0, log_amplitude=log_amp),
        dtype=np.float64,
    )
    reference = LCDM_SG(**{**screen.BASE, "r": 1.0, "A_s": 1.0})
    result = screen.log_reference_map(reference, freqs, log_amp, 0.0, 12.0)
    with (ROOT / "docs" / "research" / "results" / "subnormal_map_grid_convergence_20261115.json").open(
        encoding="utf-8"
    ) as handle:
        prior = json.load(handle)
    z10 = prior["cases"][1]["reference"]["z10"]["log_DN_gw"]
    output = {
        "schema_version": 1,
        "status": "complete_direct_reference_z12_tail_continuation",
        "scope": "Direct Cartesian Reference at cr=0,n_t=0,r=1e-315 using a unit transfer scaled in log space on the same 457-node freq_res=2 grid.",
        "parameters": {**screen.BASE, "r": screen.R},
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
        "reference_z12": result,
        "z10_to_z12_relative_map_change": math.expm1(result["log_DN_gw"] - z10),
        "source_sha256": {
            "fast_sgwb.py": sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive.py": sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "reference.py": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "screen_subnormal_map_grid_convergence.py": sha(Path(__file__).with_name("screen_subnormal_map_grid_convergence_20261115.py")),
            "script": sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "frequency_count": freqs.size,
        "z12_map": result["DN_gw_float64"],
        "z10_to_z12": output["z10_to_z12_relative_map_change"],
        "runtime_s": result["runtime_s"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
