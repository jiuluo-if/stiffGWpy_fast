"""Paired Fast screen on a CR1 reheating-temperature x kappa patch."""

from __future__ import annotations

import hashlib
import importlib
import io
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

screen_support = importlib.import_module(
    "docs.research.unified_domain_20261119.screen_exact_quarter_expanded_patch"
)

BASE_SHA = "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45"
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "screen_exact_quarter_reheat_kappa_patch.json"
BASE_PARAMS = {
    "Omega_bh2": 700_000.0, "Omega_ch2": 700_000.0,
    "H0": 67.32117, "DN_eff": 0.0, "A_s": 2.100549e-9,
    "r": 0.01, "n_t": 0.0, "cr": 1.0, "DN_re": 10.0,
}
COORDS = tuple(
    (temperature, kappa)
    for temperature in (2_000.0, 2_250.0, 2_500.0)
    for kappa in (0.01, 0.0125, 0.015)
)
HOLDOUT = (2_375.0, 0.01375)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run(root: Path, case: dict[str, object]) -> dict[str, object]:
    result = subprocess.run(
        [sys.executable, "-c", screen_support.RUNNER],
        cwd=root,
        input=json.dumps([case]),
        text=True,
        check=True,
        capture_output=True,
        timeout=30.0,
    )
    return json.loads(result.stdout)[0]


def main() -> None:
    cases = []
    specs = [(f"transition_{index:02d}", temperature, kappa)
             for index, (temperature, kappa) in enumerate(COORDS)]
    specs.append(("transition_holdout", *HOLDOUT))
    for case_id, temperature, kappa in specs:
        params = dict(BASE_PARAMS, T_re=temperature, kappa10=kappa)
        preflight = LCDM_SG(**params)
        cases.append({
            "case_id": case_id,
            "family": "independent_holdout" if case_id == "transition_holdout" else "T_re_x_kappa10",
            "params": params,
            "derived_nt": float(preflight.derived_param["nt"]),
            "derived_N_inf": float(preflight.derived_param["N_inf"]),
        })
    archive = subprocess.run(
        ["git", "archive", BASE_SHA], cwd=ROOT,
        check=True, stdout=subprocess.PIPE,
    ).stdout
    with tempfile.TemporaryDirectory(prefix="stiffgwpy_cr1_transition_base_") as temp:
        base_root = Path(temp)
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as archive_file:
            archive_file.extractall(base_root, filter="data")
        rows = []
        for case in cases:
            baseline = _run(base_root, case)
            candidate = _run(ROOT, case)
            rows.append({"input": case, "baseline_fast": baseline,
                         "candidate_fast": candidate})
    payload = {
        "status": "complete",
        "base_sha": BASE_SHA,
        "candidate_head": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT,
            check=True, capture_output=True, text=True,
        ).stdout.strip(),
        "candidate_source_sha256": {
            name: _sha(ROOT / "stiffgwpy_fast" / name)
            for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
        },
        "slice": {
            "cr": 1.0,
            "Omega_bh2": 700_000.0,
            "Omega_ch2": 700_000.0,
            "T_re": [2_000.0, 2_500.0],
            "kappa10": [0.01, 0.015],
            "lattice": {"T_re": [2_000.0, 2_250.0, 2_500.0],
                        "kappa10": [0.01, 0.0125, 0.015]},
        },
        "paired_rows": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "old_successes": sum(r["baseline_fast"]["success"] for r in rows),
        "candidate_successes": sum(r["candidate_fast"]["success"] for r in rows),
        "old_failures": sorted({r["baseline_fast"]["failure_reason"] for r in rows}),
        "candidate_failures": sorted({r["candidate_fast"]["failure_reason"] for r in rows
                                       if not r["candidate_fast"]["success"]}),
        "rows": rows,
    }, indent=2))


if __name__ == "__main__":
    main()
