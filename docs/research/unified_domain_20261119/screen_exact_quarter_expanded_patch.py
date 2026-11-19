"""Paired CR1 Fast screen on the 1.30M--1.40M density patch."""

from __future__ import annotations

import hashlib
import io
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

BASE_SHA = "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45"
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "screen_exact_quarter_expanded_patch.json"
BASE_PARAMS = {
    "H0": 67.32117, "DN_eff": 0.0, "A_s": 2.100549e-9,
    "r": 0.01, "n_t": 0.0, "cr": 1.0, "T_re": 2000.0,
    "DN_re": 10.0, "kappa10": 0.01,
}
# Rectangle corners, opposite-edge points and a withheld interior point.
COORDS = ((650_000.0, 650_000.0), (650_000.0, 700_000.0),
          (700_000.0, 650_000.0), (700_000.0, 700_000.0),
          (675_000.0, 675_000.0))

RUNNER = r'''
import contextlib, json, os, sys, time, warnings
sys.path.insert(0, os.getcwd())
import numpy as np
from stiffgwpy_fast.exact_background import H2_vec
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config
from stiffgwpy_fast.stiff_SGWB import LCDM_SG
rows=[]
for item in json.load(sys.stdin):
    params=item["params"]
    model=LCDM_SG(**params)
    start=time.perf_counter()
    try:
        with contextlib.redirect_stdout(sys.stderr), warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result=SGWB_iter_fast(model, config=resolve_config("fast",threads=2),
                                  tol=1e-6, physical_guard=True)
        success=result is model
        failure=model.fast_failure_reason
        root=float(model.DN_gw[-1]) if success else None
        finite=bool(np.isfinite(model.log10OmegaGW).all()) if success else False
        freqs=np.asarray(model.f,dtype=np.float64) if success else np.empty(0)
        spectrum=np.asarray(model.log10OmegaGW,dtype=np.float64) if success else np.empty(0)
        grid_sha=__import__("hashlib").sha256(freqs.tobytes()).hexdigest() if success else None
        spectrum_sha=__import__("hashlib").sha256(spectrum.tobytes()).hexdigest() if success else None
        warning_count=len(caught)
        exception=None
    except Exception as exc:
        success=False; failure=None; root=None; finite=False
        grid_sha=None; spectrum_sha=None; warning_count=0
        exception=type(exc).__name__+": "+str(exc)
    try:
        h2=H2_vec(model.Nv,model,params["DN_eff"])
        h2_ok=bool(np.isfinite(h2).all() and np.all(h2>0.0))
    except Exception:
        h2_ok=False
    rows.append({"case_id":item["case_id"],"params":params,
        "success":success,"failure_reason":failure,"exception":exception,
        "root":root,"finite_spectrum":finite,"frequency_sha256":grid_sha,
        "spectrum_sha256":spectrum_sha,"warning_count":warning_count,
        "h2_positive_finite":h2_ok,"runtime_seconds":time.perf_counter()-start})
print(json.dumps(rows))
'''


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    cases = []
    for index, (ob, oc) in enumerate(COORDS):
        params = dict(BASE_PARAMS, Omega_bh2=ob, Omega_ch2=oc)
        cases.append({"case_id": f"patch_{index:02d}", "params": params})
    archive = subprocess.run(["git", "archive", BASE_SHA], cwd=ROOT,
                             check=True, stdout=subprocess.PIPE).stdout
    with tempfile.TemporaryDirectory(prefix="stiffgwpy_fast_cr1_patch_base_") as tmp:
        base_root = Path(tmp)
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as bundle:
            bundle.extractall(base_root, filter="data")
        baseline = []
        for case in cases:
            result = subprocess.run([sys.executable, "-c", RUNNER], cwd=base_root,
                input=json.dumps([case]), text=True, check=True, capture_output=True,
                timeout=15.0)
            baseline.append(json.loads(result.stdout)[0])
    candidate = []
    for case in cases:
        result = subprocess.run([sys.executable, "-c", RUNNER], cwd=ROOT,
            input=json.dumps([case]), text=True, check=True, capture_output=True,
            timeout=15.0)
        candidate.append(json.loads(result.stdout)[0])
    paired = [{"input": case, "baseline_fast": old, "candidate_fast": new}
              for case, old, new in zip(cases, baseline, candidate)]
    payload = {
        "status": "complete", "base_sha": BASE_SHA,
        "candidate_head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
            check=True, capture_output=True, text=True).stdout.strip(),
        "candidate_source_sha256": {
            name: _sha(ROOT / "stiffgwpy_fast" / name)
            for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
        },
        "rectangle": {"Omega_bh2": [650_000.0, 700_000.0],
                      "Omega_ch2": [650_000.0, 700_000.0]},
        "paired_rows": paired,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"],
        "old_successes": sum(row["baseline_fast"]["success"] for row in paired),
        "candidate_successes": sum(row["candidate_fast"]["success"] for row in paired),
        "old_failures": sorted({row["baseline_fast"]["failure_reason"] for row in paired}),
        "candidate_failures": sorted({row["candidate_fast"]["failure_reason"] for row in paired
                                       if not row["candidate_fast"]["success"]}),
        "results": paired}, indent=2))


if __name__ == "__main__":
    main()
