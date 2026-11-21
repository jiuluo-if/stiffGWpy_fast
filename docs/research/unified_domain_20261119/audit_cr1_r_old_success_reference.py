"""Audit the one old-Fast success inside the CR1 log-wide candidate box."""
from __future__ import annotations

import hashlib
import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
BASE_SHA = "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45"
SCREEN = HERE / "screen_cr1_r_logwide_lattice.json"
REFERENCE_ARTIFACT = HERE / "validate_cr1_r_logwide_z8.json"
OUTPUT = HERE / "audit_cr1_r_old_success_reference.json"
CASE_ID = "rectangle_07"

OLD_RUNNER = r'''
import json, numpy as np
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config
from stiffgwpy_fast.stiff_SGWB import LCDM_SG
case=json.load(__import__('sys').stdin)
m=LCDM_SG(**case)
result=SGWB_iter_fast(m,config=resolve_config('fast',threads=2),tol=1e-6,physical_guard=True)
if result is not m:
 print(json.dumps({'success':False,'reason':m.fast_failure_reason})); raise SystemExit
print(json.dumps({'success':True,'root':float(m.DN_gw[-1]),'freqs':np.asarray(m.f).tolist(),
 'log_spectrum':np.asarray(m.log10OmegaGW).tolist(),'n_freq':len(m.f)}))
'''


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    screen = json.loads(SCREEN.read_text(encoding="utf-8"))
    row = next(item for item in screen["cases"] if item["input"]["case_id"] == CASE_ID)
    params = row["input"]["params"]
    old_recorded = row["baseline_fast"]
    current_recorded = row["candidate_fast"]
    archive = subprocess.run(["git", "archive", BASE_SHA], cwd=ROOT,
                             check=True, stdout=subprocess.PIPE).stdout
    with tempfile.TemporaryDirectory(prefix="stiffgwpy_cr1_old_success_") as tmp:
        base = Path(tmp)
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as bundle:
            bundle.extractall(base, filter="data")
        env = os.environ.copy()
        env.update({"NUMBA_NUM_THREADS":"2","NUMBA_THREADING_LAYER":"workqueue",
                    "OMP_NUM_THREADS":"1","OPENBLAS_NUM_THREADS":"1",
                    "MKL_NUM_THREADS":"1"})
        old = json.loads(subprocess.run(
            [sys.executable,"-c",OLD_RUNNER],cwd=base,env=env,
            input=json.dumps(params),text=True,capture_output=True,check=True).stdout)
    if not old["success"]:
        raise RuntimeError(f"Pinned old Fast no longer succeeds: {old}")
    if float(old["root"]) != float(old_recorded["root"]):
        raise RuntimeError("Replayed pinned old root differs from paired-screen record")

    ref_artifact = json.loads(REFERENCE_ARTIFACT.read_text(encoding="utf-8"))
    current = next(
        item for item in ref_artifact["cases"]
        if item.get("params", {}).get("r") == params["r"]
        and item.get("params", {}).get("T_re") == params["T_re"]
    )
    current_f = np.asarray(current["raw"]["freqs_log10_hz"], dtype=np.float64)
    old_f = np.asarray(old["freqs"], dtype=np.float64)
    if not np.array_equal(old_f, current_f):
        raise RuntimeError("Old and current success-point frequency grids differ")
    model = __import__("stiffgwpy_fast.stiff_SGWB", fromlist=["LCDM_SG"]).LCDM_SG(**params)
    reference = __import__("stiffgwpy_fast.reference", fromlist=["run_reference"]).run_reference
    started = time.perf_counter()
    ref = reference(model, dn_eff=float(params["DN_eff"])+float(old["root"]),
                    freq_subset=old_f, z_tail=8.0, rtol=1e-9,
                    self_consistent=False, workers=2)
    elapsed = time.perf_counter()-started
    signal = np.asarray(ref["Ogw"]) - np.asarray(ref["Oj"])
    valid = np.isfinite(signal) & (signal > 0.0)
    old_spectrum = np.asarray(old["log_spectrum"], dtype=np.float64)
    direct_log = np.full(signal.shape, np.nan)
    direct_log[valid] = np.log10(signal[valid])
    old_root = float(old["root"])
    new_root = float(current["fast_root"])
    if new_root != float(current_recorded["fast_root"]):
        raise RuntimeError("Current Reference-validation root differs from paired-screen record")
    ref_map_old = float(ref["DN_gw"])
    result = {
        "status":"complete_old_success_precision_audit",
        "reference_code_modified":False,
        "old_failure_baseline_commit":BASE_SHA,
        "reference_code_sha256":sha256(ROOT/"stiffgwpy_fast"/"reference.py"),
        "params":params,
        "old_fast_root":old_root,
        "current_fast_root":new_root,
        "old_to_current_root_relative_delta":abs(old_root-new_root)/abs(old_root),
        "current_root_relative_error_proxy_from_cached_z8_map":current["root_relative_error_proxy"],
        "old_root_reference_map":ref_map_old,
        "old_root_reference_relative_error":abs(old_root-ref_map_old)/abs(ref_map_old),
        "old_root_spectrum_max_abs_dex_from_reference_raw_linear":float(np.max(np.abs(old_spectrum[valid]-direct_log[valid]))),
        "reference_log_floor_bins":int(np.sum(np.asarray(ref["log10OmegaGW"])==-40.0)),
        "reference_seconds":elapsed,
        "frequency_grid_equal":True,
        "old_fast_n_freq":int(old["n_freq"]),
        "current_fast_n_freq":int(len(current_f)),
        "raw":{
            "freqs_log10_hz":old_f.tolist(),
            "old_fast_log10_OmegaGW":old_spectrum.tolist(),
            "reference_raw_Ogw":np.asarray(ref["Ogw"]).tolist(),
            "reference_raw_Oj":np.asarray(ref["Oj"]).tolist(),
            "reference_direct_log10_Ogw_minus_Oj":direct_log.tolist(),
            "valid_positive_reference_mask":valid.tolist(),
        },
    }
    OUTPUT.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k!="raw"},indent=2))


if __name__ == "__main__":
    main()
