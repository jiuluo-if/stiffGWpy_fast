"""Compare old/current Fast to Reference at the retained high-r success point."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
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
OLD_SHA = "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45"
DEFAULT_SCREEN = "screen_cr1_r_upperband_candidate_lattice.json"
DEFAULT_OUTPUT = "audit_cr1_r_upperband_old_success_z8.json"
DEFAULT_CASE_ID = "holdout_00"
spec = importlib.util.spec_from_file_location("old_upperband_screen_helper", HERE / "screen_h2_identity_region.py")
if spec is None or spec.loader is None:
    raise RuntimeError("could not load old Fast baseline helper")
screen_helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen_helper)

OLD_RUNNER = r'''
import json,numpy as np,sys
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast,resolve_config
from stiffgwpy_fast.stiff_SGWB import LCDM_SG
p=json.load(sys.stdin);m=LCDM_SG(**p)
out=SGWB_iter_fast(m,config=resolve_config('fast',threads=2),tol=1e-6,physical_guard=True)
if out is not m: print(json.dumps({'success':False,'reason':m.fast_failure_reason})); raise SystemExit
print(json.dumps({'success':True,'root':float(m.DN_gw[-1]),'freqs':np.asarray(m.f).tolist(),
 'log_spectrum':np.asarray(m.log10OmegaGW).tolist()}))
'''


def log_metrics(spectrum: np.ndarray, signal: np.ndarray) -> tuple[float, int]:
    valid = np.isfinite(signal) & (signal > 0.0)
    direct = np.log10(signal[valid])
    return float(np.max(np.abs(spectrum[valid]-direct))), int(valid.sum())


def main(screen_name: str = DEFAULT_SCREEN, case_id: str = DEFAULT_CASE_ID,
         output_name: str = DEFAULT_OUTPUT) -> None:
    screen_path = HERE / screen_name
    output_path = HERE / output_name
    artifact = json.loads(screen_path.read_text(encoding="utf-8"))
    case = next(row for row in artifact["cases"] if row["input"]["case_id"] == case_id)
    params = case["input"]["params"]
    old_source = subprocess.run(["git","archive",OLD_SHA],cwd=ROOT,check=True,stdout=subprocess.PIPE).stdout
    with tempfile.TemporaryDirectory(prefix="stiffgwpy_upperband_old_root_") as temp:
        old_root_dir = Path(temp)
        with tarfile.open(fileobj=io.BytesIO(old_source),mode="r:") as bundle:
            bundle.extractall(old_root_dir,filter="data")
        env=os.environ.copy()
        env.update({"NUMBA_NUM_THREADS":"2","NUMBA_THREADING_LAYER":"workqueue",
                    "OMP_NUM_THREADS":"1","OPENBLAS_NUM_THREADS":"1","MKL_NUM_THREADS":"1"})
        old = json.loads(subprocess.run([sys.executable,"-c",OLD_RUNNER],cwd=old_root_dir,env=env,
                                        input=json.dumps(params),text=True,capture_output=True,check=True).stdout)
    if not old.get("success"):
        raise RuntimeError(f"Pinned old Fast did not succeed in retained-point audit: {old}")

    fast = __import__("stiffgwpy_fast.fast_sgwb",fromlist=["resolve_config"])
    reference = __import__("stiffgwpy_fast.reference",fromlist=["run_reference"]).run_reference
    LCDM_SG = __import__("stiffgwpy_fast.stiff_SGWB",fromlist=["LCDM_SG"]).LCDM_SG
    current_model=LCDM_SG(**params)
    if fast.SGWB_iter_fast(current_model,config=fast.resolve_config("fast",threads=2),
                            tol=1e-6,physical_guard=True) is not current_model:
        raise RuntimeError(f"Current Fast failed: {current_model.fast_failure_reason}")
    old_freq=np.asarray(old["freqs"],dtype=np.float64)
    new_freq=np.asarray(current_model.f,dtype=np.float64)
    old_log=np.asarray(old["log_spectrum"],dtype=np.float64)
    new_log=np.asarray(current_model.log10OmegaGW,dtype=np.float64)
    results={}
    for label,root,spectrum,frequencies in (("old",float(old["root"]),old_log,old_freq),
                                            ("current",float(current_model.DN_gw[-1]),new_log,new_freq)):
        model=LCDM_SG(**params)
        started=time.perf_counter()
        ref=reference(model,dn_eff=float(params["DN_eff"])+root,
                      freq_subset=frequencies,z_tail=8.0,rtol=1e-9,
                      self_consistent=False,workers=2)
        elapsed=time.perf_counter()-started
        signal=np.asarray(ref["Ogw"])-np.asarray(ref["Oj"])
        spec_error,valid_count=log_metrics(spectrum,signal)
        ref_map=float(ref["DN_gw"])
        results[label]={
            "root":root,"reference_map_at_that_root":ref_map,
            "root_relative_error_proxy":abs(root-ref_map)/abs(ref_map),
            "spectrum_max_abs_dex_from_raw_linear":spec_error,
            "spectrum_valid_count":valid_count,"reference_seconds":elapsed,
            "frequency_count":int(len(frequencies)),
            "reference_log_floor_bins":int(np.sum(np.asarray(ref["log10OmegaGW"])==-40.0)),
            "raw":{"reference_Ogw":np.asarray(ref["Ogw"]).tolist(),
                   "reference_Oj":np.asarray(ref["Oj"]).tolist(),
                   "reference_direct_log10_Ogw_minus_Oj":np.log10(signal).tolist(),
                   "fast_log10_OmegaGW":spectrum.tolist(),
                   "freqs_log10_hz":frequencies.tolist()},
        }
    payload={
        "status":"complete_cr1_upperband_old_success_precision_audit",
        "reference_code_modified":False,"old_failure_baseline_commit":OLD_SHA,
        "params":params,"frequency_grid_equal":bool(np.array_equal(old_freq,new_freq)),
        "old_frequency_count":int(len(old_freq)),"current_frequency_count":int(len(new_freq)),
        "old_to_current_root_relative_delta":abs(float(old["root"])-float(current_model.DN_gw[-1]))/abs(float(old["root"])),
        "old":results["old"],"current":results["current"],
        "reference_sha256":hashlib.sha256((ROOT/"stiffgwpy_fast"/"reference.py").read_bytes()).hexdigest(),
    }
    OUTPUT = output_path
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in payload.items() if k not in ("old","current")},indent=2))
    print("old",{k:v for k,v in results["old"].items() if k!="raw"})
    print("current",{k:v for k,v in results["current"].items() if k!="raw"})


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--screen",default=DEFAULT_SCREEN)
    parser.add_argument("--case-id",default=DEFAULT_CASE_ID)
    parser.add_argument("--output",default=DEFAULT_OUTPUT)
    args=parser.parse_args()
    main(args.screen,args.case_id,args.output)
