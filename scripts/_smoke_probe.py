# -*- coding: utf-8 -*-
import os
import tempfile
from pathlib import Path

import yaml
from cobaya import mpi
from cobaya.run import run

# The smoke test is intentionally single-process.  GitHub runners may have
# mpi4py installed without a system libmpi; explicitly disabling MPI keeps
# Cobaya's serial execution path deterministic and portable.
mpi.set_mpi_disabled()

os.chdir(tempfile.gettempdir())  # different cwd (qualified-name test)
repo = str(Path(__file__).resolve().parents[1])
base = os.path.join(repo, "stiffgwpy", "cobaya")
yaml_text = open(os.path.join(base, "mcmc_compare.yaml"), encoding="utf-8").read()
run_info = yaml.safe_load(yaml_text)
run_info["theory"] = {"stiffgwpy.cobaya.stiffGW.stiffGW": {
    "engine": "fast", "fallback": True, "accuracy_mode": "production",
    "fast_threads": 8}}
for sec in ("theory", "likelihood"):
    for cfg in (run_info.get(sec) or {}).values():
        if not isinstance(cfg, dict):
            continue
        pp = cfg.get("python_path")
        if isinstance(pp, str) and not os.path.isabs(pp):
            cfg["python_path"] = os.path.join(base, pp)
        if sec == "likelihood":
            for key in cfg:
                if key.endswith("_file"):
                    v = cfg[key]
                    if isinstance(v, str) and not os.path.isabs(v):
                        cfg[key] = os.path.join(base, v)
run_info["sampler"]["mcmc"]["max_samples"] = 5
run_info["sampler"]["mcmc"]["seed"] = 20260830
run_info["output"] = os.path.join(tempfile.mkdtemp(), "out")
print("CWD:", os.getcwd())
upd, sampler = run(run_info, force=True, resume=False)
prod = sampler.products()
print("columns:", list(prod["sample"].columns))
print("n samples:", len(prod["sample"]))
print("OK")
