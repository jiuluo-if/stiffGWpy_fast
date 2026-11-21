"""Extended 1000-pair A/B for the P95-tight CR1 mid-band corner."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=HERE/"benchmark_cr1_upper_vs_fda_head_process_ab.py"
spec=importlib.util.spec_from_file_location("cr1_r_midband_low_corner_ab",SOURCE)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load persistent full-Fast benchmark")
benchmark=importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)
LATTICE=json.loads((HERE/"screen_cr1_r_midband_lattice.json").read_text(encoding="utf-8"))
CASE=next(row["input"]["params"] for row in LATTICE["cases"]
          if row["input"]["case_id"]=="rectangle_02")


if __name__=="__main__":
    benchmark.BASE_SHA="d4ee1b9c8f9fad401c5e22c2316094a608074727"
    benchmark.REPEATS=1000
    benchmark.WARMUPS=3
    benchmark.OUTPUT=HERE/"benchmark_cr1_r_midband_low_corner_ab_1000pairs.json"
    benchmark._cases=lambda:{"low_T_high_r_mid_band_corner":CASE}
    benchmark.main()
