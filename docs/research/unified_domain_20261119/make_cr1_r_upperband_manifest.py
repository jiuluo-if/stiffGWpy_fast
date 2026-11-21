"""Create source/evidence hashes for the CR1 high-r adjacent extensions."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUTPUT=HERE/"cr1_r_upperband_manifest.json"
SOURCE_FILES=["exact_background.py","fast_sgwb.py","freq_adaptive.py","_wkb_transfer.py","LCDM_stiff_Neff.py","reference.py"]
ARTIFACTS=[
    "cr1_r_upperband_summary.md","cr1_r_upperband_evidence.png","plot_cr1_r_upperband.py",
    "make_cr1_r_upperband_manifest.py",
    "screen_cr1_r_upperband.py","screen_cr1_r_upperband_lattice.json","screen_cr1_r_upperband_sobol.json",
    "validate_cr1_r_upperband_upper_probe_z8.py","validate_cr1_r_upperband_upper_probe_z8.json",
    "validate_cr1_r_upperband_edge_probes_z8.py","validate_cr1_r_upperband_edge_probes_z8.json",
    "screen_cr1_r_upperband_candidate.py","screen_cr1_r_upperband_candidate_lattice.json","screen_cr1_r_upperband_candidate_sobol.json",
    "validate_cr1_r_upperband_candidate_z8.py","validate_cr1_r_upperband_candidate_z8.json",
    "close_reference_cr1_r_upperband_z8.py","close_reference_cr1_r_upperband_z8.json",
    "close_reference_cr1_r_upperband_sobol_z8.py","close_reference_cr1_r_upperband_sobol_z8.json",
    "close_reference_cr1_r_upperband_low_corner_z8.py","close_reference_cr1_r_upperband_low_corner_z8.json",
    "audit_cr1_r_upperband_old_success_z8.py","audit_cr1_r_upperband_old_success_z8.json",
    "audit_cr1_r_upperband_old_repeatability.py","audit_cr1_r_upperband_old_repeatability_20.json",
    "audit_cr1_r_upperband_timeout_recheck.py","audit_cr1_r_upperband_timeout_recheck.json",
    "screen_cr1_r_midband.py","screen_cr1_r_midband_lattice.json","screen_cr1_r_midband_sobol.json",
    "validate_cr1_r_midband_z8.py","validate_cr1_r_midband_z8.json",
    "validate_cr1_r_mid_upper_probe_z8.py","validate_cr1_r_mid_upper_probe_z8.json",
    "close_reference_cr1_r_mid_upper_z8.py","close_reference_cr1_r_mid_upper_z8.json",
    "audit_cr1_r_midband_old_success_z8.json",
    "benchmark_cr1_r_upperband_ab.py","benchmark_cr1_r_upperband_ab_500pairs.json",
    "benchmark_cr1_r_midband_ab.py","benchmark_cr1_r_midband_ab_500pairs.json",
    "benchmark_cr1_r_midband_low_corner_1000pairs.py","benchmark_cr1_r_midband_low_corner_ab_1000pairs.json",
]


def sha256(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main()->None:
    payload={
        "status":"validated_finite_sample_cr1_upper_r_adjacent_regions",
        "accepted_regions":[
            {"T_re":[1750.0,2000.0],"r":[0.005,0.0075],"kappa10":0.01},
            {"T_re":[1875.0,2000.0],"r":[0.005,0.01],"kappa10":0.01},
        ],
        "fixed_parameters":{"cr":1.0,"Omega_bh2":700000.0,"Omega_ch2":700000.0,
                            "H0":67.32117,"DN_eff":0.0,"A_s":2.100549e-9,"DN_re":10.0},
        "derived_tilt":"n_t=-r/8; red throughout accepted boxes",
        "old_failure_baseline_commit":"872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "performance_baseline_commit":"d4ee1b9c8f9fad401c5e22c2316094a608074727",
        "reference_code_modified":False,"production_code_changed_this_tranche":False,
        "continuum_claim":False,
        "paired_screen":{"midband":{"points":19,"old_success":1,"current_success":19,
                                      "old_nonfinite_background_primitive":18},
                         "upperband_raw":{"points":19,"old_success":0,"current_success":19,
                                          "old_nonfinite_background_primitive":18,"old_resource_timeout":1},
                         "upperband_old_success_repeat_audit":{"repeats":20,"successes":20,
                                                              "one_paired_screen_run_had_failure":True}},
        "accuracy":{"reference_cases_midband":4,"reference_cases_upperband":6,
                    "max_root_relative_to_third_map":0.00019497517047322135,
                    "max_spectrum_abs_error_dex":0.000902850545459799,
                    "rejected_T1750_r001_root_error":0.0002414678255871807},
        "performance":{"upperband_p50_max":1.0362213247131074,"upperband_p95_max":1.02857910187396,
                        "midband_1500pair_corner_p50_ratio":1.0109237994154043,
                        "midband_1500pair_corner_p95_ratio":1.0261009611821417,
                        "trusted_cr0_upperband_p50":0.9972490935639965,
                        "trusted_cr0_upperband_p95":1.0030097328250753,
                        "trusted_cr0_midband_p50":0.9931174775984756,
                        "trusted_cr0_midband_p95":1.000959835818977},
        "source_sha256":{name:sha256(ROOT/"stiffgwpy_fast"/name) for name in SOURCE_FILES},
        "artifact_sha256":{name:sha256(HERE/name) for name in ARTIFACTS},
        "manifest_self_hash_included":False,
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(OUTPUT)


if __name__=="__main__":
    main()
