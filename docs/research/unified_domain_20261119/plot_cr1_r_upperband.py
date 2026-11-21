"""Plot sampled CR1 upper-r coverage and Reference precision boundaries."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def load_cases(lattice_name: str, sobol_name: str) -> list[dict]:
    left=json.loads((HERE/lattice_name).read_text(encoding="utf-8"))["cases"]
    right=json.loads((HERE/sobol_name).read_text(encoding="utf-8"))["cases"]
    return left+right


mid=load_cases("screen_cr1_r_midband_lattice.json","screen_cr1_r_midband_sobol.json")
upper=load_cases("screen_cr1_r_upperband_candidate_lattice.json","screen_cr1_r_upperband_candidate_sobol.json")
old_repeat=json.loads((HERE/"audit_cr1_r_upperband_old_repeatability_20.json").read_text(encoding="utf-8"))
mid_refs=json.loads((HERE/"validate_cr1_r_midband_z8.json").read_text(encoding="utf-8"))["cases"]
mid_refs.append(json.loads((HERE/"validate_cr1_r_mid_upper_probe_z8.json").read_text(encoding="utf-8"))["case"])
upper_refs=json.loads((HERE/"validate_cr1_r_upperband_candidate_z8.json").read_text(encoding="utf-8"))["cases"]
mid_closure=json.loads((HERE/"close_reference_cr1_r_mid_upper_z8.json").read_text(encoding="utf-8"))["case"]
upper_closure=json.loads((HERE/"close_reference_cr1_r_upperband_sobol_z8.json").read_text(encoding="utf-8"))["case"]
upper_low_closure=json.loads((HERE/"close_reference_cr1_r_upperband_low_corner_z8.json").read_text(encoding="utf-8"))["case"]
rejected=json.loads((HERE/"validate_cr1_r_upperband_upper_probe_z8.json").read_text(encoding="utf-8"))["case"]

fig,axes=plt.subplots(1,3,figsize=(15.4,4.6),constrained_layout=True)
ax=axes[0]
for rows in (mid,upper):
    for row in rows:
        p=row["input"]["params"]
        old=row["baseline_fast"]
        if old.get("failure_reason")=="resource_timeout":
            ax.scatter(p["T_re"],p["r"],marker="^",c="#9467bd",s=42,zorder=3)
        elif old.get("success"):
            ax.scatter(p["T_re"],p["r"],marker="s",c="#4c78a8",s=40,zorder=3)
        else:
            ax.scatter(p["T_re"],p["r"],marker="o",facecolors="none",
                       edgecolors="#2ca02c",linewidths=1.2,s=42,zorder=2)
# The repeated old-source audit confirms this point succeeds despite one screen-run failure.
old_p=old_repeat["params"]
ax.scatter(old_p["T_re"],old_p["r"],marker="s",facecolors="none",
           edgecolors="#1f4e79",linewidths=1.8,s=74,zorder=4)
ax.scatter(rejected["params"]["T_re"],rejected["params"]["r"],marker="x",
           c="#d62728",linewidths=2,s=62,zorder=5)
ax.add_patch(plt.Rectangle((1750,.005),250,.0025,fill=False,edgecolor="#2ca02c",
                           linewidth=1.3,linestyle="--"))
ax.add_patch(plt.Rectangle((1875,.005),125,.005,fill=False,edgecolor="#4c78a8",
                           linewidth=1.3,linestyle="-"))
ax.set(xlabel=r"$T_{re}$",ylabel=r"$r$",title="Paired Fast samples and accepted slices")
ax.set_xlim(1738,2010)
ax.set_ylim(.0045,.0107)
ax.grid(alpha=.25)
ax.plot([],[],"o",mfc="none",mec="#2ca02c",label="old numerical fail → current success")
ax.plot([],[],"s",mfc="none",mec="#1f4e79",label="old success, repeated 20/20")
ax.plot([],[],"^",c="#9467bd",label="old resource timeout")
ax.plot([],[],"x",c="#d62728",label="Reference rejected point")
ax.legend(fontsize=7,loc="best")

ax=axes[1]
labels=["mid-band\nold nonfinite","upper-band\nold nonfinite","upper-band\ntimeout","current Fast\nsuccess"]
values=[18,18,1,19]
colors=["#d62728","#d62728","#9467bd","#2ca02c"]
ax.bar(labels,values,color=colors)
ax.set(ylabel="sample count",title="Screen classifications")
ax.tick_params(axis="x",labelrotation=15)
ax.grid(axis="y",alpha=.25)
ax.text(2,1.3,"120 s recheck\nstill timed out",ha="center",fontsize=7)
ax.text(1,.3,"Old holdout repeat audit: 20/20 success",ha="center",fontsize=7)

ax=axes[2]
for rows,color,label in ((mid_refs,"#2ca02c","mid-band"),(upper_refs,"#4c78a8","upper-band")):
    for row in rows:
        root=row["root_relative_error_proxy"]
        spec=row.get("spectrum_max_abs_dex_from_raw_linear",row.get("spectrum_max_abs_dex"))
        # Use the converged third-map residual for the tightly closed holdouts.
        if row.get("case_id")==mid_closure["case_id"]:
            root=mid_closure["fast_root_relative_to_third_map"]
        if row.get("case_id")==upper_closure["case_id"]:
            root=upper_closure["fast_root_relative_to_third_map"]
        if row.get("case_id")==upper_low_closure["case_id"]:
            root=upper_low_closure["fast_root_relative_to_third_map"]
        ax.scatter(root,spec,c=color,s=34,alpha=.8)
ax.scatter(rejected["root_relative_error_proxy"],rejected["spectrum_max_abs_dex_from_raw_linear"],
           c="#d62728",marker="x",s=55,label="T=1750, r=.01 (rejected)")
ax.axvline(2e-4,color="#4c78a8",linestyle="--",linewidth=1,label="root gate")
ax.axhline(1e-3,color="#f2a541",linestyle=":",linewidth=1,label="spectrum gate")
ax.set(xlabel="root relative error",ylabel="spectrum max error (dex)",title="z8 Reference checks")
ax.grid(alpha=.25)
ax.legend(fontsize=7,loc="best")
ax.set_xlim(0, max(.00025,rejected["root_relative_error_proxy"]*1.1))
ax.set_ylim(0,.0014)

fig.suptitle("CR1 upper-r adjacent regions; finite-sample evidence")
out=HERE/"cr1_r_upperband_evidence.png"
fig.savefig(out,dpi=180)
print(out)
