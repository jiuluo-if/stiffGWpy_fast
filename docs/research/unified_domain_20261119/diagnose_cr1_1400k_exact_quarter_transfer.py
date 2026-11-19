"""Offline exact-quarter transfer prototype against saved CR1 mode traces."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.exact_background import fast_phi_s2_split  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HERE = Path(__file__).resolve().parent
TRAJECTORY = HERE / "cr1_1400k_dominant_mode_trajectories.json"
OUTPUT = HERE / "cr1_1400k_exact_quarter_transfer.json"
PARAMS = {
    "Omega_bh2": 7.0e5, "Omega_ch2": 7.0e5, "H0": 67.32117,
    "DN_eff": 0.0, "A_s": 2.100549e-9, "r": 0.01, "n_t": 0.0,
    "cr": 1.0, "T_re": 2000.0, "DN_re": 10.0, "kappa10": 0.01,
}


def _propagate_exact_quarters(Nv, Phi, Phi_mid, Phi_fine, S2inv, j0, z0,
                              freq, kink_index, kink_fraction, phi_re):
    xh, yh = 0.0, math.exp(z0) * S2inv[j0]
    log_gain = 0.0
    P0 = Phi[j0] - Nv[j0]
    samples = {int(j0): (z0, 0.0, math.exp((Phi[j0]-Nv[j0]-P0)+log_gain)*yh)}
    z = z0
    for k in range(int(j0), len(Nv)-1):
        step = float(Nv[k+1]-Nv[k])
        if z >= 5.0:
            break
        z_mid = float(z0 + Phi_mid[k]-Phi[j0])
        z_end_exact = float(z0 + Phi[k+1]-Phi[j0])
        if k == kink_index and 0.0 < kink_fraction < 1.0:
            z_break = float(z0 + phi_re-Phi[j0])
            xh,yh=FS._phase_segment(xh,yh,z,z_break,step*kink_fraction,0.25)
            xh,yh=FS._phase_segment(xh,yh,z_break,z_end_exact,step*(1.0-kink_fraction),0.25)
        elif abs(z_end_exact-(2.0*z_mid-z)) > 1.0e-5:
            q1=float(z0+Phi_fine[4*k+1]-Phi[j0])
            mid=float(z0+Phi_fine[4*k+2]-Phi[j0])
            q3=float(z0+Phi_fine[4*k+3]-Phi[j0])
            bounds=(z,q1,mid,q3,z_end_exact)
            for j in range(4):
                xh,yh=FS._phase_segment(xh,yh,bounds[j],bounds[j+1],0.25*step,0.25)
        else:
            xh,yh=FS._phase_segment(xh,yh,z,2.0*z_mid-z,step,0.25)
        size=max(abs(xh),abs(yh))
        if size>2.0**450 or (0.0<size<2.0**-450):
            xh/=size; yh/=size; log_gain+=math.log(size)
        k+=1
        z=float(z0+Phi[k]-Phi[j0])
        if k%8==0 or k==len(Nv)-1:
            scale=math.exp((Phi[k]-Nv[k]-P0)+log_gain)
            samples[k]=(z,scale*xh,scale*yh)
        if z>=5.0:
            break
    ids=np.asarray(sorted(samples),dtype=np.int64)
    return ids,np.asarray([samples[i][0] for i in ids]),np.asarray([samples[i][1] for i in ids]),np.asarray([samples[i][2] for i in ids])


def main() -> None:
    model=LCDM_SG(**PARAMS)
    if FS.SGWB_iter_fast(model,config=FS.resolve_config("fast",threads=2),tol=1e-6,physical_guard=True) is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    dn=float(model.DN_gw[-1]); Nv=np.asarray(model.Nv,dtype=np.float64)
    trace=json.loads(TRAJECTORY.read_text(encoding="utf-8"))
    freq=np.asarray([r["frequency"] for r in trace["modes"] if r["variant"]=="base"],dtype=np.float64)
    _,_,j0s,z0s,_=FS.prep_frequency_only(model,Nv,freq)
    Phi,Pm,_,S2inv,ki,kf,phi_re=fast_phi_s2_split(model,Nv,dn,sigma_nodes=model.sigma,relative_scale=True)
    fine_nodes=np.sort(np.concatenate((Nv,Nv[:-1]+0.25*np.diff(Nv),Nv[:-1]+0.5*np.diff(Nv),Nv[:-1]+0.75*np.diff(Nv))))
    Phi_fine=fast_phi_s2_split(model,fine_nodes,dn,relative_scale=True)[0]
    base_rows=[r for r in trace["modes"] if r["variant"]=="base"]
    results=[]
    for i,(fr,j,zstart,base) in enumerate(zip(freq,j0s,z0s,base_rows)):
        ids,zf,xf,yf=_propagate_exact_quarters(Nv,Phi,Pm,Phi_fine,S2inv,int(j),float(zstart),float(fr),int(ki),float(kf),float(phi_re))
        # Match the saved Reference arrays at the common coarse nodes.
        full_ids=np.asarray(base["sample_N"],float)
        refx=np.asarray(base["reference_x"],float); refy=np.asarray(base["reference_y"],float)
        base_x=np.asarray(base["fast_x"],float); base_y=np.asarray(base["fast_y"],float)
        if not np.array_equal(Nv[ids],full_ids):
            keep=np.isin(Nv[ids],full_ids); ids=ids[keep];zf=zf[keep];xf=xf[keep];yf=yf[keep]
        if not np.array_equal(Nv[ids],full_ids):
            raise RuntimeError("candidate/Reference shared-node sets differ")
        scale=np.maximum(np.maximum(np.abs(refx),np.abs(refy)),1e-250)
        e=np.maximum(np.abs(xf-refx),np.abs(yf-refy))/scale
        eb=np.maximum(np.abs(base_x-refx),np.abs(base_y-refy))/scale
        results.append({'frequency':float(fr),'max_state_error_base':float(np.max(eb)),'max_state_error_exact_quarter':float(np.max(e)),'p95_state_error_base':float(np.quantile(eb,.95)),'p95_state_error_exact_quarter':float(np.quantile(e,.95)),'shared_nodes':int(len(ids))})
    payload={'params':PARAMS,'DN_eff':dn,'curvature_gate_abs_delta_z':1e-5,'panel_transfer':'four quarter segments using exact primitive z values','modes':results,'source_sha256':{n:hashlib.sha256((ROOT/'stiffgwpy_fast'/n).read_bytes()).hexdigest() for n in ('fast_sgwb.py','exact_background.py','reference.py')}}
    OUTPUT.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8');print(json.dumps(payload,indent=2))


if __name__ == "__main__":
    main()
