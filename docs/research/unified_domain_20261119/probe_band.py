import json, math, time, sys, warnings
import numpy as np
from stiffgwpy_fast.stiff_SGWB import LCDM_SG
from stiffgwpy_fast import global_param as gp
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config

base=LCDM_SG()
f0=base.derived_param['H_0']/(2*math.pi)
print('ratio', gp.f_piv/f0,flush=True)
for kind, nt, width in [('normal',2,1),('narrow',1,1.01),('narrow',1,1.05),('narrow',1.5,1.05),('narrow',2,1.1),('narrow',1,1.15),('narrow',2,1.3)]:
    p=dict(cr=0., T_re=2e3, kappa10=.01, n_t=nt, r=.01, DN_re=740.)
    if kind=='narrow':p['A_s']=(gp.f_piv/(width*f0))**nt/p['r']
    row={'params':p,'runs':[]}
    for rg,res,h in [('goal',1,.005),('grid_independent',4,.0025),('grid_independent',8,.00125)]:
        m=LCDM_SG(**p);t=time.perf_counter()
        try:
            with warnings.catch_warnings(record=True) as ww:
                ok=SGWB_iter_fast(m,config=resolve_config('fast',threads=2,h=h),freq_grid=rg,tol=1e-6,freq_res=res)
            v=None if ok is None else (m.DN_gw[-1],len(m.f),[float(m.f[-1]),float(m.f[0])],np.isfinite(m.log10OmegaGW).all().item())
            row['runs'].append([v,m.fast_failure_reason,round(time.perf_counter()-t,3),[str(r.message) for r in ww][:3]])
        except Exception as ex: row['runs'].append([repr(ex)])
    print(json.dumps(row),flush=True)

