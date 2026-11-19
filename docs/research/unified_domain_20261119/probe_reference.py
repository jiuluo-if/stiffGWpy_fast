import json, math, time
import numpy as np
from stiffgwpy_fast.stiff_SGWB import LCDM_SG
from stiffgwpy_fast import global_param as gp
from stiffgwpy_fast.reference import run_reference
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast,resolve_config

for nt,width,duration in [(1.,1.01,740.),(1.,1.006,800.),(1.35,1.007,705.)]:
    f0=LCDM_SG().derived_param['H_0']/(2*math.pi)
    p=dict(cr=0.,T_re=2000.,kappa10=.01,n_t=nt,r=.01,DN_re=duration,A_s=(gp.f_piv/(width*f0))**nt/.01)
    m=LCDM_SG(**p);SGWB_iter_fast(m,tol=1e-6,config=resolve_config('fast',threads=2))
    root=float(m.DN_gw[-1]); nf=len(m.f)
    dense=np.linspace(m.f[0],m.f[-1],(nf-1)*2+1)
    rec={'params':p,'fast':root,'logf':m.f.tolist(),'log10OmegaGW':m.log10OmegaGW.tolist(),'maps':[]}
    dn=root
    for _ in range(4):
        t=time.perf_counter()
        ref=run_reference(LCDM_SG(**p),dn_eff=dn,freq_subset=dense,z_tail=5.,workers=1,rtol=1e-11)
        lg=np.log10(ref['Ogw']-ref['Oj'])
        rn=float(ref['DN_gw'])
        rec['maps'].append(dict(dn=dn,map=rn,log10OmegaGW=lg.tolist(),freq=dense.tolist(),raw_max_err=float(np.max(abs(lg[::2]-m.log10OmegaGW))),seconds=time.perf_counter()-t))
        print(nt,width,duration,dn,rn,rec['maps'][-1]['raw_max_err'],flush=True)
        if abs(rn-dn)<abs(dn)*1e-8:break
        dn=rn
    with open(f'probe_reference_{duration:.0f}.json','w') as fh:json.dump(rec,fh)

