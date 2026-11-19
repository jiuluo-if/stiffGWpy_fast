from __future__ import annotations
import hashlib,json,math,os,subprocess,sys,time
from pathlib import Path
import psutil
root=Path(__file__).resolve().parents[3];sys.path.insert(0,str(root));p=psutil.Process();aff=p.cpu_affinity()[:2];p.cpu_affinity(aff)
for k,v in {'NUMBA_NUM_THREADS':'2','NUMBA_THREADING_LAYER':'workqueue','FAST_THREADS':'2','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'}.items():os.environ[k]=v
import numpy as np
from scipy.special import jv,jvp,yv,yvp
from numba import set_num_threads
from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast._wkb_transfer import wkb_transfer_embedded

def exact_bessel(z0,z1,h):
 q=(z1-z0)/h;alpha=q+1+.25*q*q;order=math.sqrt(alpha)/abs(q)
 def basis(z,n):
  w=math.exp(z);arg=w/abs(q);scale=math.exp(.5*q*n)
  jj=float(jv(order,arg));yy=float(yv(order,arg));djj=float(jvp(order,arg,1));dyy=float(yvp(order,arg,1))
  xj=scale*jj;xy=scale*yy;xnj=scale*(.5*q*jj+q*arg*djj);xny=scale*(.5*q*yy+q*arg*dyy)
  return np.array([[xj,xy],[-(xnj+xj)/w,-(xny+xy)/w]])
 return basis(z1,h)@np.linalg.solve(basis(z0,0.),np.eye(2))

def calc_matrix(fn,z0,z1,h):
 c0=fn(1.,0.,z0,z1,h);c1=fn(0.,1.,z0,z1,h)
 return np.array([[c0[0],c1[0]],[c0[1],c1[1]]])

def phase(x,y,z0,z1,h):return FS._phase_segment_base(x,y,z0,z1,h,.25)

data=json.loads((root/'docs/research/unified_domain_20261119/diagnose_cr0_phase_substep_work.json').read_text())
segs=[];seen=set()
for case in data['cases']:
 if case['case'] not in ('cr0_rmin_nt50','cr0_rmax_nt50'):continue
 for row in case['largest_phase_panels']:
  key=(round(row['z_start'],10),round(row['z_end'],10),round(row['h_segment'],12))
  if row['n_sub']>=10000 and key not in seen:segs.append((case['case'],row));seen.add(key)
set_num_threads(2);rows=[]
for case,row in segs[:8]:
 z0=row['z_start'];z1=row['z_end'];h=row['h_segment'];q=(z1-z0)/h;nsub=FS._phase_substeps(h,.5*(z0+z1),.25)
 exact=exact_bessel(z0,z1,h)
 def wfn(x,y,*args):
  a,b,e,ok=wkb_transfer_embedded(x,y,*args,1e-3)
  if not ok:raise RuntimeError(f'WKB rejected {e}')
  return a,b
 wm=calc_matrix(wfn,z0,z1,h);pm=calc_matrix(phase,z0,z1,h)
 wtime=[];ptime=[]
 for _ in range(8):
  t=time.perf_counter();calc_matrix(wfn,z0,z1,h);wtime.append(time.perf_counter()-t)
  t=time.perf_counter();calc_matrix(phase,z0,z1,h);ptime.append(time.perf_counter()-t)
 errw=np.linalg.norm(wm-exact,ord=2)/max(np.linalg.norm(exact,ord=2),1e-300)
 errp=np.linalg.norm(pm-exact,ord=2)/max(np.linalg.norm(exact,ord=2),1e-300)
 poww=np.max(np.abs(np.sum(wm*wm,axis=0)-np.sum(exact*exact,axis=0)))
 defect=max(wkb_transfer_embedded(1.,0.,z0,z1,h,1e-3)[2],wkb_transfer_embedded(0.,1.,z0,z1,h,1e-3)[2])
 rows.append({'case':case,'panel':row['panel'],'mode':row['mode'],'q':q,'n_phase_substeps':nsub,'wkb3_defect_basis':defect,'matrix_rel_error_vs_bessel':errw,'phase_matrix_rel_error_vs_bessel':errp,'max_column_power_absdiff_vs_bessel':float(poww),'wkb3_transfer_p50_s':float(np.median(wtime)),'phase_substeps_p50_s':float(np.median(ptime)),'speed_ratio':float(np.median(ptime)/np.median(wtime))})
for r in rows:print(r)
out={'status':'complete_local_panel_screen','source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'source_sha256':{name:hashlib.sha256((root/'stiffgwpy_fast'/name).read_bytes()).hexdigest() for name in ('fast_sgwb.py','_wkb_transfer.py')},'resources':{'affinity':aff,'threads':2},'oracle':'exact constant-q Bessel fundamental matrix; phase chain is original production _phase_segment_base with phase_max=.25','cases':rows}
path=root/'docs/research/unified_domain_20261119/benchmark_wkb3_high_phase_panels.json';path.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(path)
