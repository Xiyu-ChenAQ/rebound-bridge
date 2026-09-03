#!/usr/bin/env python3
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def main():
 p=argparse.ArgumentParser(); p.add_argument('--csv',required=True); p.add_argument('--outdir',default='.'); a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
 d={}; times={}
 with open(a.csv) as f:
  for r in csv.DictReader(f): d[(r['method'],r['time_years'],r['body'])]=np.array([float(r['x_au']),float(r['y_au']),float(r['z_au'])]); times[r['method']]=float(r['elapsed_seconds'])
 bodies=sorted({k[2] for k in d}); ts=sorted({k[1] for k in d},key=float); errs={}
 for m in ('WHFast','BridgeLF(4,2)'):
  vals=[]
  for t in ts:
   for b in bodies:
    ref=d.get(('IAS15',t,b)); cur=d.get((m,t,b));
    if ref is not None and cur is not None: vals.append(np.linalg.norm(cur-ref)/max(np.linalg.norm(ref),1e-12))
  errs[m]={'max':float(np.max(vals)),'rms':float(np.sqrt(np.mean(np.square(vals))))}
 plt.figure();
 for m in errs: plt.semilogy([float(t) for t in ts],[max(np.linalg.norm(d[(m,t,b)]-d[('IAS15',t,b)])/max(np.linalg.norm(d[('IAS15',t,b)]),1e-12) for b in bodies) for t in ts],label=m)
 plt.xlabel('Years'); plt.ylabel('Max relative position error'); plt.legend(); plt.tight_layout(); plt.savefig(out/'errors_timeseries.png',dpi=160); plt.close()
 plt.figure(); plt.bar(list(errs),[errs[m]['max'] for m in errs]); plt.yscale('log'); plt.ylabel('Maximum relative error'); plt.tight_layout(); plt.savefig(out/'errors_summary.png',dpi=160); plt.close()
 speed={m:times['IAS15']/v for m,v in times.items()}; plt.figure(); plt.bar(list(speed.keys()),list(speed.values())); plt.ylabel('Speedup vs IAS15'); plt.tight_layout(); plt.savefig(out/'speedup.png',dpi=160); plt.close()
 json.dump({'timings_seconds':times,'speedup':speed,'errors':errs},open(out/'summary.json','w'),indent=2)
if __name__=='__main__': main()
