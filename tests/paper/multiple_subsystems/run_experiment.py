"""Run and audit the small two-subsystem benchmark (no timing claims).

python tests/paper/multiple_subsystems/run_experiment.py --exe build/paper/Release/multiple_subsystems.exe --output tests/paper/results/multiple_subsystems
"""
from pathlib import Path
import argparse,csv,hashlib,json,math,platform,subprocess,sys,time
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
G=4*math.pi**2
A=np.array([.02,.025])

def git_metadata(root):
    """Keep source archives runnable even without a Git checkout or command."""
    unavailable={'git_head':None,'git_status':None}
    if not (root/'.git').exists():return unavailable
    try:
        return {'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True,stderr=subprocess.PIPE).strip(),
                'git_status':subprocess.check_output(['git','status','--short'],cwd=root,text=True,stderr=subprocess.PIPE)}
    except (OSError,subprocess.CalledProcessError):
        return unavailable


def read_state(path):
    a=np.genfromtxt(path,delimiter=',',names=True)
    assert len(a)==325 and np.all(np.isfinite(np.column_stack([a[k] for k in a.dtype.names])))
    assert np.array_equal(a['body'].reshape(-1,5),np.tile(np.arange(5),(65,1)))
    x=np.stack([a[k] for k in ['x','y','z']],axis=-1).reshape(-1,5,3)
    v=np.stack([a[k] for k in ['vx','vy','vz']],axis=-1).reshape(-1,5,3)
    m=a['m'][:5]
    assert np.allclose(a['t_days'].reshape(-1,5)[:,0],np.arange(65)*2,rtol=0,atol=1e-13)
    return x,v,m

def com(x,m,indices):
    return np.sum(x[...,indices,:]*m[indices,None],axis=-2)/m[indices].sum()

def error(x,ref):
    rel=x[:,[2,4]]-x[:,[1,3]]
    rref=ref[:,[2,4]]-ref[:,[1,3]]
    return float(np.max(np.linalg.norm(rel-rref,axis=-1)/A))

def conservation(x,v,m,coupling):
    P=np.sum(m[None,:,None]*v,axis=1)
    L=np.sum(m[None,:,None]*np.cross(x,v),axis=1)
    E=.5*np.sum(m[None,:]*np.sum(v*v,axis=-1),axis=1)
    for i in range(5):
        for j in range(i+1,5):E-=G*m[i]*m[j]/np.linalg.norm(x[:,i]-x[:,j],axis=-1)
    if coupling=='monopole':
        ca,cb=com(x,m,[1,2]),com(x,m,[3,4]);ma,mb=m[1:3].sum(),m[3:5].sum()
        Ufull=np.zeros(len(x));Uapp=G*ma*mb/np.linalg.norm(ca-cb,axis=-1)
        for i in [1,2]:
            Uapp-=G*m[i]*mb/np.linalg.norm(x[:,i]-cb,axis=-1)
            for j in [3,4]:Ufull-=G*m[i]*m[j]/np.linalg.norm(x[:,i]-x[:,j],axis=-1)
        for j in [3,4]:Uapp-=G*ma*m[j]/np.linalg.norm(ca-x[:,j],axis=-1)
        E+=Uapp-Ufull
    return {'max_relative_energy_error':float(np.max(np.abs(E-E[0]))/abs(E[0])),
            'max_scaled_momentum_error':float(np.max(np.linalg.norm(P-P[0],axis=-1))/np.sum(m*np.linalg.norm(v[0],axis=-1))),
            'max_relative_angular_momentum_error':float(np.max(np.linalg.norm(L-L[0],axis=-1))/np.linalg.norm(L[0]))}

def spatial_checks(initial,kick_file,out):
    mp.mp.dps=70
    x,v,m=initial
    masses=list(map(lambda a:mp.mpf(str(a)),m))
    pos=[[mp.mpf(str(a)) for a in row] for row in x[0]]
    gm=mp.mpf(str(G))
    norm=lambda a:mp.sqrt(sum(z*z for z in a))
    sub=lambda a,b:[q-r for q,r in zip(a,b)]
    def centres(p):
        return [[sum(masses[i]*p[i][k] for i in ids)/sum(masses[i] for i in ids) for k in range(3)] for ids in [[1,2],[3,4]]]
    ma,mb=sum(masses[1:3]),sum(masses[3:5])
    def potentials(p):
        ca,cb=centres(p)
        full=-gm*sum(masses[i]*masses[j]/norm(sub(p[i],p[j])) for i in [1,2] for j in [3,4])
        u00=-gm*ma*mb/norm(sub(ca,cb))
        us0=-gm*mb*sum(masses[i]/norm(sub(p[i],cb)) for i in [1,2])
        u0t=-gm*ma*sum(masses[j]/norm(sub(ca,p[j])) for j in [3,4])
        return full,us0+u0t-u00,u00
    def B(p,resolved):
        ca,cb=centres(p);full,app,mono=potentials(p)
        star=-gm*masses[0]*sum(masses[i]/norm(sub(p[i],p[0])) for i in range(1,5))
        star_mono=-gm*masses[0]*(ma/norm(sub(ca,p[0]))+mb/norm(sub(cb,p[0])))
        return star-star_mono+(full if resolved else app)-mono
    # Derivatives of the written Hamiltonians; independent of the C force formulas.
    kicks=list(csv.DictReader(kick_file.open(encoding='utf-8')))
    gradient_errors=[]
    for row in kicks:
        i,k=int(row['body']),int(row['component'])
        def function(z):
            q=[r.copy() for r in pos];q[i][k]=z
            return B(q,row['coupling']=='resolved')
        expected=-mp.diff(function,pos[i][k])/masses[i]
        gradient_errors.append(abs(float(expected)-float(row['actual_B'])))
    assert max(gradient_errors)<1e-10,max(gradient_errors)
    ca,cb=centres(pos)
    offsets=[[mp.mpf(0)]*3]+[sub(pos[i],ca if i<3 else cb) for i in range(1,5)]
    rs2=sum(masses[i]*norm(offsets[i])**2 for i in [1,2])/ma
    rt2=sum(masses[i]*norm(offsets[i])**2 for i in [3,4])/mb
    rows=[]
    for exponent in range(10):
        D=mp.mpf('.4')*mp.mpf('1.5')**exponent
        p=[r.copy() for r in offsets]
        for i in [3,4]:p[i][0]+=D
        full,app,mono=potentials(p)
        delta=full-app;scale=gm*ma*mb*rs2*rt2/D**5
        rows.append({'D_au':float(D),'delta_U':float(delta),'relative_to_monopole':float(abs(delta/mono)),
                     'leading_scale':float(scale),'coefficient':float(delta/scale)})
    with (out/'spatial_scaling.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    # Fixed a priori fit window: five most separated snapshots.
    slope=float(np.polyfit(np.log([r['D_au'] for r in rows[-5:]]),np.log([abs(r['delta_U']) for r in rows[-5:]]),1)[0])
    return rows,{'potential_gradient_kick_max_abs_error':max(gradient_errors),'delta_U_log_slope_last5':slope,
                 'r_s_au':float(mp.sqrt(rs2)),'r_t_au':float(mp.sqrt(rt2)), 'arithmetic_decimal_digits':mp.mp.dps}

def analyse(out,manifest,data_dir=None):
    data_dir=out if data_dir is None else data_dir
    refs={}
    for c in ['resolved','monopole']:
        for e in ['1e-13','1e-14']:refs[c,e]=read_state(data_dir/f'reference_{c}_{e}.csv')
    base=refs['resolved','1e-14']
    ref_uncertainty={c:error(refs[c,'1e-13'][0],refs[c,'1e-14'][0]) for c in ['resolved','monopole']}
    floor=error(refs['monopole','1e-14'][0],base[0])
    component_checks=[]
    for run in manifest['runs']:
        if run['kind']=='bridge' and run['epsilon']==1e-14:
            loose=data_dir/f"bridge_{run['coupling']}_{run['method']}_{run['h_days']:g}.csv"
            component_checks.append({'coupling':run['coupling'],'h_days':run['h_days'],
                                     'trajectory_change':error(read_state(loose)[0],read_state(data_dir/run['file'])[0])})
    metrics=[]
    for run in manifest['runs']:
        if run['kind']!='bridge':continue
        x,v,m=read_state(data_dir/run['file']);same=refs[run['coupling'],'1e-14'][0]
        assert np.max(abs(x[0]-base[0][0]))<1e-14,'initial mismatch'
        metrics.append({'coupling':run['coupling'],'method':run['method'],'epsilon':run['epsilon'],
                        'h_days':run['h_days'],'error_same_model':error(x,same),'error_full_model':error(x,base[0]),
                        **conservation(x,v,m,run['coupling'])})
    with (out/'convergence.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=metrics[0].keys());w.writeheader();w.writerows(metrics)
    spatial,checks=spatial_checks(base,data_dir/'kick_check.csv',out)
    slopes=[]
    for c in ['resolved','monopole']:
        for method in ['kdk','lf4_2']:
            series=sorted([r for r in metrics if r['coupling']==c and r['method']==method and r['epsilon']==1e-12],key=lambda r:r['h_days'],reverse=True)
            for a,b in zip(series,series[1:]):
                slopes.append({'coupling':c,'method':method,'h_coarse_days':a['h_days'],'h_fine_days':b['h_days'],
                               'local_order_same_model':math.log(a['error_same_model']/b['error_same_model'])/math.log(a['h_days']/b['h_days'])})
    with (out/'local_orders.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=slopes[0].keys());w.writeheader();w.writerows(slopes)
    # Geometric diagnostics establish that the fixed hierarchy remains separated.
    x,v,m=base;ca,cb=com(x,m,[1,2]),com(x,m,[3,4]);D=np.linalg.norm(ca-cb,axis=-1)
    diagnostics={'min_subsystem_separation_au':float(D.min()),'max_subsystem_separation_au':float(D.max()),
                 'internal_separation_au_min':np.linalg.norm(x[:,[2,4]]-x[:,[1,3]],axis=-1).min(axis=0).tolist(),
                 'internal_separation_au_max':np.linalg.norm(x[:,[2,4]]-x[:,[1,3]],axis=-1).max(axis=0).tolist()}
    summary={'reference_tolerance_check':ref_uncertainty,'component_tolerance_check':component_checks,
             'default_vs_full_reference_model_gap':floor,
             **checks,**diagnostics,'metrics':metrics,'local_orders':slopes,
             'definitions':{'error':'max over 65 epochs and two subsystems of |(r_moon-r_host)-(r_moon-r_host)_ref| / initial a',
                            'conservation':'sampled maxima; energy evaluated with the selected Hamiltonian',
                            'reference':'independent direct IAS15; default model uses differentiated proxy potential'},
             'limitations':['Synthetic fixed hierarchy, not ephemeris reconstruction.',
                            'IAS15 component flows isolate outer splitting; this does not test WHFast cache or establish exact symplecticity.',
                            'Finite-size sweep changes only subsystem separation at frozen internal geometry.',
                            'No timing or general scalability claim.']}
    (out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    plot(out,summary,spatial)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['metrics','local_orders']},indent=2),flush=True)

def plot(out,summary,spatial):
    plt.rcParams.update({'font.family':'DejaVu Serif','font.size':8,'axes.labelsize':8,'legend.fontsize':6.5,
                         'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(7.2,2.9),layout='constrained')
    specs=[('resolved','kdk','error_full_model','Full coupling, KDK','#547896','o','-'),
           ('resolved','lf4_2','error_full_model','Full coupling, LF(4,2)','#202a34','s','-'),
           ('monopole','lf4_2','error_full_model','Default coupling vs full IAS15','#b34f44','^','-'),
           ('monopole','lf4_2','error_same_model','Default coupling vs its own IAS15','#b34f44','^',':')]
    for c,m,key,label,color,marker,ls in specs:
        rows=sorted([r for r in summary['metrics'] if r['coupling']==c and r['method']==m and r['epsilon']==1e-12],key=lambda r:r['h_days'])
        ax[0].loglog([r['h_days'] for r in rows],[r[key] for r in rows],marker=marker,ms=3,color=color,ls=ls,label=label)
    ax[0].axhline(summary['default_vs_full_reference_model_gap'],color='#b34f44',ls='--',lw=.8,label='Difference between the two IAS15 models')
    ax[0].axhline(max(summary['reference_tolerance_check'].values()),color='.5',ls='-.',lw=.7,label='IAS15 tolerance check')
    ax[0].set(xlabel='Outer timestep (day)',ylabel='Maximum normalized internal-position error',title='(a) Time error and model error')
    ax[0].legend(loc='best',frameon=False);ax[0].grid(alpha=.15,which='major')
    d=np.array([r['D_au'] for r in spatial]);u=np.array([abs(r['delta_U']) for r in spatial])
    ax[1].loglog(d,u,'o-',color='#547896',ms=3,label=r'$|U^{\rm full}-U^{\rm app}|$ (70-digit arithmetic)')
    ax[1].loglog(d,u[-1]*(d/d[-1])**-5,'--',color='.3',label=r'$D^{-5}$ guide')
    ax[1].set(xlabel='Frozen subsystem separation D (au)',ylabel=r'Omitted potential $|\Delta U|$ ($M_\odot\,\mathrm{au}^2\,\mathrm{yr}^{-2}$)',title='(b) Spatial approximation')
    ax[1].legend(frameon=False);ax[1].grid(alpha=.15,which='major')
    ax[1].text(.06,.08,f"Last-five-point slope: {summary['delta_U_log_slope_last5']:.3f}",transform=ax[1].transAxes,fontsize=8)
    for suffix in ['pdf','png']:fig.savefig(out/f'multiple_subsystems.{suffix}',dpi=240,bbox_inches='tight')
    plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',type=Path,help='Benchmark executable; required unless --analyse-only is used.')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--data-dir',type=Path,help='Saved raw runs for --analyse-only; defaults to --output.')
    p.add_argument('--pilot',action='store_true');p.add_argument('--analyse-only',action='store_true');a=p.parse_args()
    if not a.analyse_only and a.exe is None:p.error('--exe is required unless --analyse-only is used')
    if not a.analyse_only and a.data_dir is not None:p.error('--data-dir requires --analyse-only')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    if a.analyse_only:
        data_dir=a.data_dir.resolve() if a.data_dir is not None else out
        return analyse(out,json.loads((data_dir/'manifest.json').read_text(encoding='utf-8')),data_dir)
    exe=a.exe.resolve()
    assert exe.is_file(),exe
    assert not (out/'manifest.json').exists(),'Use a fresh output folder; existing runs are preserved.'
    hvalues=[.5,.0625,.0078125] if a.pilot else [2.**(1-i) for i in range(10)]
    source_files=[Path(__file__),Path(__file__).with_name('multiple_subsystems.c')]+list((ROOT/'src').glob('*.c'))+[ROOT/'src/include/rebound_bridge.h']
    manifest={'experiment':'two_subsystem_fivebody','duration_days':128,'sample_days':2,'h_values_days':hvalues,
              'initial_conditions':{'star_mass':1,'host_masses':[.001,.0008],'moon_masses':[.00005,.00004],
                'outer_COM_radii_au':[1,1.45],'outer_phase':[0,0],'outer_velocity':'astrocentric circular using star + own subsystem mass',
                'inner_a_au':[.02,.025],'inner_e':[.05,.04],'inner_inc_rad':[.1,-.12],'inner_true_anomaly_rad':[.3,1.8],
                'inner_Omega_rad':[0,0],'inner_omega_rad':[0,0],'frame':'whole-system centre of mass'},
              'component_integrator':'IAS15','component_epsilon':1e-12,'python':sys.version,'platform':platform.platform(),
              **git_metadata(ROOT),
              'source_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in source_files},
              'exe_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'runs':[]}
    def execute(command):
        start=time.perf_counter();r=subprocess.run(command,capture_output=True,text=True,timeout=120)
        with (out/'execution.log').open('a',encoding='utf-8') as log:log.write(str(command)+'\n'+r.stdout+r.stderr+'\n')
        r.check_returncode();return time.perf_counter()-start
    execute([str(exe),'--kick-check',str(out/'kick_check.csv')])
    for c in ['resolved','monopole']:
        for eps in ['1e-13','1e-14']:
            name=f'reference_{c}_{eps}.csv'
            elapsed=execute([str(exe),'reference',c,'kdk','.5',eps,'128','2',str(out/name)])
            manifest['runs'].append({'kind':'reference','coupling':c,'epsilon':float(eps),'file':name,'elapsed_seconds':elapsed})
    for c in ['resolved','monopole']:
        for method in ['kdk','lf4_2']:
            for h in hvalues:
                name=f'bridge_{c}_{method}_{h:g}.csv'
                elapsed=execute([str(exe),'bridge',c,method,str(h),'1e-12','128','2',str(out/name)])
                manifest['runs'].append({'kind':'bridge','coupling':c,'method':method,'h_days':h,'epsilon':1e-12,'file':name,'elapsed_seconds':elapsed})
                (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
                print(c,method,h,round(elapsed,3),flush=True)
    # Tighten numerical component errors at one intermediate and the finest step.
    if not a.pilot:
        for c in ['resolved','monopole']:
            for h in [.125,hvalues[-1]]:
                name=f'bridge_{c}_lf4_2_{h:g}_tight.csv'
                elapsed=execute([str(exe),'bridge',c,'lf4_2',str(h),'1e-14','128','2',str(out/name)])
                manifest['runs'].append({'kind':'bridge','coupling':c,'method':'lf4_2','h_days':h,'epsilon':1e-14,'file':name,'elapsed_seconds':elapsed})
    manifest['output_sha256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in out.glob('*.csv')}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    analyse(out,manifest)

if __name__=='__main__':main()
