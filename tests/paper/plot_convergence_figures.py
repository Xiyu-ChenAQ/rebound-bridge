"""Redraw manuscript convergence and subsystem figures from saved data; never runs simulations."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
import numpy as np

BASE = Path(__file__).resolve().parent
BLUE, RED, ORANGE = "#245c73", "#a03c30", "#bd601e"

def read(path):
    with path.open(encoding="utf-8", newline="") as f:
        return [{k: (v if k in ("method","sweep") else float(v)) for k,v in row.items()}
                for row in csv.DictReader(f)]

def ticks(ax, values):
    ax.xaxis.set_major_locator(FixedLocator(values))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{v:g}" for v in values]))
    ax.xaxis.set_minor_locator(NullLocator())

def refinement(ax, y=-.28):
    ax.annotate("Finer steps", xy=(.03,y), xytext=(.39,y),
                xycoords="axes fraction", textcoords="axes fraction",
                va="center", ha="center", fontsize=7, color=".3",
                arrowprops=dict(arrowstyle="->", color=".4",lw=.75),
                annotation_clip=False)

def save(fig, out, name):
    fig.savefig(out/(name+".pdf"))
    fig.savefig(out/(name+".png"),dpi=220)
    plt.close(fig)

def galilean(rows, out):
    fig,axs=plt.subplots(2,2,figsize=(7.03,4.85))
    fig.subplots_adjust(left=.10,right=.99,top=.89,bottom=.13,wspace=.29,hspace=.65)
    fig.text(.55,.967,r"Refine the outer step; keep $h/h_{\rm sub}=200$",
             ha="center",fontsize=9)
    fig.text(.55,.505,r"Refine the inner step only; keep $h=0.5$ d",
             ha="center",fontsize=9)
    audit={}
    for col,(key,ylabel) in enumerate([("moon_error_km","Position error [km]"),
                                      ("laplace_error_deg","Laplace-angle error [deg]")]):
        for row,sweep in enumerate(["outer","inner"]):
            ax=axs[row,col]
            ax.set_xscale("log",base=2); ax.set_yscale("log")
            ax.set_ylabel(ylabel)
            ax.set_title(f"({chr(97+2*row+col)}) "+("Position" if col==0 else "Orbital phase"),
                         loc="left",pad=5)
            for method,label,color,marker in [
                ("KDK / order2","KDK",RED,"s"),("EOS LF4_2","LF(4,2)",BLUE,"o")]:
                group=sorted([r for r in rows if r["method"]==method and r["sweep"]==sweep],
                             key=lambda r:r["dt_outer_days" if row==0 else "dt_inner_days"])
                x=np.array([r["dt_outer_days" if row==0 else "dt_inner_days"] for r in group])
                y=np.array([r[key] for r in group])
                legend=label
                if row==0:
                    use=x<=.5
                    slope=np.polyfit(np.log(x[use]),np.log(y[use]),1)[0]
                    legend+=rf" ($p={slope:.2f}$)"
                    audit[f"{method}:{key}"]=float(slope)
                ax.plot(x,y,color=color,marker=marker,label=legend,lw=1.35,ms=4)
            if row==0:
                ax.axvspan(.125,.5,color=".3",alpha=.055,zorder=-5,lw=0)
                ax.set_xlim(.112,1.12)
                ticks(ax,[.125,.25,.5,1])
                ax.set_xlabel(r"Outer step $h$ [d]")
                ax.legend(loc="upper left",fontsize=7)
                ax.text(.97,.035,"Shaded: slope-fit interval",transform=ax.transAxes,
                        ha="right",fontsize=6.5,color=".35")
            else:
                ax.set_xlim(.00125,.125)
                ticks(ax,[.0015625,.00625,.025,.1])
                # Compact numeric labels avoid binary exponents and crowded tick text.
                ax.xaxis.set_major_formatter(FixedFormatter(["0.00156","0.00625","0.025","0.1"]))
                ax.set_xlabel(r"Inner step $h_{\rm sub}$ [d]")
                ax.legend(loc="upper left",fontsize=7)
                finest_kdk=min((r for r in rows if r["method"]=="KDK / order2"
                                and r["sweep"]=="inner"),key=lambda r:r["dt_inner_days"])
                ax.annotate("Outer-coupling floor",
                            xy=(finest_kdk["dt_inner_days"],finest_kdk[key]),
                            xytext=(.42,.60),textcoords="axes fraction",fontsize=7,color=".25",
                            arrowprops=dict(arrowstyle="->",color=".4",lw=.7))
                refinement(ax,y=-.32)
    save(fig,out,"galilean_two_timestep_convergence")
    return audit

def crossover(rows,out):
    fig,(top,bottom)=plt.subplots(2,1,figsize=(3.48,5.05),sharex=True,
                                  gridspec_kw={"height_ratios":[1.25,1]})
    fig.subplots_adjust(left=.18,right=.98,top=.935,bottom=.15,hspace=.32)
    settings=[(.01,BLUE,"o","Weak"),(.1,"#9babb3","v",""),
              (1.,"#b8afa0","D",""),(5.,ORANGE,"s","Strong")]
    audit={}
    for scale,color,marker,name in settings:
        g=sorted([r for r in rows if r["method"]=="lf42" and r["coupling_scale"]==scale],
                 key=lambda r:r["dt_days"])
        x=np.array([r["dt_days"] for r in g]); y=np.array([r["max_position_error_au"] for r in g])
        threshold=10*g[0]["reference_floor_au"]; keep=y>threshold
        eps=g[0]["effective_epsilon"]
        exponent=int(np.floor(np.log10(eps))); mantissa=eps/10**exponent
        label=(name+": " if name else "")+rf"$\epsilon={mantissa:.1f}\times10^{{{exponent}}}$"
        top.plot(x,y,color=color,lw=1.4 if name else .85,
                 linestyle="-" if name else "--",label=label,zorder=3 if name else 2)
        top.plot(x[keep],y[keep],ls="none",marker=marker,color=color,ms=4,zorder=4)
        top.plot(x[~keep],y[~keep],ls="none",marker=marker,mfc="white",mec=color,ms=4,zorder=5)
        if name:
            mid=np.sqrt(x[:-1]*x[1:])
            p=np.diff(np.log(y))/np.diff(np.log(x)); ok=keep[:-1]&keep[1:]
            # Horizontal spans identify the two actual steps used to compute each slope.
            bottom.errorbar(mid[ok],p[ok],
                            xerr=np.vstack((mid[ok]-x[:-1][ok],x[1:][ok]-mid[ok])),
                            color=color,marker=marker,ms=4,lw=1.2,elinewidth=.55,
                            capsize=1.5,label=name,zorder=4)
            audit[name]={"steps":x.tolist(),"errors":y.tolist(),"threshold":threshold,
                         "slope_steps":mid[ok].tolist(),"slopes":p[ok].tolist(),
                         "excluded_steps":x[~keep].tolist()}
    top.set_xscale("log",base=2);top.set_yscale("log")
    top.set_xlim(.21,9.5);top.set_ylim(3e-13,1.8e-2)
    top.set_ylabel("Position error [au]")
    top.set_title("(a) LF(4,2): vary coupling strength",loc="left",pad=6)
    # Arrange strong/weak in the same column, retain intermediate values in the other.
    h,l=top.get_legend_handles_labels()
    top.legend([h[i] for i in [0,3,1,2]],[l[i] for i in [0,3,1,2]],
               loc="upper left",ncol=2,fontsize=6.1,handlelength=1.1,
               columnspacing=.65,handletextpad=.35)
    # Reference slopes show scaling only; they are not fits to the data.
    x=np.array([.3,.65]);y=3e-7*(x/.3)**2
    top.plot(x,y,ls=":",color=".35",lw=.9);top.text(.40,6e-7,r"$h^2$",fontsize=8,color=".25")
    x=np.array([1.1,2.3]);y=1.3e-11*(x/1.1)**4
    top.plot(x,y,ls=":",color=".35",lw=.9);top.text(1.4,1.05e-11,r"$h^4$",fontsize=8,color=".25")
    bottom.set_ylim(1.1,5.15)
    bottom.set_yticks([2,3,4,5]);bottom.set_ylabel("Slope between adjacent steps")
    bottom.set_title("(b) Weak and strong coupling",loc="left",pad=6)
    for order in [2,4]:
        bottom.axhline(order,color=".5",ls="--",lw=.7,zorder=1)
        bottom.text(.235,order+.07,f"Order {order}",fontsize=6.5,color=".35",
                    bbox=dict(facecolor="white",edgecolor="none",pad=.2))
    bottom.legend(loc="center right",fontsize=7)
    ticks(bottom,[.25,.5,1,2,4,8])
    bottom.set_xlabel(r"Outer step $h$ [d]")
    refinement(bottom,y=-.32)
    save(fig,out,"generalized_order_crossover_column")
    return audit


def multiple(data_dir,out):
    summary=json.loads((data_dir/"summary.json").read_text(encoding="utf-8"))
    spatial=read(data_dir/"spatial_scaling.csv")
    fig,axes=plt.subplots(1,3,figsize=(7.03,2.95),gridspec_kw={"width_ratios":[1,1,1.03]})
    fig.subplots_adjust(left=.08,right=.992,bottom=.25,top=.82,wspace=.48)
    def curve(ax,coupling,method,key,color,label,marker="o",ls="-"):
        rows=sorted([r for r in summary["metrics"] if r["coupling"]==coupling
                     and r["method"]==method and r["epsilon"]==1e-12],
                    key=lambda r:r["h_days"])
        ax.loglog([r["h_days"] for r in rows],[r[key] for r in rows],
                  color=color,label=label,marker=marker,ms=2.8,ls=ls,lw=1.05)
    a,b,c=axes
    curve(a,"resolved","lf4_2","error_full_model",BLUE,"Full coupling")
    curve(a,"monopole","lf4_2","error_full_model",ORANGE,"Default coupling","s")
    gap=summary["default_vs_full_reference_model_gap"]
    a.axhline(gap,color=ORANGE,ls="--",lw=.8)
    a.text(.0045,gap*2.5,"Model difference",fontsize=6.5,color=ORANGE)
    curve(b,"resolved","lf4_2","error_same_model",BLUE,"Full, LF(4,2)")
    curve(b,"monopole","lf4_2","error_same_model",ORANGE,"Default, LF(4,2)","s","--")
    curve(b,"resolved","kdk","error_same_model",".5","Full, KDK","^",":")
    b.axhline(max(summary["reference_tolerance_check"].values()),color=".55",ls=":",lw=.7)
    for ax in [a,b]:
        ax.set_xlim(.003,2.5);ax.set_ylim(3e-14,.3)
        ticks(ax,[.00390625,.03125,.25,2])
        ax.xaxis.set_major_formatter(FixedFormatter(["1/256","1/32","1/4","2"]))
        ax.set_xlabel(r"Outer step $h$ [d]")
        ax.set_yticks([1e-12,1e-9,1e-6,1e-3])
        ax.legend(loc="upper left",fontsize=6.1,handlelength=1.4)
        refinement(ax,y=-.34)
    a.set_ylabel("Normalized position error")
    b.tick_params(axis="y",labelsize=6.5)
    a.set_title("(a) Compare with full gravity\nLF(4,2) in both modes",fontsize=7.5,loc="left",pad=7)
    b.set_title("(b) Compare with each model\nIntegration error",fontsize=7.5,loc="left",pad=7)
    d=np.array([r["D_au"] for r in spatial]);u=np.abs([r["delta_U"] for r in spatial])
    c.loglog(d,u,color=BLUE,lw=1,marker="o",ms=2.8,label=r"$|U_{\rm full}-U_{\rm app}|$")
    c.loglog(d,u[-1]*(d/d[-1])**-5,color=".4",ls="--",lw=.9,label=r"$D^{-5}$ guide")
    c.axvspan(d[-5],d[-1],color=".4",alpha=.06,lw=0)
    c.set_title("(c) Vary separation only\nOmitted interaction",fontsize=7.5,loc="left",pad=7)
    c.set_xlabel(r"Separation $D$ [au]")
    c.set_ylabel(r"$|\Delta U|$ [$M_\odot\,{\rm au}^2\,{\rm yr}^{-2}$]",labelpad=2)
    c.legend(loc="upper right",fontsize=6.1,handlelength=1.4)
    ticks(c,[.4,1,4,10])
    c.text(.06,.09,f"Fitted slope: {summary['delta_U_log_slope_last5']:.3f}",
           transform=c.transAxes,fontsize=6.5)
    c.text(.06,.025,"Shaded: fit interval",transform=c.transAxes,fontsize=6,color=".4")
    save(fig,out,"multiple_subsystems")
    return {"model_gap":gap,"potential_slope":summary["delta_U_log_slope_last5"],
            "source":str(data_dir),"spatial_points":len(spatial)}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="+",
                        choices=["galilean", "generalized_order", "multiple_subsystems"],
                        default=["galilean", "generalized_order", "multiple_subsystems"],
                        help="Render only the selected figures; by default render all three.")
    parser.add_argument("--multiple-data",type=Path,default=BASE/"results/multiple_subsystems_20261005_verified")
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--galilean-summary",type=Path,default=BASE/"results/two_timestep_convergence/galilean_two_timestep_convergence_summary.csv")
    parser.add_argument("--crossover-csv",type=Path,default=BASE/"results/generalized_order/generalized_order_crossover.csv")
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    plt.style.use(BASE/"mnras.mplstyle")
    result={}
    if "galilean" in args.only:
        result["galilean_slopes"]=galilean(read(args.galilean_summary),args.output_dir)
    if "generalized_order" in args.only:
        result["crossover"]=crossover(read(args.crossover_csv),args.output_dir)
    if "multiple_subsystems" in args.only:
        result["multiple"]=multiple(args.multiple_data,args.output_dir)
    selected=[name for name in ["galilean", "generalized_order", "multiple_subsystems"]
              if name in args.only]
    audit_name="plot_audit.json" if len(selected)==3 else "plot_audit_"+"_".join(selected)+".json"
    (args.output_dir/audit_name).write_text(json.dumps(result,indent=2),encoding="utf-8")
    if "galilean_slopes" in result:
        print(json.dumps(result["galilean_slopes"],indent=2))
    for k,v in result.get("crossover",{}).items(): print(k, "slopes:",v["slopes"],"excluded:",v["excluded_steps"])
if __name__=="__main__":main()
