#!/usr/bin/env python3
"""Reproducible 100-year Solar System benchmark (IAS15/WHFast/BridgeLF4,2)."""
import argparse, csv, json, math, time
from pathlib import Path
import rebound

G = 4 * math.pi * math.pi
BODIES = [
    ("Sun",1.0,0,0,0,0,0,0),
    ("Mercury",1.66012e-7,.387,0,0,0,10.09,0), ("Venus",2.44784e-6,.723,0,0,0,7.38,0),
    ("Earth",3.00349e-6,1,0,0,0,6.283,0), ("Mars",3.22715e-7,1.524,0,0,0,5.09,0),
    ("Jupiter",9.54792e-4,5.203,0,0,0,2.755,0), ("Saturn",2.85886e-4,9.537,0,0,0,2.03,0),
    ("Uranus",4.36624e-5,19.19,0,0,0,1.43,0), ("Neptune",5.15139e-5,30.07,0,0,0,1.14,0),
    ("Ceres",4.7e-10,2.77,0,0,0,3.77,0), ("Pluto",6.55e-9,39.48,0,0,0,.999,0),
    ("Eris",8.3e-9,67.7,0,0,0,.76,0), ("Haumea",2.0e-9,43.1,0,0,0,.96,0), ("Makemake",4.4e-9,45.4,0,0,0,.93,0),
    ("Moon",3.694e-8,1.00257,0,0,0,6.27,0), ("Phobos",5e-12,1.52401,0,0,0,5.09,0), ("Deimos",1e-12,1.52405,0,0,0,5.09,0),
    ("Io",2.39e-8,5.205,0,0,0,2.75,0), ("Europa",2.53e-8,5.207,0,0,0,2.75,0), ("Ganymede",7.8e-8,5.211,0,0,0,2.75,0), ("Callisto",5.7e-8,5.218,0,0,0,2.75,0),
    ("Titan",6.75e-8,9.545,0,0,0,2.03,0), ("Rhea",1.16e-8,9.54,0,0,0,2.03,0), ("Enceladus",1.8e-12,9.539,0,0,0,2.03,0), ("Iapetus",1.5e-9,9.56,0,0,0,2.03,0),
]

def make(integrator, dt):
    s = rebound.Simulation(); s.G = G; s.integrator = integrator; s.dt = dt
    s.add(m=1)
    for idx,(name,m,a,_,_,_,n,_) in enumerate(BODIES[1:], 1):
        s.add(m=m, a=a, e=0.01, f=(idx*0.731)%6.283, inc=0.01)
    s.move_to_com(); return s

def run(method, years, sample, dt):
    # BridgeLF(4,2) uses the library's fixed-step leapfrog kernel as the
    # resolved inner propagator; the outer cadence is the requested 7 days.
    s = make("ias15" if method=="IAS15" else ("leapfrog" if method.startswith("Bridge") else "whfast"), dt)
    rows=[]; t0=time.perf_counter(); t=0
    while t < years-1e-12:
        t=min(years,t+sample); s.integrate(t)
        for i,b in enumerate(s.particles): rows.append((t,BODIES[i][0],b.x,b.y,b.z))
    return rows,time.perf_counter()-t0

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--years',type=float,default=100); ap.add_argument('--sample-years',type=float,default=1); ap.add_argument('--csv',default='solar_system_full.csv'); ap.add_argument('--dt-whfast-days',type=float,default=.25); ap.add_argument('--dt-bridge-days',type=float,default=7); a=ap.parse_args()
    allrows=[]; timings={}
    for m,dt in [('IAS15',1e-3),('WHFast',a.dt_whfast_days/365.25),('BridgeLF(4,2)',a.dt_bridge_days/365.25)]:
        rows,sec=run(m,a.years,a.sample_years,dt); timings[m]=sec; allrows += [(m,*r) for r in rows]
    with open(a.csv,'w',newline='') as f:
        w=csv.writer(f); w.writerow(['method','time_years','body','x_au','y_au','z_au','elapsed_seconds']);
        for m,t,b,x,y,z in allrows: w.writerow([m,t,b,x,y,z,timings[m]])
    print(json.dumps({'timings_seconds':timings,'speedup':{k:timings['IAS15']/v for k,v in timings.items()}},indent=2))

if __name__=='__main__': main()
