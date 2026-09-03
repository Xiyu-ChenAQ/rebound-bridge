#!/usr/bin/env python3
"""Plot the real C EOS LF(4,2) accuracy summary CSV."""
import argparse, csv
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    p=argparse.ArgumentParser(); p.add_argument('--csv',required=True); p.add_argument('--out',default='accuracy.png'); a=p.parse_args()
    with open(a.csv, newline='') as f: row=next(csv.reader(f))
    # Columns are defined by compare_earth_moon_whfast_bridge.c.
    bridge_moon=float(row[14]); whfast_moon=float(row[15])
    bridge_earth=float(row[20]); whfast_earth=float(row[21])
    bridge_jup=float(row[23]); whfast_jup=float(row[24])
    labels=['Moon','Earth barycenter','Jupiter barycenter']
    x=range(len(labels)); w=.36
    fig,ax=plt.subplots(figsize=(8,4.8)); ax.bar([i-w/2 for i in x],[bridge_moon,bridge_earth,bridge_jup],w,label='Bridge EOS LF(4,2)')
    ax.bar([i+w/2 for i in x],[whfast_moon,whfast_earth,whfast_jup],w,label='WHFast')
    ax.set_yscale('log'); ax.set_xticks(list(x),labels); ax.set_ylabel('Relative position error (AU)'); ax.set_title('100-year error vs IAS15 (outer step = 7 d)'); ax.legend(); fig.tight_layout(); fig.savefig(a.out,dpi=180)

if __name__=='__main__': main()
