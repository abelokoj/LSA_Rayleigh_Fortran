"""
compare_confined.py
===================
Compares every measured interface rate with the two-fluid viscous relation without a wall and with the no-slip wall each run actually had (confined_dispersion.py), and writes the table to CSV.

    python compare_confined.py SUMMARY_CSV OUT_CSV
"""
import os
import csv, sys
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from confined_dispersion import growth_rate
S=sys.argv[1]
R={}
for r in csv.DictReader(open(S)):
    r={a.strip():b.strip() for a,b in r.items()}; R[r['label']]=r
cases=[('sw040_16',4,2.6e-4),('sw050_16',4,2.6e-4),('sw060_16',4,2.6e-4),('sw070_16',4,2.6e-4),('sw080_16',4,2.6e-4),('dom06_16',6,2.6e-4),
       ('conf30_16',3,2.6e-4),('conf20_16',2,2.6e-4),('conf15_16',1.5,2.6e-4),('conf15_32',1.5,2.6e-4),
       ('conf15k050_16',1.5,2.6e-4),('conf15k060_16',1.5,2.6e-4),('conf15k080_16',1.5,2.6e-4),
       ('visc010_16',4,0.002),('visc030_16',4,0.006),('visc100_16',4,0.02),('visc100r6_16',6,0.02),('visc100r8_16',8,0.02)]
out=[]
print(f"{'case':14}{'k':>5}{'Rw':>5}{'mu2/mu1':>8}{'unconf.':>10}{'confined':>10}{'wall effect':>12}{'measured':>10}{'vs unconf.':>11}{'vs confined':>12}{'LS DMD vs conf.':>16}")
for l,Rw,mu2 in cases:
    r=R[l]; k=float(r['kr0']); s=float(r['rate_interface']); su=float(r['visc_two_fluid']); ls=r['dmd_levelset_median']
    sc=growth_rate(k,0.02,1.225e-3,mu2,Rw)
    lsd='' if not ls else f"{100*(float(ls)/sc-1):+.2f}"
    print(f"{l:14}{k:5.2f}{Rw:5}{mu2/0.02:8.3f}{su:10.6f}{sc:10.6f}{100*(sc/su-1):+11.3f}%{s:10.6f}{100*(s/su-1):+10.2f}%{100*(s/sc-1):+11.2f}%{lsd:>16}")
    out.append(dict(label=l,kr0=k,Rw=Rw,mu_ratio=mu2/0.02,unconfined=su,confined=sc,measured=s,dev_unconfined=100*(s/su-1),dev_confined=100*(s/sc-1),ls_dmd=ls))
with open(sys.argv[2],'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
