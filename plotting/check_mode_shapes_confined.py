"""
check_mode_shapes_confined.py
=============================
Compares the leading velocity DMD mode of each run (archives written by scripts/run_studies.sh modes) with the two-fluid viscous eigenfunction, both without a wall (Tomotika 1935) and with the no-slip wall at r = Rw that the runs have (confined_dispersion.py). One complex scale factor is fitted over the whole domain; departures are rms over each region relative to the peak of the fitted theory.

    python check_mode_shapes_confined.py ANALYSIS_DIR [--rw 4]

ANALYSIS_DIR is a snapshot's analysis/ folder holding <case>/studies/modes. Uses read_archive from plot_mode_shape_archive.py (same folder) and checks first that the unconfined eigenfunction of confined_dispersion.py equals that script's.
"""
import argparse, os
import sys, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import confined_dispersion as cd
import importlib.util
spec=importlib.util.spec_from_file_location('m',os.path.join(HERE,'plot_mode_shape_archive.py')); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
ap=argparse.ArgumentParser(); ap.add_argument("analysis"); ap.add_argument("--rw",type=float,default=4.0); args=ap.parse_args()
A=args.analysis
oh,rho2,mu2,Rw=0.02,1.225e-3,2.6e-4,args.rw
# sanity: unconfined eigenfunction from this code equals the existing script's, up to scale
r=np.linspace(0.01,3.99,200); k=0.7
s=cd.growth_rate(k,oh,rho2,mu2); a,b,_=cd.eigenfunction(r,k,s,oh,rho2,mu2)
a0,b0=m.eigenfunction(r,k,s,oh,rho2,mu2/oh)
t=np.concatenate([a0,b0]); u=np.concatenate([a,b]); al=np.dot(t,u)/np.dot(t,t)
print('unconfined eigenfunction vs existing script: max rel diff', np.max(np.abs(u-al*t))/np.max(np.abs(u)))
for lab,k in (('sw040_16',0.4),('sw050_16',0.5),('sw060_16',0.6),('sw070_16',0.7),('meas32',0.7)):
    comps,nx,ny,weight,rate=m.read_archive(f"{A}/{lab}/studies/modes",0)
    r=(np.arange(nx)+0.5)*Rw/nx; z=(np.arange(ny)+0.5)*(2*np.pi/k)/ny; proj=np.exp(-1j*k*z)[:,None]
    pr,pz=[2*np.mean((c/np.sqrt(r)[None,:] if weight=='volume' else c)*proj,axis=0) for c in comps[:2]]
    res={}
    for name,rw in (('unbounded',None),('wall 4R0',Rw)):
        s=cd.growth_rate(k,oh,rho2,mu2,rw); ur,uz,cond=cd.eigenfunction(r,k,s,oh,rho2,mu2,rw)
        t=np.concatenate([ur,1j*uz]); mm=np.concatenate([pr,pz]); al=np.vdot(t,mm)/np.vdot(t,t); fit=al*t; peak=np.max(np.abs(fit))
        n=len(r)
        reg=lambda sel: 100*np.sqrt(np.mean(np.abs((mm-fit)[np.concatenate([sel,sel])])**2))/peak
        res[name]=(reg(r<=1),reg((r>1)&(r<=2)),reg(r>2),reg(r<=4.0))
    print(f"{lab:9} k={k}: " + " | ".join(f"{nm}: liquid {v[0]:.2f}%, gas 1-2 {v[1]:.2f}%, gas 2-4 {v[2]:.2f}%, all {v[3]:.2f}%" for nm,v in res.items()))
