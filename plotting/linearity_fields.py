"""
linearity_fields.py
===================
A stronger linearity test than the agreement of scaled norms (proportionality_mod.f90). Equal norms do not imply equal directions, so here I compare the twin-difference fields themselves.

For runs that differ only in the amplitude eps, each snapshot is divided by its eps, q_i(t) = dS_i(t) / eps_i. At every time common to all runs I compare each q_i with the smallest-amplitude run q_ref:

    rel_diff(t) = ||q_i - q_ref|| / ||q_ref||       (all differences, size and direction together)
    angle(t)    = arccos( <q_i, q_ref> / (||q_i|| ||q_ref||) )   in degrees (direction only)
    norm_ratio  = ||q_i|| / ||q_ref|| - 1           (size only; this is what the norm test measures)

Linear response means rel_diff, angle and norm_ratio all stay small and do not grow with time inside the admissible window. A difference that grows with the amplitude reached, epsilon * exp(s t), is the quadratic remainder of Appendix A.1; one that is large at early times and decays is round-off or reconstruction noise.

The norms use the axisymmetric volume weight (rows times sqrt(r)), as the decomposition does, unless the export is already weighted. Inputs are the lossless exports written by `run_lsa export` (perturbation.bin, perturbation.txt, times.txt), which are kept on the cluster.

    python3 linearity_fields.py --tmax T  EPS1=EXPORT_DIR1  EPS2=EXPORT_DIR2  [EPS3=EXPORT_DIR3 ...]

Example for the 16-cell set at k R0 = 0.7 (run on the cluster from LSA_Rayleigh_Fortran):

    python3 plotting/linearity_fields.py --tmax 13.0 \
        2.5e-3=runs/analysis/sw070_16/export 5e-3=runs/analysis/lin050_16/export 1e-2=runs/analysis/lin100_16/export

Use the same --tmax as the norm test (the time at which the largest amplitude reaches a/R0 = 0.1) and repeat with the level-set exports (export_ls). Prints one table per run and writes the table to --csv (default linearity_fields.csv). Needs numpy only.
"""
import argparse, os, sys
import numpy as np


def keyvals(path):
    d = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def load(d, band=0.0, r0=1.0):
    info = keyvals(os.path.join(d, "perturbation.txt"))
    ns, nt = int(info["n_space"]), int(info["n_snap"])
    nx, ny, nc = (int(v) for v in info["grid_nx_ny_ncomp"].split())
    lo_r, lo_z, hi_r, hi_z = (float(v) for v in info["domain_lo_r_lo_z_hi_r_hi_z"].split())
    X = np.fromfile(os.path.join(d, "perturbation.bin"), "<f8").reshape((ns, nt), order="F")
    t = np.loadtxt(os.path.join(d, "times.txt"), comments="#", ndmin=1)
    t = t if t.ndim == 1 else t[:, -1]
    if info.get("weight", "none") == "none":
        r = (np.arange(nx) + 0.5) * (hi_r - lo_r) / nx + lo_r
        w = np.sqrt(np.tile(np.tile(r, ny), nc))
        X = X * w[:, None]
    if band > 0:
        # Keep only the cells within band cells of the unperturbed interface r = r0, as the narrow-band level-set
        # decomposition does. Far from the interface the level-set difference is not the displacement of a signed
        # distance but whatever the solver's redistancing leaves there, which need not scale with the amplitude.
        dr = (hi_r - lo_r) / nx
        rc = lo_r + (np.arange(nx) + 0.5) * dr
        keep = np.tile(np.tile(np.abs(rc - r0) < band * dr, ny), nc)
        X = X[keep]
    return X, t[:nt]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tmax", type=float, default=np.inf)
    ap.add_argument("--csv", default="linearity_fields.csv", help="output table (default linearity_fields.csv)")
    ap.add_argument("--band", type=float, default=0.0, help="keep only cells with |r - r0| < BAND cells (0: whole domain); use 2 for the level set, as the narrow-band decomposition does")
    ap.add_argument("runs", nargs="+", help="EPS=EXPORT_DIR")
    a = ap.parse_args()
    runs = []
    for item in a.runs:
        e, d = item.split("=", 1)
        X, t = load(d, a.band)
        runs.append((float(e), X, t, d))
    runs.sort(key=lambda r: r[0])
    eref, Xref, tref, _ = runs[0]
    rows = []
    for eps, X, t, d in runs[1:]:
        if X.shape[0] != Xref.shape[0]:
            sys.exit(f"{d}: n_space {X.shape[0]} differs from the reference {Xref.shape[0]}")
        print(f"\n== eps = {eps:g} against eps = {eref:g}   ({d})")
        print(f"{'t':>10} {'rel_diff %':>11} {'angle deg':>10} {'norm_ratio %':>13}")
        for j, tj in enumerate(t):
            if tj > a.tmax:
                break
            i = np.flatnonzero(np.abs(tref - tj) <= 1e-8 * max(1.0, abs(tj)))
            if i.size == 0:
                continue
            q, qr = X[:, j] / eps, Xref[:, i[0]] / eref
            nr = np.linalg.norm(qr)
            if nr == 0.0:
                continue
            rel = np.linalg.norm(q - qr) / nr
            cosang = np.clip(np.dot(q, qr) / (np.linalg.norm(q) * nr), -1.0, 1.0)
            ang = np.degrees(np.arccos(cosang))
            ratio = np.linalg.norm(q) / nr - 1.0
            rows.append((eps, tj, rel, ang, ratio))
            print(f"{tj:10.4f} {100 * rel:11.4f} {ang:10.4f} {100 * ratio:13.4f}")
        sel = [r for r in rows if r[0] == eps]
        if sel:
            m = np.array(sel)
            print(f"   max rel_diff {100 * m[:, 2].max():.3f} %, max angle {m[:, 3].max():.3f} deg, max |norm_ratio| {100 * np.abs(m[:, 4]).max():.3f} %")
    with open(a.csv, "w") as f:
        f.write("eps,t,rel_diff,angle_deg,norm_ratio\n")
        for r in rows:
            f.write(",".join(f"{v:.10g}" for v in r) + "\n")


if __name__ == "__main__":
    main()
