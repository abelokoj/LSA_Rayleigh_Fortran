"""
weighting_test.py
=================
Does the weighting between velocity and level set change the DMD growth rate? (Table 6.3 of the dissertation)

For each case, stacks the velocity twin difference (export/) and the level-set twin difference in the narrow band |r - R0| < 2 dx (export_ls/) into one observation vector, weighted by the energy norm of Chapter 4 of the dissertation (eq:4_energy_norm),

    ||q'||^2 = sum_cells V_i |u'|^2  +  lambda * sum_band (V_i / delta) phi'^2,      V_i ∝ r_i,  delta = 4 dx,

and runs the same DMD for lambda = 0 (velocity only), 0.01, 0.1, 1, 10, 100 and the level set alone. The snapshot window (skip, number of snapshots) is read from the case's own dmd_ls/summary.txt, so it is the window the pipeline used.

DMD: snapshots X1 = [x_0 .. x_{N-2}], X2 = [x_1 .. x_{N-1}]; rank-r SVD of X1; reduced operator U* X2 V S^-1; eigenvalues mu_j, rates ln(mu_j)/dt. Mode selection as in the thesis (Ch. 4): among growing, near-stationary eigenvalues (|arg mu| < 1e-2), the one with the largest amplitude b_j ||phi_j|| (b = least-squares fit of x_0). Reported: the median and interquartile range of that rate over ranks 2..12.

    python3 weighting_test.py ANALYSIS_DIR CASE [CASE ...]

ANALYSIS_DIR is runs/analysis on the cluster. Needs numpy only.
"""
import os, sys
import numpy as np

def keyvals(path):
    d = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1); d[k.strip()] = v.strip()
    return d

def load_export(d):
    info = keyvals(os.path.join(d, "perturbation.txt"))
    ns, nt = int(info["n_space"]), int(info["n_snap"])
    nx, ny, nc = (int(v) for v in info["grid_nx_ny_ncomp"].split())
    lo_r, lo_z, hi_r, hi_z = (float(v) for v in info["domain_lo_r_lo_z_hi_r_hi_z"].split())
    if info.get("weight", "none") != "none":
        sys.exit(f"{d}: expected an unweighted export, found weight = {info['weight']}")
    X = np.fromfile(os.path.join(d, "perturbation.bin"), "<f8").reshape((ns, nt), order="F")
    return X, nx, ny, nc, hi_r, float(info["dt"])

def dmd_rate(X, dt, rank):
    X1, X2 = X[:, :-1], X[:, 1:]
    U, s, Vh = np.linalg.svd(X1, full_matrices=False)
    r = min(rank, int(np.sum(s > s[0] * 1e-12)))
    U, s, Vh = U[:, :r], s[:r], Vh[:r]
    A = U.conj().T @ X2 @ Vh.conj().T / s
    mu, W = np.linalg.eig(A)
    Phi = X2 @ Vh.conj().T / s @ W                       # exact DMD modes
    b = np.linalg.lstsq(Phi, X[:, 0], rcond=None)[0]
    amp = np.abs(b) * np.linalg.norm(Phi, axis=0)
    ok = (np.abs(mu) > 1) & (np.abs(np.angle(mu)) < 1e-2)
    if not ok.any():
        return np.nan
    j = np.argmax(np.where(ok, amp, -1))
    return np.log(np.abs(mu[j])) / dt

def rank_scan(X, dt, ranks=range(2, 13)):
    v = np.array([dmd_rate(X, dt, r) for r in ranks])
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    return med, q3 - q1

def case(adir, label, lambdas=(0.0, 0.01, 0.1, 1.0, 10.0, 100.0)):
    Xu, nx, ny, nc, rmax, dt = load_export(os.path.join(adir, label, "export"))
    Xl, nxl, nyl, ncl, rmaxl, dtl = load_export(os.path.join(adir, label, "export_ls"))
    assert (nx, ny, rmax) == (nxl, nyl, rmaxl) and abs(dt - dtl) < 1e-12 and Xu.shape[1] == Xl.shape[1]
    summ = keyvals(os.path.join(adir, label, "dmd_ls", "summary.txt"))
    skip, n = int(summ["skip"]), int(summ["n_snapshots"])
    sl = slice(skip, skip + n)
    dr = rmax / nx
    r = (np.arange(nx) + 0.5) * dr                       # cell centres; layout component-major, then z, then r fastest
    Rr = np.tile(r, ny)
    wu = np.sqrt(np.tile(Rr, nc))                        # sqrt of V_i (constant factors 2 pi dr dz dropped)
    band = np.abs(Rr - 1.0) < 2 * dr
    delta = 4 * dr
    wl = np.sqrt(Rr[band] / delta)
    U = Xu[:, sl] * wu[:, None]
    L = Xl[band][:, sl] * wl[:, None]
    out = []
    for lam in lambdas:
        X = U if lam == 0 else np.vstack([U, np.sqrt(lam) * L])
        out.append((f"lambda = {lam:g}" if lam else "velocity only", *rank_scan(X, dt)))
    out.append(("level set only (band)", *rank_scan(L, dt)))
    # relative size of the two blocks at lambda = 1, to show what the weighting does
    ratio = np.linalg.norm(L) / np.linalg.norm(U)
    return out, ratio, skip, n, dt

if __name__ == "__main__":
    adir = sys.argv[1]
    for label in sys.argv[2:]:
        out, ratio, skip, n, dt = case(adir, label)
        print(f"== {label}: snapshots {skip}..{skip + n - 1} (dt = {dt:.6g}); ||level-set block|| / ||velocity block|| at lambda = 1: {ratio:.3g}")
        for name, med, iqr in out:
            print(f"   {name:24s} rate {med:.6f}   IQR {iqr:.2e}")
