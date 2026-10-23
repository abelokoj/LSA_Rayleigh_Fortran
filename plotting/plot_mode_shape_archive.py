#!/usr/bin/env python3
"""Leading DMD velocity mode against the analytic eigenfunction, from a DMD archive (no plotfiles, no yt).

Reads modes.bin / modes.txt / summary.txt written by `run_lsa analyse --archive DIR` (or scripts/run_studies.sh modes) on the velocity fields (x_velocity = u_r, y_velocity = u_z), projects the chosen mode onto exp(-i k z) along the axis, and compares the radial profiles of u_r and u_z with

  * the exact two-fluid viscous eigenfunction (the null vector of the Tomotika 1935 boundary-condition matrix at the growing root s, the same matrix as dispersion_mod::two_fluid_viscous_det): inside u_r = A1 k I1(kr) + C1 I1(l1 r), u_z = i[k A1 I0(kr) + (l1/k) C1 I0(l1 r)]; outside u_r = -A2 k K1(kr) + C2 K1(l2 r), u_z = i[k A2 K0(kr) - (l2/k) C2 K0(l2 r)], with l^2 = k^2 + s/nu;
  * the inviscid single-fluid eigenfunction u_r ∝ I1(kr), u_z ∝ i I0(kr) inside the column.

One complex scale factor, fitted jointly to both components by least squares, aligns the mode with the theory (DMD modes have arbitrary amplitude and phase); the reported departure is the rms residual relative to the largest theoretical value, over r <= --r-max.

Usage:
    python3 plotting/plot_mode_shape_archive.py runs/analysis/sw070_16/studies/modes --kr0 0.7 --outdir figures
    python3 plotting/plot_mode_shape_archive.py ARCHIVE --kr0 0.7 --rmax 4 --zlen 8.976 --mu 0.02 --rho-g 1.225e-3 --mu-g 2.6e-4
"""
import argparse
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import ive, kve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def read_keyvals(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if "=" in line:
                k, v = line.split("=", 1)
                out[k.strip()] = v.strip()
    return out


def read_archive(d, which):
    info = read_keyvals(os.path.join(d, "modes.txt"))
    summ = read_keyvals(os.path.join(d, "summary.txt"))
    nsp = int(info["n_space"])
    nm = int(info["n_modes"])
    nx, ny, nc = (int(v) for v in info["grid_nx_ny_ncomp"].split())
    m = np.fromfile(os.path.join(d, "modes.bin"), "<f8").reshape(nm, 2, nsp)
    phi = m[which, 0] + 1j * m[which, 1]
    ev = np.loadtxt(os.path.join(d, "eigenvalues.txt"), comments="#", ndmin=2)
    rate = ev[which, 1]
    comps = [phi[c * nx * ny:(c + 1) * nx * ny].reshape(ny, nx) for c in range(nc)]
    return comps, nx, ny, summ.get("weight", "none"), rate


def ratio_i(z):
    return ive(1, z) / ive(0, z)


def ratio_k(z):
    return kve(1, z) / kve(0, z)


def tomotika_matrix(s, x, oh1, rr, mr):
    """Scaled boundary-condition matrix (columns A1, C1, A2, C2 divided by I1(x), I1(l1), K1(x), K1(l2)); capillary units."""
    mu1, mu2 = oh1, mr * oh1
    qk, pk = ratio_i(x), ratio_k(x)
    l1 = np.sqrt(x * x + s / mu1)
    q1 = ratio_i(l1)
    nu2 = mu2 / rr
    l2 = np.sqrt(x * x + s / nu2)
    p2 = ratio_k(l2)
    M = np.array([
        [x, 1.0, x, -1.0],
        [x / qk, (l1 / x) / q1, -x / pk, (l2 / x) / p2],
        [2 * mu1 * x * x, mu1 * (x + l1 * l1 / x), 2 * mu2 * x * x, -mu2 * (x + l2 * l2 / x)],
        [s / qk + 2 * mu1 * x * x * (1 / qk - 1 / x) + (x * x - 1) * x / s,
         2 * mu1 * l1 * (1 / q1 - 1 / l1) + (x * x - 1) / s,
         -rr * s / pk - 2 * mu2 * x * x * (1 / pk + 1 / x),
         2 * mu2 * l2 * (1 / p2 + 1 / l2)]])
    return M, l1, l2


def growth_rate_two_fluid(x, oh1, rr, mr):
    """Growing root of the two-fluid viscous relation (capillary units), bracketed below the two-fluid inviscid rate."""
    up = np.sqrt(x * (1 - x * x) / (1 / ratio_i(x) + rr / ratio_k(x))) * (1 + 1e-9)
    f = lambda s: np.linalg.det(tomotika_matrix(s, x, oh1, rr, mr)[0])
    grid = up * 10.0 ** (-14 * np.arange(0, 801) / 800)
    vals = [f(g) for g in grid]
    for i in range(1, len(grid)):
        if vals[i - 1] * vals[i] <= 0:
            return brentq(f, grid[i], grid[i - 1], xtol=1e-15, rtol=1e-14)
    return np.nan


def eigenfunction(r, x, s, oh1, rr, mr):
    """Radial profiles (u_r, u_z/i) of the two-fluid viscous eigenmode at radius r (R0 = 1)."""
    M, l1, l2 = tomotika_matrix(s, x, oh1, rr, mr)
    _, _, vt = np.linalg.svd(M)
    a1s, c1s, a2s, c2s = vt[-1]
    # undo the column scaling with exponentially scaled Bessel functions: A1 = a1s/I1(x) etc.
    ur = np.empty_like(r)
    uz = np.empty_like(r)
    inn = r <= 1.0
    ri, ro = r[inn], r[~inn]
    ur[inn] = (a1s * x * ive(1, x * ri) / ive(1, x) * np.exp(x * (ri - 1))
               + c1s * ive(1, l1 * ri) / ive(1, l1) * np.exp(l1 * (ri - 1)))
    uz[inn] = (a1s * x * ive(0, x * ri) / ive(1, x) * np.exp(x * (ri - 1))
               + c1s * (l1 / x) * ive(0, l1 * ri) / ive(1, l1) * np.exp(l1 * (ri - 1)))
    ur[~inn] = (-a2s * x * kve(1, x * ro) / kve(1, x) * np.exp(-x * (ro - 1))
                + c2s * kve(1, l2 * ro) / kve(1, l2) * np.exp(-l2 * (ro - 1)))
    uz[~inn] = (a2s * x * kve(0, x * ro) / kve(1, x) * np.exp(-x * (ro - 1))
                - c2s * (l2 / x) * kve(0, l2 * ro) / kve(1, l2) * np.exp(-l2 * (ro - 1)))
    return ur, uz


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("archive", help="DMD archive directory (velocity fields)")
    ap.add_argument("--kr0", type=float, required=True)
    ap.add_argument("--mode", type=int, default=1, help="mode index in the archive (1 = largest amplitude)")
    ap.add_argument("--rmax", type=float, default=4.0, help="outer radius of the domain")
    ap.add_argument("--zlen", type=float, default=None, help="axial length (default one wavelength 2 pi / kr0)")
    ap.add_argument("--mu", type=float, default=0.02, help="liquid viscosity (Oh, with rho = sigma = r0 = 1)")
    ap.add_argument("--rho-g", type=float, default=1.225e-3, dest="rho_g")
    ap.add_argument("--mu-g", type=float, default=2.6e-4, dest="mu_g")
    ap.add_argument("--r-max", type=float, default=2.0, dest="r_max", help="largest radius compared")
    ap.add_argument("--outdir", default=None, help="write mode_shape_<tag> figures here")
    a = ap.parse_args()

    comps, nx, ny, weight, rate = read_archive(a.archive, a.mode - 1)
    k = a.kr0
    zlen = a.zlen if a.zlen else 2 * np.pi / k
    r = (np.arange(nx) + 0.5) * a.rmax / nx
    z = (np.arange(ny) + 0.5) * zlen / ny
    proj = np.exp(-1j * k * z)[:, None]
    prof = []
    for c in comps[:2]:
        u = c / np.sqrt(r)[None, :] if weight == "volume" else c
        prof.append(2.0 * np.mean(u * proj, axis=0))
    pr, pz = prof

    s = growth_rate_two_fluid(k, a.mu, a.rho_g, a.mu_g / a.mu)
    ur_t, uz_t = eigenfunction(r, k, s, a.mu, a.rho_g, a.mu_g / a.mu)
    sel = r <= a.r_max
    t = np.concatenate([ur_t[sel], 1j * uz_t[sel]])
    m = np.concatenate([pr[sel], pz[sel]])
    alpha = np.vdot(t, m) / np.vdot(t, t)
    fit = alpha * t
    rel = np.sqrt(np.mean(np.abs(m - fit) ** 2)) / np.max(np.abs(fit))
    n = sel.sum()
    rel_r = np.sqrt(np.mean(np.abs(m[:n] - fit[:n]) ** 2)) / np.max(np.abs(fit[:n]))
    rel_z = np.sqrt(np.mean(np.abs(m[n:] - fit[n:]) ** 2)) / np.max(np.abs(fit[n:]))
    # inviscid single-fluid shape inside the column, for comparison
    inside = sel & (r <= 1.0)
    ti = np.concatenate([ive(1, k * r[inside]) * np.exp(k * r[inside]), 1j * ive(0, k * r[inside]) * np.exp(k * r[inside])])
    mi = np.concatenate([pr[inside], pz[inside]])
    bi = np.vdot(ti, mi) / np.vdot(ti, ti)
    rel_i = np.sqrt(np.mean(np.abs(mi - bi * ti) ** 2)) / np.max(np.abs(bi * ti))

    print(f"archive          {a.archive}")
    print(f"mode {a.mode}: growth rate {rate:.6f}; two-fluid viscous theory s = {s:.6f}")
    print(f"rms departure, r <= {a.r_max}: {100 * rel:.2f} % of the peak (u_r {100 * rel_r:.2f} %, u_z {100 * rel_z:.2f} %)")
    print(f"inviscid I1/I0 shape inside the column: {100 * rel_i:.2f} %")

    if a.outdir:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from pub_style import apply_style, cm_both, savefig_all
        apply_style()
        fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.8))
        rr_ = r[sel]
        scale = np.max(np.abs(fit[:n]))
        for ax, meas, th, lab in [(axs[0], m[:n], fit[:n], r"$u_r$"), (axs[1], m[n:], fit[n:], r"$u_z$")]:
            ph = np.exp(-1j * np.angle(th[np.argmax(np.abs(th))]))
            ax.plot(rr_, np.real(th * ph) / scale, "k-", lw=1.0, label="two-fluid viscous theory")
            ax.plot(rr_, np.real(meas * ph) / scale, "o", ms=3, mfc="none", label="DMD mode")
            ax.axvline(1.0, color="0.5", lw=0.6, ls=":")
            ax.set_xlim(0, a.r_max)
            ax.set_xlabel(r"$r/R_0$")
            ax.set_ylabel(lab + " (normalized)")
            cm_both(ax)
        axs[0].legend(frameon=False, loc="best")
        fig.tight_layout()
        savefig_all(fig, f"mode_shape_kr0_{a.kr0:g}".replace(".", "p"), outdir=a.outdir, close=True)
        print(f"figure written to {a.outdir}/")


if __name__ == "__main__":
    main()
