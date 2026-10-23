"""
make_dmd_evidence_figures.py
============================
I use this script to draw the two figures of Section 5.12 that show the decomposition itself, not only the growth rate it returns:

  ch5_dmd_evidence   k R0 = 0.7, 16 cells per radius: (a) singular values of the velocity and level-set snapshot matrices, (b) DMD eigenvalues in the complex plane at rank 8 and at full rank, with the exact value exp(s tau), (c) growth rate against truncation rank for every window in the convergence table, (d) residual of the reduced map against rank.
  ch5_mode_shapes    leading velocity DMD mode (rank 8) at k R0 = 0.4 to 0.7 against the two-fluid viscous eigenfunction with the no-slip wall at 4 R0 and without it.

Inputs are the files the analysis code writes for each run (dmd/, dmd_ls/, studies/modes/, converge.csv) in a results snapshot, and the confined two-fluid relation of confined_dispersion.py (Appendix D). Fonts are Computer Modern to match the thesis.

    python make_dmd_evidence_figures.py ANALYSIS_DIR CONFINED_DISPERSION_DIR OUTDIR
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A, CD, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, CD)
import confined_dispersion as cd

plt.rcParams.update({"font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
                     "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
                     "font.size": 10, "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
OH, RHO2, MU2, RW = 0.02, 1.225e-3, 2.6e-4, 4.0
GARNET, BLUE, GREY, GOLD = "#782F40", "#1f4e9c", "#666666", "#B39B5A"

def keyvals(p):
    d = {}
    for line in open(p):
        if "=" in line:
            k, v = line.split("=", 1); d[k.strip()] = v.strip()
    return d

def eigs(d):
    return np.loadtxt(os.path.join(d, "eigenvalues.txt"), comments="#", ndmin=2)

def svals(d):
    return np.loadtxt(os.path.join(d, "singular_values.txt"), comments="#", ndmin=2)

def leading_mode(d):
    info = keyvals(os.path.join(d, "modes.txt")); summ = keyvals(os.path.join(d, "summary.txt"))
    nsp, nm = int(info["n_space"]), int(info["n_modes"])
    nx, ny, nc = (int(v) for v in info["grid_nx_ny_ncomp"].split())
    m = np.fromfile(os.path.join(d, "modes.bin"), "<f8").reshape(nm, 2, nsp)
    phi = m[0, 0] + 1j * m[0, 1]
    comps = [phi[c * nx * ny:(c + 1) * nx * ny].reshape(ny, nx) for c in range(nc)]
    return comps, nx, ny, summ.get("weight", "none")

# ---------------------------------------------------------------- figure 1
case = os.path.join(A, "sw070_16")
k = 0.7
s_exact = cd.growth_rate(k, OH, RHO2, MU2, RW)
dt = float(keyvals(os.path.join(case, "dmd", "summary.txt"))["dt"])
fig, ax = plt.subplots(2, 2, figsize=(7.2, 6.6))

a = ax[0, 0]
sv_u = svals(os.path.join(case, "dmd")); sv_l = svals(os.path.join(case, "dmd_ls"))
a.semilogy(sv_u[:, 0], sv_u[:, 2], "o", ms=3, color=GARNET, label="velocity")
a.semilogy(sv_l[:, 0], sv_l[:, 2], "s", ms=3, mfc="none", color=BLUE, label="level set (band)")
a.axvline(8, color=GREY, lw=0.8, ls=":")
a.set_xlabel(r"index $i$"); a.set_ylabel(r"$\sigma_i/\sigma_1$"); a.set_title("(a) singular values", loc="left", fontsize=10)
a.legend(frameon=False, fontsize=8)

a = ax[0, 1]
th = np.linspace(0, 2 * np.pi, 400); a.plot(np.cos(th), np.sin(th), color=GREY, lw=0.6)
ef = eigs(os.path.join(case, "dmd")); e8 = eigs(os.path.join(case, "studies", "modes"))
for e, mk, lab, col in ((ef, "x", "full rank (49)", GREY), (e8, "o", "rank 8", GARNET)):
    mu = e[:, 3] * np.exp(1j * e[:, 4])
    a.plot(mu.real, mu.imag, mk, ms=4, mfc="none" if mk == "o" else None, color=col, label=lab)
mu_ex = np.exp(s_exact * dt)
a.plot([mu_ex], [0], "*", ms=9, color=GOLD, label=r"exact $e^{s\tau}$")
a.set_xlim(0.6, 1.25); a.set_ylim(-0.45, 0.45); a.set_aspect("equal")
a.set_xlabel(r"Re $\mu$"); a.set_ylabel(r"Im $\mu$"); a.set_title(r"(b) eigenvalues $\mu_j$", loc="left", fontsize=10)
a.legend(frameon=False, fontsize=7, loc="lower left")

rows = np.genfromtxt(os.path.join(case, "converge.csv"), delimiter=",", names=True)
rows = rows[rows["stride"] == 1]          # stride study shown separately (Table 5.6)
a = ax[1, 0]; b = ax[1, 1]
wins = sorted(set(zip(rows["start"].astype(int), rows["n_snap"].astype(int))))
cmap = plt.get_cmap("tab20")
handles = []
for iw, (st, n) in enumerate(wins):
    sel = (rows["start"] == st) & (rows["n_snap"] == n)
    rk = rows["rank"][sel]; g = rows["growth"][sel]; res = rows["residual"][sel]
    o = np.argsort(rk)
    adopted = (st == 111 and n == 50)
    col = GARNET if adopted else cmap(iw % 20)
    kw = dict(color=col, lw=2.0 if adopted else 0.8, marker="o", ms=4 if adopted else 2.5, zorder=3 if adopted else 2)
    h, = a.plot(rk[o], 100 * (g[o] / s_exact - 1), **kw)
    b.semilogy(rk[o], res[o], **kw)
    handles.append((h, f"start {st}, {n} snapshots" + (" (used)" if adopted else "")))
a.axhline(0, color="k", lw=0.6); a.set_xscale("log")
a.set_xlabel("truncation rank"); a.set_ylabel(r"departure from theory (%)")
a.set_title("(c) growth rate, all windows", loc="left", fontsize=10)
b.set_xscale("log"); b.set_xlabel("truncation rank"); b.set_ylabel("reduced-map residual")
b.set_title("(d) residual", loc="left", fontsize=10)
fig.legend([h for h, _ in handles], [l for _, l in handles], loc="lower center", ncol=3, fontsize=6.5, frameon=False)
fig.tight_layout(rect=[0, 0.13, 1, 1])
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"ch5_dmd_evidence.{ext}"), dpi=300)
plt.close(fig)

# ---------------------------------------------------------------- figure 2
fig, ax = plt.subplots(2, 4, figsize=(7.4, 4.0), sharex=True)
report = []
for j, (lab, k) in enumerate((("sw040_16", 0.4), ("sw050_16", 0.5), ("sw060_16", 0.6), ("sw070_16", 0.7))):
    comps, nx, ny, weight = leading_mode(os.path.join(A, lab, "studies", "modes"))
    r = (np.arange(nx) + 0.5) * RW / nx; z = (np.arange(ny) + 0.5) * (2 * np.pi / k) / ny
    proj = np.exp(-1j * k * z)[:, None]
    pr, pz = [2 * np.mean((c / np.sqrt(r)[None, :] if weight == "volume" else c) * proj, axis=0) for c in comps[:2]]
    mm = np.concatenate([pr, pz])
    out = {}
    for name, rw in (("wall", RW), ("unbounded", None)):
        s = cd.growth_rate(k, OH, RHO2, MU2, rw); ur, uz, _ = cd.eigenfunction(r, k, s, OH, RHO2, MU2, rw)
        t = np.concatenate([ur, 1j * uz]); al = np.vdot(t, mm) / np.vdot(t, t); fit = al * t
        peak = np.max(np.abs(fit)); rms = 100 * np.sqrt(np.mean(np.abs(mm - fit) ** 2)) / peak
        out[name] = (fit, al, peak, rms)
    fit, al, peak, rms = out["wall"]; fu = out["unbounded"][0]
    ph = np.exp(-1j * np.angle(al))           # rotate so that the fitted theory is real
    # u_z is a quarter period out of phase with u_r (u_z = i W), so its panel shows the imaginary part
    for i, (dat, th_w, th_u, nm, q) in enumerate(((pr, fit[:nx], fu[:nx], r"$\hat u_r$", 1.0), (pz, fit[nx:], fu[nx:], r"$\hat u_z / i$", 1j))):
        a = ax[i, j]
        a.plot(r, (th_w * ph / q).real / peak, color=GARNET, lw=1.4, label=r"theory, wall at $4R_0$")
        a.plot(r, (th_u * ph / q).real / peak, color=GREY, lw=1.0, ls="--", label="theory, unbounded")
        a.plot(r[::3], (dat * ph / q).real[::3] / peak, "o", ms=2.4, mfc="none", color=BLUE, label="DMD mode")
        a.axvline(1, color="k", lw=0.5, ls=":")
        if j == 0: a.set_ylabel(nm)
        if i == 1: a.set_xlabel(r"$r/R_0$")
    ax[0, j].set_title(fr"$kR_0={k}$  ({rms:.2f}% rms)", fontsize=9)
    report.append((lab, k, out["wall"][3], out["unbounded"][3]))
h, l = ax[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8)
fig.tight_layout(rect=[0, 0.06, 1, 1])
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"ch5_mode_shapes.{ext}"), dpi=300)
plt.close(fig)
for lab, k, w, u in report:
    print(f"{lab} k={k}: rms with wall {w:.2f} %, unbounded {u:.2f} %")
print("exact confined rate at k=0.7:", s_exact, " exp(s tau) =", np.exp(s_exact * dt))
