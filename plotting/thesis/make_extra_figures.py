"""
make_extra_figures.py
=====================
I use this script to draw the figures of Chapters 5 and 6 that show the raw evidence behind the tables:

  ch5_window          interface amplitude against time with the fitted exponential and the admissible window, and the local growth rate against time, for x = 0.7, for x = 0.9 near the cut-off, and for a preliminary run at eps = 0.02 whose window lies inside the start-up transient
  ch5_proportionality interface amplitude divided by eps for the three amplitudes at x = 0.7 (16 cells), and the departure of each from the smallest amplitude
  ch5_profiles        interface profiles r(z) - R0 of the run at x = 0.7 at several times in the linear range
  ch5_mesh            departures of the three estimates at x = 0.7 against the number of cells per radius (Table 5.4)
  ch5_stride          departure of the velocity-decomposition rate against the sampling stride (Table 5.6)
  ch6_weighting       departure of the joint decomposition against the weight lambda (Table 6.3)

Inputs: the interface-amplitude records (records/*.dat), the window table of each record (analysis/records/*.csv) and the profile bundle of one results snapshot, the raw output of the weighting test, and the confined two-fluid relation of confined_dispersion.py (Appendix D). Computer Modern fonts to match the thesis.

    python make_extra_figures.py SNAPSHOT_DIR WEIGHTING_OUTPUT_TXT CONFINED_DISPERSION_DIR OUTDIR
"""
import os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SNAP, WTXT, CD, OUT = sys.argv[1:5]
sys.path.insert(0, CD)
import confined_dispersion as cd

plt.rcParams.update({"font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
                     "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
                     "font.size": 10, "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True, "legend.frameon": False})
OH, RHO2, MU2, RW = 0.02, 1.225e-3, 2.6e-4, 4.0
GARNET, BLUE, GREY, GOLD, GREEN = "#782F40", "#1f4e9c", "#666666", "#B39B5A", "#2E6B3A"


def save(fig, name):
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=300)
    plt.close(fig)


def record(label):
    d = np.loadtxt(os.path.join(SNAP, "records", f"{label}.dat"), comments="#")
    return d[:, 0], d[:, 1]


def window(label):
    import csv
    r = next(csv.DictReader(open(os.path.join(SNAP, "analysis", "records", f"{label}.csv"))))
    return float(r["t_lo"]), float(r["t_hi"]), float(r["rate"])


def local_rate(t, a, width=0.8):
    tc, sc = [], []
    la = np.log(a)
    for i in range(len(t)):
        sel = np.abs(t - t[i]) <= width / 2
        if sel.sum() >= 3 and t[i] - width / 2 >= t[0] and t[i] + width / 2 <= t[-1]:
            p = np.polyfit(t[sel], la[sel], 1)
            tc.append(t[i]); sc.append(p[0])
    return np.array(tc), np.array(sc)


# ------------------------------------------------------------ ch5_window
cases = (("sw070_16", 0.7, r"$x = 0.7$, $\varepsilon = 2.5\times 10^{-3}$", GARNET),
         ("sw090_16", 0.9, r"$x = 0.9$, $\varepsilon = 2.5\times 10^{-3}$", BLUE),
         ("res08", 0.7, r"$x = 0.7$, $\varepsilon = 0.02$ (preliminary)", GOLD))
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.2))
for lab, k, name, col in cases:
    t, a = record(lab)
    tlo, thi, rate = window(lab)
    keep = t <= thi + 1.0            # stop shortly after the linear range closes
    t, a = t[keep], a[keep]
    s_th = cd.growth_rate(k, OH, RHO2, MU2, RW)
    ax[0].semilogy(t, a, color=col, lw=1.2, label=name)
    sel = (t >= tlo) & (t <= thi)
    p = np.polyfit(t[sel], np.log(a[sel]), 1)
    tt = np.linspace(tlo, thi, 50)
    ax[0].semilogy(tt, np.exp(np.polyval(p, tt)), color="k", lw=2.4, alpha=0.35)
    ax[0].axvspan(tlo, thi, color=col, alpha=0.10)
    tc, sc = local_rate(t, a)
    ax[1].plot(tc, sc / s_th, color=col, lw=1.2)
    ax[1].axvspan(tlo, thi, color=col, alpha=0.10)
ax[0].axhline(0.1, color=GREY, lw=0.8, ls="--"); ax[0].text(0.3, 0.115, r"$a/R_0 = 0.1$", fontsize=8, color=GREY)
ax[0].set_xlabel(r"$t/t_\sigma$"); ax[0].set_ylabel(r"$a/R_0$"); ax[0].set_title("(a) interface amplitude", loc="left", fontsize=10)
ax[0].legend(fontsize=7, loc="lower right")
ax[1].axhline(1, color="k", lw=0.7); ax[1].axhspan(0.995, 1.005, color=GREY, alpha=0.15)
ax[1].set_ylim(0.0, 1.15); ax[1].set_xlabel(r"$t/t_\sigma$"); ax[1].set_ylabel(r"local rate / theory")
ax[1].set_title("(b) local growth rate", loc="left", fontsize=10)
save(fig, "ch5_window")

# ------------------------------------------------------------ ch5_proportionality
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.0))
runs = (("sw070_16", 2.5e-3, GARNET), ("lin050_16", 5e-3, BLUE), ("lin100_16", 1e-2, GOLD))
t0, a0 = record("sw070_16")
thi_big = window("lin100_16")[1]
for lab, eps, col in runs:
    t, a = record(lab)
    ax[0].semilogy(t, a / eps, color=col, lw=1.2, label=rf"$\varepsilon = {eps:g}$")
    if lab != "sw070_16":
        n = min(len(t), len(t0)); m = np.abs(t[:n] - t0[:n]) < 1e-8
        ax[1].plot(t[:n][m], 100 * ((a[:n] / eps)[m] / (a0[:n] / 2.5e-3)[m] - 1), color=col, lw=1.2, label=rf"$\varepsilon = {eps:g}$")
for x in ax:
    x.axvline(thi_big, color=GREY, lw=0.8, ls="--")
ax[0].set_xlabel(r"$t/t_\sigma$"); ax[0].set_ylabel(r"$a/\varepsilon$"); ax[0].legend(fontsize=8)
ax[0].set_title("(a) amplitude scaled by $\\varepsilon$", loc="left", fontsize=10)
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlabel(r"$t/t_\sigma$"); ax[1].set_ylabel(r"departure from $\varepsilon = 2.5\times10^{-3}$ (%)")
ax[1].set_title("(b) departure of the scaled amplitude", loc="left", fontsize=10); ax[1].legend(fontsize=8)
save(fig, "ch5_proportionality")

# ------------------------------------------------------------ ch5_profiles
def read_bundle(path):
    lines = open(path).read().split("\n")
    z, blocks, cur, mode = [], [], None, None
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith("#") or s.startswith("WAVELENGTH") or s.startswith("NZ") or s.startswith("NBLOCKS"):
            continue
        if s == "END" or s.startswith("END"):
            continue
        if s == "Z":
            mode = "z"; continue
        if s.startswith("BLOCK"):
            cur = [float(s.split()[2]), []]; blocks.append(cur); mode = "b"; continue
        vals = [float(v) for v in s.split()]
        if mode == "z":
            z.extend(vals)
        else:
            cur[1].extend(vals)
    return np.array(z), blocks

z, blocks = read_bundle(os.path.join(SNAP, "bundles", "sw070_16", "profiles.dat"))
nz = len(z)
fig, ax = plt.subplots(figsize=(5.6, 3.0))
times = (0.0, 4.0, 8.0, 10.0, 12.0)
cmap = plt.get_cmap("copper")
for i, tw in enumerate(times):
    b = min(blocks, key=lambda bl: abs(bl[0] - tw))
    r = np.array(b[1][:nz])
    ax.plot(z, (r - 1.0) * 100, color=cmap(i / (len(times) - 1) * 0.85), lw=1.3, label=rf"$t = {b[0]:.1f}\,t_\sigma$")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel(r"$z/R_0$"); ax.set_ylabel(r"$(r - R_0)/R_0$ ($\times 100$)"); ax.set_xlim(0, z[-1] + z[0])
ax.legend(fontsize=7, ncol=3, loc="upper center")
save(fig, "ch5_profiles")

# ------------------------------------------------------------ ch5_mesh  (Table 5.4)
N = np.array([8, 16, 32])
interface = np.array([0.16, -0.45, -0.43]); spread = np.array([0.92, 0.75, 0.83])
levelset = np.array([-1.00, -0.86, -0.94]); velocity = np.array([0.39, -1.68, 0.63])
fig, ax = plt.subplots(figsize=(4.4, 3.0))
ax.errorbar(N, interface, yerr=spread / 2, fmt="o-", color=GARNET, capsize=3, label=r"interface fit ($\pm$ half spread)")
ax.plot(N, levelset, "s-", color=BLUE, mfc="none", label="level-set decomposition")
ax.plot(N[[0, 2]], velocity[[0, 2]], "^", color=GREEN, mfc="none", label="velocity decomposition")
ax.plot(N[1], velocity[1], "^", color=GREEN, alpha=0.35, label="velocity, flagged")
ax.axhline(0, color="k", lw=0.6); ax.set_xscale("log", base=2); ax.set_xticks(N); ax.set_xticklabels([str(n) for n in N])
ax.set_xlabel(r"cells per radius $N_R$"); ax.set_ylabel(r"departure from Eqn. (3.10) (%)"); ax.legend(fontsize=7)
save(fig, "ch5_mesh")

# ------------------------------------------------------------ ch5_stride  (Table 5.6)
xs = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
tab = np.array([[0.52, 0.21, 0.11], [0.19, 0.41, 0.16], [0.34, 0.31, 0.72], [1.08, 0.45, 0.71], [0.38, 0.64, 0.93], [0.22, -0.55, 1.46]])
fig, ax = plt.subplots(figsize=(4.4, 3.0))
for i, x in enumerate(xs):
    ax.plot([1, 2, 4], tab[i], "o-", lw=1.0, ms=3.5, color=plt.get_cmap("viridis")(i / 5), label=rf"$x = {x}$")
ax.axhline(0, color="k", lw=0.6); ax.set_xscale("log", base=2); ax.set_xticks([1, 2, 4]); ax.set_xticklabels(["1", "2", "4"])
ax.set_xlabel("stride (stored snapshots)"); ax.set_ylabel(r"departure from Eqn. (3.10) (%)"); ax.legend(fontsize=7, ncol=2)
save(fig, "ch5_stride")

# ------------------------------------------------------------ ch6_weighting (raw output of the weighting test)
th = {"sw040_16": (0.4, RW), "sw050_16": (0.5, RW), "sw060_16": (0.6, RW), "sw070_16": (0.7, RW),
      "sw080_16": (0.8, RW), "meas16": (0.7, RW), "meas32": (0.7, RW), "conf15_32": (0.7, 1.5), "visc100r8_16": (0.7, 8.0)}
names = {"sw040_16": r"$x=0.4$", "sw050_16": r"$x=0.5$", "sw060_16": r"$x=0.6$", "sw070_16": r"$x=0.7$",
         "meas32": r"$x=0.7$, 32 cells", "conf15_32": r"$R_w = 1.5R_0$, 32 cells", "visc100r8_16": r"$\mu_2=\mu_1$, $R_w=8R_0$",
         "sw080_16": r"$x=0.8$ (flagged)", "meas16": r"$x=0.7$, 16 cells (flagged)"}
txt = open(WTXT).read()
lams = [1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3]   # end points stand for velocity only and level set only
fig, ax = plt.subplots(figsize=(6.2, 3.4))
order = ["sw040_16", "sw050_16", "sw060_16", "sw070_16", "meas32", "conf15_32", "visc100r8_16", "sw080_16", "meas16"]
cm = plt.get_cmap("tab10")
for j, lab in enumerate(order):
    blk = txt.split(f"== {lab}:")[1].split("==")[0]
    rates = [float(v) for v in re.findall(r"rate ([\d.]+)", blk)]
    k, rw = th[lab]
    mu2 = OH * 1.0 if lab == "visc100r8_16" else MU2
    s = cd.growth_rate(k, OH, RHO2, mu2, rw)
    dev = 100 * (np.array(rates) / s - 1)
    flagged = lab in ("sw080_16", "meas16")
    ax.plot(lams, dev, "o-" if not flagged else "x--", ms=3.5, lw=1.0, color=cm(j), alpha=0.5 if flagged else 1, label=names[lab])
ax.axhline(0, color="k", lw=0.6); ax.axvline(1, color=GREY, lw=0.8, ls=":")
ax.set_xscale("log"); ax.set_xticks(lams); ax.set_xticklabels(["vel. only", "0.01", "0.1", "1", "10", "100", "l.s. only"])
ax.set_ylim(-3, 3); ax.set_xlabel(r"weight $\lambda$ of the level-set block"); ax.set_ylabel("departure from theory (%)")
ax.legend(fontsize=6.5, ncol=1, loc="center left", bbox_to_anchor=(1.0, 0.5))
save(fig, "ch6_weighting")
print("ok")
