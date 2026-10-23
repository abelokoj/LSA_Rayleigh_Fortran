"""
plot_followup_ch6.py
====================
Chapter 6 figures from the October 2026 follow-up runs: the confined column (wall at R_w) and the viscous ambient (gas-to-liquid viscosity ratio), read from a results snapshot's analysis/summary.csv.

    python plot_followup_ch6.py "<snapshot>/analysis/summary.csv" --outdir <dir>

Measured rates are the interface-amplitude fits (rate_interface) and the narrow-band level-set DMD (dmd_levelset_median). The reference everywhere is the two-fluid viscous relation of Tomotika (1935) for each run's own parameters (visc_two_fluid). The theory curves are read from tables written by the tested Fortran code, `run_lsa theory`, with the campaign's fluid properties passed explicitly (r0 = sigma = rho_l = 1, mu_l = 0.02, rho_g = 1.225e-3, mu_g = 2.6e-4, i.e. Oh = 0.02, mu_2/mu_1 = 0.013): theory_vs_k.csv (k r0 = 0.05 to 0.995 in steps of 0.005) and theory_vs_mu_ratio.csv (k r0 = 0.7, 44 viscosity ratios from 0.01 to 1.3, including the four that were run). They agree with the visc_two_fluid column of summary.csv to within 5e-10 at every sampled point, which the script checks before plotting. The confined curves (no-slip wall at r = Rw, the condition the solver applies at the outer boundary, hi_bc = 5) come from confined_vs_k_Rw1p5.csv, confined_vs_Rw_k0p7.csv and confined_vs_mu_Rw4.csv, written by make_confined_tables.py with confined_dispersion.py, which reproduces the Fortran unconfined values to 5e-13 and is checked again there. Panels (b) show departures from the unconfined relation, so that the measured points and the confined theory are compared on the same axis.

    python plot_followup_ch6.py SUMMARY --theory THEORY_DIR --outdir OUTDIR
"""
import argparse, csv, os
import numpy as np
import matplotlib.pyplot as plt
from pub_style import apply_style, savefig_all, cm_both, cm_y, PALETTE


def load(path):
    rows = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            r = {k.strip(): v.strip() for k, v in r.items()}
            num = lambda key: float(r[key]) if r.get(key) else np.nan
            rows[r["label"]] = dict(k=num("kr0"), s=num("rate_interface"), ls=num("dmd_levelset_median"), th=num("visc_two_fluid"), ch=num("visc_exact"), ef=num("efolds"))
    return rows

def load_theory(path, key):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    x = np.array([float(r[key]) for r in rows])
    cols = {c: np.array([float(r[c]) for r in rows]) for c in rows[0] if c != key}
    return x, cols

def check_theory(R, tk, tm, tol=1e-8):
    """Stop if the tabulated theory disagrees with the reference rate used in the comparisons."""
    k, ck = tk
    for lab, r in R.items():
        if lab.startswith(("sw0", "conf", "dom06")) and np.isfinite(r["th"]):
            i = np.argmin(abs(k - r["k"]))
            if abs(k[i] - r["k"]) < 1e-9 and abs(ck["two_fluid_viscous"][i] - r["th"]) > tol:
                raise SystemExit(f"theory table disagrees with summary for {lab}")
    m, cm = tm
    for lab, mu in (("sw070_16", 0.013), ("visc010_16", 0.1), ("visc030_16", 0.3), ("visc100_16", 1.0)):
        i = np.argmin(abs(m - mu))
        if lab in R and abs(cm["two_fluid_viscous"][i] - R[lab]["th"]) > tol:
            raise SystemExit(f"theory table disagrees with summary for {lab}")

def dev(row, key="s"):
    return 100.0 * (row[key] / row["th"] - 1.0)

def confinement(R, tk, outdir):
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.7), constrained_layout=True)
    k, ck = tk
    a.plot(k, ck["two_fluid_viscous"], color=PALETTE[3], lw=0.9, label=r"theory, unconfined")
    kc, cc = TH["k15"]
    a.plot(kc, cc["confined"], color=PALETTE[1], lw=0.9, ls="--", label=r"theory, $R_w = 1.5R_0$")
    sweep = [l for l in ("sw020_16", "sw030_16", "sw040_16", "sw050_16", "sw060_16", "sw070_16", "sw080_16") if l in R]
    a.plot([R[l]["k"] for l in sweep], [R[l]["s"] for l in sweep], "o", mfc="none", color=PALETTE[0], ms=4.5, label=r"$R_w = 4R_0$")
    conf = [l for l in ("conf15k050_16", "conf15k060_16", "conf15_16", "conf15k080_16") if l in R]
    a.plot([R[l]["k"] for l in conf], [R[l]["s"] for l in conf], "s", color=PALETTE[1], ms=4.0, label=r"$R_w = 1.5R_0$")
    if "conf15_32" in R:
        a.plot([R["conf15_32"]["k"]], [R["conf15_32"]["s"]], "x", color=PALETTE[3], ms=6.0, mew=1.0, label=r"$R_w = 1.5R_0$, 32 cells")
    a.set_xlim(0.15, 1.0); a.set_ylim(0.10, 0.35)
    a.set_xlabel(r"$kR_0$"); a.set_ylabel(r"$s\,t_\sigma$")
    a.legend(frameon=False, loc="lower center"); cm_both(a)
    a.text(0.03, 0.95, r"(a)", transform=a.transAxes, va="top")

    walls = [("conf15_16", 1.5), ("conf20_16", 2.0), ("conf30_16", 3.0), ("sw070_16", 4.0), ("dom06_16", 6.0)]
    walls = [(l, w) for l, w in walls if l in R]
    b.axhline(0.0, color="0.5", lw=0.6)
    rw, cr = TH["rw"]; sel = np.isclose(cr["mu_ratio"], 0.013)
    b.plot(rw[sel], 100 * (cr["confined"][sel] / cr["unconfined"][sel] - 1), color=PALETTE[3], lw=0.9, label=r"theory, wall at $R_w$")
    b.plot([w for _, w in walls], [dev(R[l]) for l, _ in walls], "o", color=PALETTE[0], ms=4.0, label=r"interface fit")
    wl = [(l, w) for l, w in walls if np.isfinite(R[l]["ls"])]
    b.plot([w for _, w in wl], [dev(R[l], "ls") for l, _ in wl], "s", color=PALETTE[2], mfc="none", ms=4.0, label=r"level-set DMD")
    if "conf15_32" in R:
        b.plot([1.5], [dev(R["conf15_32"])], "x", color=PALETTE[3], ms=6.0, mew=1.0, label=r"interface fit, 32 cells")
    b.set_xlim(1.25, 6.25); b.set_ylim(-3.75, 0.5)
    b.set_xlabel(r"$R_w/R_0$"); b.set_ylabel(r"departure from unconfined theory (%)")
    b.legend(frameon=False, loc="lower right"); cm_both(b)
    b.text(0.03, 0.95, r"(b) $kR_0 = 0.7$", transform=b.transAxes, va="top")
    return savefig_all(fig, "ch6_confinement", outdir=outdir, close=True)

def viscous(R, tm, outdir):
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.7), constrained_layout=True)
    cases = [("sw070_16", 0.013), ("visc010_16", 0.1), ("visc030_16", 0.3), ("visc100_16", 1.0)]
    cases = [(l, m) for l, m in cases if l in R]
    m = [m for _, m in cases]
    mu, cm = tm
    a.plot(mu, cm["chandrasekhar"], color="0.5", lw=0.8, ls="--", label=r"Chandrasekhar (no gas)")
    a.plot(mu, cm["two_fluid_viscous"], "-", color=PALETTE[3], lw=0.9, label=r"two-fluid viscous, unconfined")
    mc, c4 = TH["mu4"]
    a.plot(mc, c4["confined"], "--", color=PALETTE[0], lw=0.9, label=r"two-fluid viscous, $R_w = 4R_0$")
    a.plot(m, [R[l]["s"] for l, _ in cases], "o", mfc="none", color=PALETTE[0], ms=4.5, label=r"$R_w = 4R_0$")
    for lab, w, mk, c in (("visc100r6_16", 6, "^", PALETTE[2]), ("visc100r8_16", 8, "s", PALETTE[1])):
        if lab in R:
            a.plot([1.0], [R[lab]["s"]], mk, color=c, ms=4.0, label=rf"$R_w = {w}R_0$")
    a.set_xscale("log"); a.set_xlim(0.01, 1.3); a.set_ylim(0.295, 0.337)
    a.set_xlabel(r"$\mu_2/\mu_1$"); a.set_ylabel(r"$s\,t_\sigma$")
    a.legend(frameon=False, loc="lower left"); cm_y(a)
    a.text(0.03, 0.97, r"(a) $kR_0 = 0.7$", transform=a.transAxes, va="top")

    walls = [(l, w) for l, w in (("visc100_16", 4.0), ("visc100r6_16", 6.0), ("visc100r8_16", 8.0)) if l in R]
    b.axhline(0.0, color="0.5", lw=0.6)
    rw, cr = TH["rw"]; sel = np.isclose(cr["mu_ratio"], 1.0) & (rw >= 3.5)
    b.plot(rw[sel], 100 * (cr["confined"][sel] / cr["unconfined"][sel] - 1), color=PALETTE[3], lw=0.9, label=r"theory, wall at $R_w$")
    b.plot([w for _, w in walls], [dev(R[l]) for l, _ in walls], "o", color=PALETTE[0], ms=4.0, label=r"interface fit")
    b.plot([w for _, w in walls], [dev(R[l], "ls") for l, _ in walls], "s", color=PALETTE[2], mfc="none", ms=4.0, label=r"level-set DMD")
    b.set_xlim(3.5, 8.5); b.set_ylim(-2.6, 0.5)
    b.set_xlabel(r"$R_w/R_0$"); b.set_ylabel(r"departure from unconfined theory (%)")
    b.legend(frameon=False, loc="lower right"); cm_both(b)
    b.text(0.03, 0.95, r"(b) $\mu_2/\mu_1 = 1$", transform=b.transAxes, va="top")
    return savefig_all(fig, "ch6_viscous_ambient", outdir=outdir, close=True)

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    p.add_argument("summary"); p.add_argument("--theory", required=True, help="folder holding theory_vs_k.csv and theory_vs_mu_ratio.csv")
    p.add_argument("--outdir", default="figs")
    a = p.parse_args()
    apply_style()
    R = load(a.summary)
    tk = load_theory(os.path.join(a.theory, "theory_vs_k.csv"), "kr0")
    tm = load_theory(os.path.join(a.theory, "theory_vs_mu_ratio.csv"), "mu_ratio")
    check_theory(R, tk, tm)
    global TH
    TH = {"k15": load_theory(os.path.join(a.theory, "confined_vs_k_Rw1p5.csv"), "kr0"),
          "rw": load_theory(os.path.join(a.theory, "confined_vs_Rw_k0p7.csv"), "Rw"),
          "mu4": load_theory(os.path.join(a.theory, "confined_vs_mu_Rw4.csv"), "mu_ratio")}
    for path in confinement(R, tk, a.outdir) + viscous(R, tm, a.outdir):
        print("wrote", path)
