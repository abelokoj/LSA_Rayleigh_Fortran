#!/usr/bin/env python3
"""The full figure set of the October 2026 campaign, the successor of the September `all_figures.pdf`.

The September collection (8 cells per radius, epsilon = 0.02) was drawn by `plot_results.py`, `plot_pinchoff.py`, `plot_ohnesorge.py`, `plot_schematic.py` and `animate.py` from that campaign's runs. This script draws the same kinds of figure from a results snapshot of the October campaign (16 cells per radius, epsilon = 0.04 dx), with the fit windows of the final analysis and with both Chandrasekhar's one-fluid and the two-fluid viscous (Tomotika) rates, and then calls the older scripts for the figures that need no change. Everything is read from files on the workstation: the snapshot (records, profile bundles, analysis CSVs), the two-fluid records of 2026-10-16, and optionally a velocity twin-difference export for the field figures. No plotfile is needed.

Figures written to --outdir (PNG, PDF, SVG, house style):

    gallery_growth_<case>        interface amplitude a(t) on a log scale, the fit window of the final analysis, the fitted rate, and the slopes of Chandrasekhar's and the two-fluid viscous rates
    gallery_local_rate_<case>    local rate d(ln a)/dt over sliding windows, against the three theories
    gallery_interface_<case>     interface r(z) at six instants up to the end of the linear range
    gallery_sweep_growth         a(t)/a(0) for the whole 16-cell sweep, k r0 = 0.2 to 0.9
    gallery_linearity            a(t)/epsilon for three amplitudes at k r0 = 0.7, and their departure from the mean
    gallery_field_<case>         axial velocity of the twin difference (the perturbation) in the (r, z) plane at four instants, with the interface (needs --export)
    anim_field_<case>.gif        the same, every snapshot (needs --export)
    field_y_velocity_kr0_*       full axial velocity of the perturbed run with the interface, as in September (only if the bundle holds field slices, i.e. on the cluster)
    anim_field_full_<case>.gif   the same as an animation (cluster only)
    anim_interface_<case>.gif/.html, pinchoff_*, ohnesorge_*, schematic_twin_run, campaign_*: from the older scripts
    all_figures.pdf              every PNG above, one per page

Usage (from the top of this repository):
    python3 plotting/plot_campaign_gallery.py results/campaign_2026-10-09 --two-fluid results/two_fluid_2026-10-16/records --export results/exports/meas16 --export-case meas16 --outdir figures/gallery
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import bundle as bundle_io
import dispersion
from pub_style import apply_style, savefig_all, PALETTE
from plot_mode_shape_archive import growth_rate_two_fluid

MU, RHO_G, MU_G = 0.02, 1.225e-3, 2.6e-4      # campaign properties; rho_l = sigma = r0 = 1
SWEEP = ["sw020_16", "sw030_16", "sw040_16", "sw050_16", "sw060_16", "sw070_16", "sw080_16", "sw090_16"]
LINEARITY = [("sw070_16", 2.5e-3), ("lin050_16", 5e-3), ("lin100_16", 1e-2)]


def tex(case):
    """Case label safe for the mathtext/LaTeX renderer of the house style."""
    return case.replace("_", r"\_")


def theory(kr0):
    return {"inviscid": float(dispersion.growth_rate(kr0)),
            "Chandrasekhar": float(dispersion.growth_rate_viscous_exact(kr0, mu=MU)),
            "two-fluid viscous": float(growth_rate_two_fluid(kr0, MU, RHO_G, MU_G / MU))}


def read_record(path):
    d = np.loadtxt(path, comments="#")
    return d[:, 0], d[:, 1]


def read_fit(csv):
    """Fit window and rate of the final analysis, from the record CSV of run_lsa record."""
    with open(csv) as fh:
        head = [h.strip() for h in fh.readline().split(",")]
        row = [v.strip() for v in fh.readline().split(",")]
    d = dict(zip(head, row))
    return float(d["t_lo"]), float(d["t_hi"]), float(d["rate"])


def local_rate(t, a, width=1.0):
    tc, sc = [], []
    for i in range(len(t)):
        m = (t >= t[i] - width / 2) & (t <= t[i] + width / 2)
        if m.sum() >= 4 and t[i] - width / 2 >= t[0] and t[i] + width / 2 <= t[-1]:
            tc.append(t[i])
            sc.append(np.polyfit(t[m], np.log(a[m]), 1)[0])
    return np.array(tc), np.array(sc)


def fig_growth(case, kr0, rec, csv, outdir):
    t, a = read_record(rec)
    t_lo, t_hi, rate = read_fit(csv)
    th = theory(kr0)
    win = (t >= t_lo) & (t <= t_hi)
    c = np.polyfit(t[win], np.log(a[win]), 1)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.axvspan(t_lo, t_hi, color="0.92", zorder=0, label="fit window")
    ax.semilogy(t, a, "o", ms=2.5, mfc="none", color="k", label="measured")
    tt = np.linspace(t[0], t[-1], 200)
    ax.semilogy(tt, np.exp(np.polyval(c, tt)), "-", color=PALETTE[0], label=f"fit, $s = {rate:.4f}$")
    a0 = np.exp(np.polyval(c, t_lo))
    for (name, ls, col) in [("Chandrasekhar", "--", PALETTE[1]), ("two-fluid viscous", "-.", PALETTE[2])]:
        ax.semilogy(tt, a0 * np.exp(th[name] * (tt - t_lo)), ls, color=col, label=f"{name}, $s = {th[name]:.4f}$")
    ax.set_xlabel("$t$")
    ax.set_ylabel("interface mode amplitude $a(t)$")
    ax.set_xlim(t[0], t[-1])
    ax.set_title(f"Interface growth, $kR_0 = {kr0:g}$ ({tex(case)})")
    ax.legend(loc="lower right")
    savefig_all(fig, f"gallery_growth_{case}", outdir, close=True)


def fig_local_rate(case, kr0, rec, csv, outdir):
    t, a = read_record(rec)
    t_lo, t_hi, rate = read_fit(csv)
    tc, sc = local_rate(t, a)
    th = theory(kr0)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.axvspan(t_lo, t_hi, color="0.92", zorder=0, label="fit window")
    ax.plot(tc, sc, "-", color=PALETTE[0], label=r"local $\mathrm{d}\ln a/\mathrm{d}t$")
    for (name, ls, col) in [("inviscid", ":", "0.4"), ("Chandrasekhar", "--", PALETTE[1]), ("two-fluid viscous", "-.", PALETTE[2])]:
        ax.axhline(th[name], ls=ls, color=col, label=f"{name}, {th[name]:.4f}")
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$\mathrm{d}\ln a/\mathrm{d}t$")
    ax.set_xlim(tc[0], tc[-1])
    ax.set_ylim(0, 1.15 * th["inviscid"])
    ax.set_title(f"Local growth rate, $kR_0 = {kr0:g}$ ({tex(case)}); fit {rate:.4f}")
    ax.legend(loc="lower right")
    savefig_all(fig, f"gallery_local_rate_{case}", outdir, close=True)


def fig_interface(case, kr0, prof, outdir, n=6, bound=0.1):
    p = bundle_io.read_profiles(prof)
    z, times, R = p["z"], np.asarray(p["times"]), np.asarray(p["radius"])
    amp = np.nanmax(np.abs(R - 1.0), axis=1)
    last = np.where(amp <= bound)[0]
    last = last[-1] if last.size else len(times) - 1
    idx = np.unique(np.linspace(0, last, n).round().astype(int))
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for j, i in enumerate(idx):
        ax.plot(z, R[i], color=PALETTE[j % len(PALETTE)], label=f"$t = {times[i]:.2f}$")
    ax.set_xlabel("$z$")
    ax.set_ylabel(r"$r_\mathrm{interface}$")
    ax.set_xlim(z[0], z[-1])
    ax.set_title(f"Interface shape, $kR_0 = {kr0:g}$ ({tex(case)})")
    ax.legend(loc="upper right", ncol=2)
    savefig_all(fig, f"gallery_interface_{case}", outdir, close=True)


def fig_sweep(records, outdir):
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for j, case in enumerate(SWEEP):
        path = os.path.join(records, case + ".dat")
        if not os.path.exists(path):
            continue
        t, a = read_record(path)
        ax.semilogy(t, a / a[0], color=PALETTE[j % len(PALETTE)], label=f"$kR_0 = {int(case[2:5]) / 100:g}$")
    ax.axhline(0.1 / 2.5e-3, color="0.5", ls=":", lw=0.8)
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$a(t)/a(0)$")
    ax.set_ylim(1, 60)
    ax.set_title("Interface growth over the 16-cell sweep (dotted: $a = 0.1 R_0$)")
    ax.legend(loc="upper left", ncol=2)
    savefig_all(fig, "gallery_sweep_growth", outdir, close=True)


def fig_linearity(records, outdir):
    data = []
    for case, eps in LINEARITY:
        path = os.path.join(records, case + ".dat")
        if os.path.exists(path):
            t, a = read_record(path)
            data.append((case, eps, t, a))
    if len(data) < 2:
        return
    t_end = min(d[2][-1] for d in data)
    tg = np.linspace(0, t_end, 300)
    curves = [np.interp(tg, t, a) / eps for _, eps, t, a in data]
    mean = np.mean(curves, axis=0)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 6.0), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    for j, ((case, eps, _, _), c) in enumerate(zip(data, curves)):
        ax1.semilogy(tg, c, color=PALETTE[j], label=fr"$\epsilon = {eps:g}$ ({tex(case)})")
        ax2.plot(tg, 100 * (c / mean - 1), color=PALETTE[j])
    ax1.set_ylabel(r"$a(t)/\epsilon$")
    ax1.set_title(r"Linearity at $kR_0 = 0.7$, 16 cells per radius")
    ax1.legend(loc="upper left")
    ax2.axhline(0, color="0.5", lw=0.8)
    ax2.set_xlabel("$t$")
    ax2.set_ylabel("departure from mean (%)")
    ax2.set_xlim(0, t_end)
    savefig_all(fig, "gallery_linearity", outdir, close=True)


def read_export(d):
    info = {}
    with open(os.path.join(d, "perturbation.txt")) as fh:
        for line in fh:
            if "=" in line:
                k, v = line.split("=", 1)
                info[k.strip()] = v.strip()
    nsp, nsn = int(info["n_space"]), int(info["n_snap"])
    nx, ny, nc = (int(v) for v in info["grid_nx_ny_ncomp"].split())
    lo_r, lo_z, hi_r, hi_z = (float(v) for v in info["domain_lo_r_lo_z_hi_r_hi_z"].split())
    fields = info.get("fields", "").split()
    X = np.fromfile(os.path.join(d, "perturbation.bin"), "<f8").reshape((nsp, nsn), order="F")
    t = np.loadtxt(os.path.join(d, "times.txt"))
    comps = {fields[c] if c < len(fields) else str(c): X[c * nx * ny:(c + 1) * nx * ny].reshape(ny, nx, nsn) for c in range(nc)}
    r = lo_r + (np.arange(nx) + 0.5) * (hi_r - lo_r) / nx
    z = lo_z + (np.arange(ny) + 0.5) * (hi_z - lo_z) / ny
    return comps, r, z, t


def _interface_at(p, time):
    times = np.asarray(p["times"])
    i = int(np.argmin(np.abs(times - time)))
    return np.asarray(p["radius"])[i]


def fig_field(case, export, prof, outdir, field="y_velocity", r_max=2.0, gif=True):
    comps, r, z, t = read_export(export)
    if field not in comps:
        print(f"  export has no {field}; fields: {list(comps)}")
        return
    F = comps[field]
    p = bundle_io.read_profiles(prof) if prof and os.path.exists(prof) else None
    keep = r <= r_max
    live = np.where(np.abs(F[:, keep, :]).max(axis=(0, 1)) > 0)[0]
    idx = np.unique(np.linspace(live[0] if live.size else 0, len(t) - 1, 4).round().astype(int))
    fig, axs = plt.subplots(1, len(idx), figsize=(9.0, 5.2), sharey=True)
    for ax, i in zip(axs, idx):
        f = F[:, keep, i]
        m = np.abs(f).max() or 1.0
        im = ax.pcolormesh(r[keep], z, f / m, cmap="RdBu_r", vmin=-1, vmax=1, shading="auto", rasterized=True)
        if p is not None:
            ax.plot(_interface_at(p, t[i]), np.asarray(p["z"]), "k-", lw=1.0)
        ax.set_xlabel("$r$")
        ax.set_title(f"$t = {t[i]:.2f}$\n" + r"$\max|u_z'| = $" + f"{m:.1e}", fontsize=9)
        ax.set_xlim(0, r_max)
    axs[0].set_ylabel("$z$")
    fig.colorbar(im, ax=axs, shrink=0.85, label=r"$u_z' / \max|u_z'|$")
    fig.suptitle(f"Axial velocity of the perturbation (twin difference) with the interface, {tex(case)}", fontsize=10)
    savefig_all(fig, f"gallery_field_{case}", outdir, close=True)
    if not gif:
        return
    fig, ax = plt.subplots(figsize=(3.2, 5.6))
    frames = list(range(int(idx[0]), len(t), 2))
    f0 = F[:, keep, 0]
    im = ax.pcolormesh(r[keep], z, f0 * 0, cmap="RdBu_r", vmin=-1, vmax=1, shading="auto")
    line, = ax.plot([], [], "k-", lw=1.0)
    ax.set_xlim(0, r_max)
    ax.set_xlabel("$r$")
    ax.set_ylabel("$z$")
    title = ax.set_title("", fontsize=8)
    fig.tight_layout()
    fig.subplots_adjust(top=0.93)

    def draw(k):
        i = frames[k]
        f = F[:, keep, i]
        m = np.abs(f).max() or 1.0
        im.set_array((f / m).ravel())
        if p is not None:
            line.set_data(_interface_at(p, t[i]), np.asarray(p["z"]))
        title.set_text(f"$t = {t[i]:.2f}$, " + r"$\max|u_z'|=$" + f"{m:.1e}")
        return im, line, title

    anim = animation.FuncAnimation(fig, draw, frames=len(frames), blit=False)
    anim.save(os.path.join(outdir, f"anim_field_{case}.gif"), writer=animation.PillowWriter(fps=12), dpi=110)
    plt.close(fig)


def run(cmd):
    print("  $ " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:])


def combine(outdir, order):
    from PIL import Image
    pngs = [os.path.join(outdir, f + ".png") for f in order if os.path.exists(os.path.join(outdir, f + ".png"))]
    pages = []
    for f in pngs:
        im = Image.open(f).convert("RGB")
        s = 1600 / max(im.size)
        if s < 1:
            im = im.resize((int(im.width * s), int(im.height * s)))
        pages.append(im)
    if pages:
        pages[0].save(os.path.join(outdir, "all_figures.pdf"), save_all=True, append_images=pages[1:], resolution=150)
        print(f"  all_figures.pdf: {len(pages)} pages")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("snapshot", help="results snapshot, e.g. results/campaign_2026-10-09")
    ap.add_argument("--two-fluid", dest="two_fluid", default=None, help="record CSVs with t_lo/t_hi (results/two_fluid_2026-10-16/records); default: <snapshot>/analysis/records")
    ap.add_argument("--case", default="sw070_16", help="case for the growth, local-rate and interface figures")
    ap.add_argument("--kr0", type=float, default=0.7)
    ap.add_argument("--export", default=None, help="velocity twin-difference export directory (perturbation.bin/txt, times.txt)")
    ap.add_argument("--export-case", dest="export_case", default="meas16", help="case of the export, for the interface overlay")
    ap.add_argument("--pinch", default="pinch16", help="bundle of the break-up run")
    ap.add_argument("--no-gif", action="store_true", dest="no_gif")
    ap.add_argument("--outdir", default="figures")
    a = ap.parse_args()
    apply_style()
    os.makedirs(a.outdir, exist_ok=True)
    snap = a.snapshot
    records = os.path.join(snap, "records")
    csvdir = a.two_fluid or os.path.join(snap, "analysis", "records")
    bundles = os.path.join(snap, "bundles")
    case = a.case

    print("figures from the October data")
    fig_growth(case, a.kr0, os.path.join(records, case + ".dat"), os.path.join(csvdir, case + ".csv"), a.outdir)
    fig_local_rate(case, a.kr0, os.path.join(records, case + ".dat"), os.path.join(csvdir, case + ".csv"), a.outdir)
    fig_interface(case, a.kr0, os.path.join(bundles, case, "profiles.dat"), a.outdir)
    fig_sweep(records, a.outdir)
    fig_linearity(records, a.outdir)
    if a.export:
        fig_field(a.export_case, a.export, os.path.join(bundles, a.export_case, "profiles.dat"), a.outdir, gif=not a.no_gif)

    print("figures from the older scripts")
    py = sys.executable
    run([py, os.path.join(HERE, "plot_campaign_results.py"), snap, "--outdir", a.outdir])
    run([py, os.path.join(HERE, "plot_pinchoff.py"), "--bundle", os.path.join(bundles, a.pinch), "--kr0", "0.7", "--cells-per-r0", a.pinch[-2:], "--outdir", a.outdir])
    run([py, os.path.join(HERE, "plot_ohnesorge.py"), "--outdir", a.outdir])
    run([py, os.path.join(HERE, "plot_schematic.py"), "--outdir", a.outdir])
    # Full-field figure (September's field_y_velocity_kr0_*): only where the bundle holds field slices, which
    # run_all.sh writes on the cluster but which are not in results/ (they are tens of MB).
    has_field = os.path.exists(os.path.join(bundles, case, "field_y_velocity.dat"))
    if has_field:
        run([py, os.path.join(HERE, "plot_results.py"), "--bundle", os.path.join(bundles, case), "--kr0", str(a.kr0), "--outdir", a.outdir])
    if not a.no_gif:
        for b in (case, a.pinch):
            field_opts = ["--field", "y_velocity"] if (b == case and has_field) else ["--no-field"]
            run([py, os.path.join(HERE, "animate.py"), "--bundle", os.path.join(bundles, b), "--outdir", a.outdir] + field_opts
                + (["--r-min", "0.0", "--r-max", "2.5", "--linear-bound", "100", "--stride", "2", "--dpi", "150", "--html-dpi", "100"] if b == a.pinch else []))
            for ext in ("gif", "html"):
                src = os.path.join(a.outdir, f"anim_interface.{ext}")
                if os.path.exists(src):
                    os.replace(src, os.path.join(a.outdir, f"anim_interface_{b}.{ext}"))
                src = os.path.join(a.outdir, f"anim_field_y_velocity.{ext}")
                if os.path.exists(src):
                    os.replace(src, os.path.join(a.outdir, f"anim_field_full_{b}.{ext}"))

    combine(a.outdir, ["schematic_twin_run", "campaign_dispersion", f"gallery_growth_{case}", f"gallery_local_rate_{case}",
                       f"gallery_interface_{case}", "gallery_sweep_growth", "gallery_linearity", f"gallery_field_{a.export_case}", f"field_y_velocity_kr0_{str(a.kr0).replace('.', 'p')}",
                       "campaign_stable_band", "pinchoff_sequence", "pinchoff_extrema", "pinchoff_extrapolation",
                       "ohnesorge_contour", "ohnesorge_surface", "ohnesorge_peak_drift"])


if __name__ == "__main__":
    main()
