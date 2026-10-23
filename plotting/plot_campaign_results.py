#!/usr/bin/env python3
"""Figures for the October 2026 campaign, drawn from one results snapshot.

Reads one results snapshot (results/<name>/: analysis/summary.csv and records/) and draws:

    campaign_dispersion   measured growth rates of the 16-cell sweep (interface fit, narrow-band level-set DMD, velocity DMD) against Rayleigh's inviscid and Chandrasekhar's exact viscous relations, with the departure from the exact viscous rate below
    campaign_stable_band  the signed interface coefficient of stab110_16 (k r0 = 1.1) with its damped-cosine fit and the exact viscous prediction

Usage:
    python3 plotting/plot_campaign_results.py results/campaign_2026-10-09 --outdir figures
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

import dispersion  # noqa: E402
from pub_style import apply_style, cm_both, savefig_all  # noqa: E402

MU = 0.02          # liquid viscosity of the campaign; with rho = sigma = r0 = 1 this is also Oh
SWEEP = ["sw020_16", "sw030_16", "sw040_16", "sw050_16", "sw060_16",
         "sw070_16", "sw080_16", "sw090_16"]


def read_summary(path):
    """Rows of runs/analysis/summary.csv keyed by label, numbers as floats (NaN when empty)."""
    rows = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            row = {}
            for k, v in r.items():
                try:
                    row[k] = float(v)
                except (TypeError, ValueError):
                    row[k] = v if k == "label" else np.nan
            rows[r["label"]] = row
    return rows


def viscous_complex_root(x, oh, guess):
    """Complex root s of Chandrasekhar's relation for k r0 = x > 1 (a damped oscillation), by Newton iteration on mpmath's complex Bessel functions. Without mpmath, the same relation is solved with scipy's complex Bessel functions; the two agree to better than 1e-9."""
    try:
        import mpmath as mp
    except ImportError:
        from scipy.optimize import fsolve
        from scipy.special import iv

        def fs(s):
            l = np.sqrt(x * x + s / oh)
            qx, ql = iv(1, x) / iv(0, x), iv(1, l) / iv(0, l)
            br = (1 - qx / x) - (2 * x * l / (x * x + l * l)) * qx * (1 / ql - 1 / l)
            return s * s + 2 * oh * x * x * s * br - ((l * l - x * x) / (l * l + x * x)) * x * (1 - x * x) * qx

        re_, im_ = fsolve(lambda v: (lambda z: [z.real, z.imag])(fs(complex(v[0], v[1]))), list(guess), xtol=1e-12)
        return complex(re_, im_)

    x, oh = mp.mpf(x), mp.mpf(oh)

    def q(z):
        return mp.besseli(1, z) / mp.besseli(0, z)

    def f(s):
        l = mp.sqrt(x**2 + s / oh)
        qx, ql = q(x), q(l)
        br = (1 - qx / x) - (2 * x * l / (x**2 + l**2)) * qx * (1 / ql - 1 / l)
        return (s**2 + 2 * oh * x**2 * s * br
                - ((l**2 - x**2) / (l**2 + x**2)) * x * (1 - x**2) * qx)

    return complex(mp.findroot(f, mp.mpc(*guess)))


def figure_dispersion(rows, outdir):
    """Measured dispersion relation of the 16-cell sweep, with departures."""
    kr = np.array([rows[c]["kr0"] for c in SWEEP])
    visc = np.array([rows[c]["visc_exact"] for c in SWEEP])
    series = [("rate_interface", "interface fit", "o"),
              ("dmd_levelset_median", "level-set DMD (narrow band)", "s"),
              ("dmd_rankscan_median", "velocity DMD", "^")]

    x = np.linspace(0.02, 1.0, 400)
    inv = dispersion.growth_rate(x)
    vex = dispersion.growth_rate_viscous_exact(x, mu=MU)

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(5.6, 5.4), sharex=True,
                                 gridspec_kw={"height_ratios": [2.2, 1.0]})
    ax.plot(x, inv, "k--", lw=1.0, label="Rayleigh (inviscid)")
    ax.plot(x, vex, "k-", lw=1.0,
            label=rf"Chandrasekhar (viscous, Oh $= {MU:g}$)")
    for key, label, mk in series:
        y = np.array([rows[c][key] for c in SWEEP])
        ax.plot(kr, y, mk, ms=4.5, mfc="none", mew=0.9, label=label)
        bx.plot(kr, 100.0 * (y / visc - 1.0), mk + "-", ms=4.5, mfc="none",
                mew=0.9, lw=0.7)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.36)
    ax.set_ylabel(r"growth rate $\sigma\, t_c$")
    ax.legend(frameon=False, loc="lower center")
    cm_both(ax)

    bx.axhline(0.0, color="k", lw=0.7)
    bx.set_xlim(0.0, 1.0)
    bx.set_ylim(-3.0, 11.0)
    bx.set_xlabel(r"wavenumber $k r_0$")
    bx.set_ylabel(r"departure (%)")
    cm_both(bx)
    fig.tight_layout()
    savefig_all(fig, "campaign_dispersion", outdir=outdir, close=True)


def figure_stable_band(record, outdir):
    """Signed interface coefficient of the stable case with fit and theory (theory carries the fitted amplitude and phase)."""
    d = np.loadtxt(record, comments="#", ndmin=2)
    d = d[d[:, 1] > 0]
    t, ac, asn = d[:, 0], d[:, 3], d[:, 4]
    phi0 = np.arctan2(asn[0], ac[0])
    c = ac * np.cos(phi0) + asn * np.sin(phi0)

    def damped(tt, amp, delta, omega, phase):
        return amp * np.exp(-delta * tt) * np.cos(omega * tt + phase)

    p, _ = curve_fit(damped, t, c, p0=[c[0], 0.03, 0.33, 0.0], maxfev=20000)
    s = viscous_complex_root(1.1, MU, (-0.03, 0.33))

    tt = np.linspace(t.min(), t.max(), 600)
    fig, ax = plt.subplots(figsize=(5.6, 2.9))
    ax.plot(t[::3], 100 * c[::3], "o", ms=3.2, mfc="none", mew=0.8,
            label=r"measured, $k r_0 = 1.1$")
    ax.plot(tt, 100 * damped(tt, *p), "-", lw=1.0,
            label=rf"fit: $\omega = {abs(p[2]):.4f}$, $\delta = {p[1]:.4f}$")
    # Theory drawn with the fitted amplitude and phase, so that only the frequency and damping are compared.
    ax.plot(tt, 100 * damped(tt, p[0], -s.real, s.imag, p[3]), "k--",
            lw=0.9, label=rf"viscous theory: $\omega = {s.imag:.4f}$, "
                          rf"$\delta = {-s.real:.4f}$")
    ax.set_xlim(t.min(), t.max())
    ax.set_ylim(-1.15, 1.15)
    ax.set_xlabel(r"time $t / t_c$")
    ax.set_ylabel(r"$c(t) / r_0$ ($\times 100$)")
    ax.legend(frameon=False, loc="lower left")
    cm_both(ax)
    fig.tight_layout()
    savefig_all(fig, "campaign_stable_band", outdir=outdir, close=True)
    return p, s


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("snapshot", help="results snapshot directory")
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args()
    apply_style()
    rows = read_summary(os.path.join(args.snapshot, "analysis", "summary.csv"))
    figure_dispersion(rows, args.outdir)
    p, s = figure_stable_band(os.path.join(args.snapshot, "records",
                                           "stab110_16.dat"), args.outdir)
    print(f"stable band: omega {abs(p[2]):.6f} (theory {s.imag:.6f}), "
          f"delta {p[1]:.6f} (theory {-s.real:.6f})")
    print(f"figures written to {args.outdir}/")


if __name__ == "__main__":
    main()
