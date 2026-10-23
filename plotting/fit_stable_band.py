#!/usr/bin/env python3
"""Frequency and damping of a stable-band case (k r0 > 1) from the signed interface coefficient.

The September 2026 evidence report pointed out that the magnitude record a(t) of a stable case cannot be fitted with a single exponential: the magnitude folds a damped oscillation, so its logarithm mixes the oscillation phase with the damping. Since 2026-10-07, bin/dump_amplitude also writes the signed coefficients of cos(kz) and sin(kz) (columns 4 and 5). This script projects them onto the direction of the initial perturbation, c(t) = a_cos cos(phi0) + a_sin sin(phi0), and fits

    c(t) = A exp(-delta t) cos(omega t + phi)

by nonlinear least squares, reporting omega and delta beside the inviscid frequency sqrt(-s^2) of the Rayleigh relation. The record must cover a good part of a period for omega to be meaningful; the script says how much of one it covers.

    python plotting/fit_stable_band.py runs/records/stab110_16.dat --kr0 1.1 [--t-min 0.5] [--plot fig.pdf]
"""

import argparse
import sys

import numpy as np
from scipy.optimize import curve_fit

import dispersion


def damped(t, amp, delta, omega, phase):
    return amp * np.exp(-delta * t) * np.cos(omega * t + phase)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("record")
    ap.add_argument("--kr0", type=float, required=True)
    ap.add_argument("--r0", type=float, default=1.0)
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--rho", type=float, default=1.0)
    ap.add_argument("--t-min", type=float, default=0.0, help="discard samples before this time (startup)")
    ap.add_argument("--plot", default=None, help="write a figure of the fit to this path")
    args = ap.parse_args(argv)

    d = np.loadtxt(args.record, comments="#", ndmin=2)
    if d.shape[1] < 5:
        sys.exit(f"{args.record}: {d.shape[1]} columns; the signed coefficients (columns 4-5) are written by bin/dump_amplitude since 2026-10-07, so this record must be regenerated from its plotfiles")
    ok = d[:, 1] > 0
    t, ac, as_ = d[ok, 0], d[ok, 3], d[ok, 4]
    phi0 = np.arctan2(as_[0], ac[0])
    c = ac * np.cos(phi0) + as_ * np.sin(phi0)
    keep = t >= args.t_min
    t, c = t[keep], c[keep]

    s2 = float(dispersion.growth_rate_squared(args.kr0 / args.r0, r0=args.r0, sigma=args.sigma, rho=args.rho))
    omega_inv = np.sqrt(-s2) if s2 < 0 else np.nan
    guess = [c[0], 0.1, omega_inv if np.isfinite(omega_inv) else 1.0, 0.0]
    popt, pcov = curve_fit(damped, t, c, p0=guess, maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    amp, delta, omega, phase = popt
    resid = c - damped(t, *popt)
    rel = np.sqrt(np.mean(resid**2)) / max(np.max(np.abs(c)), 1e-300)
    period = 2.0 * np.pi / abs(omega) if omega != 0 else np.inf

    print(f"record           {args.record}")
    print(f"samples          {t.size}, t = {t[0]:.4f} .. {t[-1]:.4f}")
    print(f"omega            {abs(omega):.6f} +- {perr[2]:.2e}   (inviscid {omega_inv:.6f})")
    print(f"damping delta    {delta:.6f} +- {perr[1]:.2e}")
    print(f"record covers    {(t[-1] - t[0]) / period:.2f} of a period")
    print(f"rms residual     {100 * rel:.3f} % of the peak")
    if (t[-1] - t[0]) < 0.5 * period:
        print("WARNING: less than half a period is covered; omega and delta are strongly correlated and should not be quoted separately")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from pub_style import apply_style
        apply_style()
        fig, ax = plt.subplots(figsize=(3.4, 2.4))
        ax.plot(t, c, ".", ms=2.5, label="signed coefficient")
        tf = np.linspace(t[0], t[-1], 400)
        ax.plot(tf, damped(tf, *popt), "-", lw=1.0, label="damped cosine fit")
        ax.set_xlabel(r"$t/t_c$")
        ax.set_ylabel(r"$c(t)/r_0$")
        ax.legend(frameon=False)
        fig.savefig(args.plot)
        print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
