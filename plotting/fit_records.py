"""Fit a growth rate to each amplitude record and tabulate the results.

The sweep writes one record per case, holding time, interface mode amplitude and mean radius. This reads those records and reports, for each case, the fitted growth rate over a window, the coefficient of determination, the drift of the mean radius, and the departure from the inviscid and viscous predictions. It reads plain text only, so it runs on a login node or a workstation without the plotfiles and without yt.

The window matters more than anything else in this file. A fit over the whole record includes the startup transient and returns a substantial underestimate, so the default excludes the first half of the linear part of the record (transient.DEFAULT_KEEP, shared with every other script that fits a rate) and the value of t_min actually used is reported alongside every rate. Where a record is short, as in the amplitude sweep, the fraction retained is the same but the absolute window is shorter. When no sound window exists the fallback taken is printed beneath the table rather than applied silently.

Usage, from the working directory:

    python3 plotting/fit_records.py                         # every record
    python3 plotting/fit_records.py --frac 0.6              # discard 60 per cent of the linear part
    python3 plotting/fit_records.py --t-min 4.0             # an explicit window
    python3 plotting/fit_records.py runs/records/res08.dat  # named records only
"""

import argparse
import glob
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))

try:
    import dispersion
    HAVE_DISPERSION = True
except ImportError:
    HAVE_DISPERSION = False

try:
    import report
    HAVE_REPORT = True
except ImportError:
    HAVE_REPORT = False

# transient imports SciPy only inside the transient fit, so the window rule shared by every script is always available; only the transient-corrected rate below needs SciPy, and assess() reports rather than raises when the fit cannot be made.
import transient  # noqa: E402

#: Columns written by bin/dump_amplitude: time, interface mode amplitude, mean radius.
N_COLUMNS = 3


def kr0_of(label, override=None):
    """Recover the wavenumber a case was run at from its label.

    The sweep labels encode the wavenumber only for the dispersion cases; every other case is at the reference value, which is what the convention in full_analysis.txt fixes. An explicit override is accepted for records whose label does not follow that convention, since inferring a wavenumber from a name and being silently wrong is worse than being told.
    """
    if override is not None:
        return override
    if label.startswith("disp"):
        digits = "".join(c for c in label[4:] if c.isdigit())
        if digits:
            return int(digits) / 100.0
    return 0.7


def fit(t, a, t_min):
    """Log-linear fit of a(t) over t >= t_min, returning the rate and r squared."""
    m = (t >= t_min) & (a > 0.0) & np.isfinite(a)
    if m.sum() < 3:
        return float("nan"), float("nan"), int(m.sum())
    x, y = t[m], np.log(a[m])
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), r2, int(m.sum())


def main(argv=None):
    ap = argparse.ArgumentParser(description="growth rates from amplitude records")
    ap.add_argument("records", nargs="*",
                    help="record files; defaults to runs/records/*.dat")
    ap.add_argument("--frac", type=float, default=transient.DEFAULT_KEEP,
                    help="fraction of the linear part of each record to "
                         "discard as transient (default %(default)s, shared "
                         "with the plotting scripts)")
    ap.add_argument("--t-min", type=float, default=None,
                    help="explicit start of the fit window, overriding --frac")
    ap.add_argument("--r0", type=float, default=1.0,
                    help="column radius; NOT the wavenumber, see --kr0")
    ap.add_argument("--kr0", type=float, default=None,
                    help="wavenumber k*r0, overriding the value inferred "
                         "from each record's label")
    ap.add_argument("--linear-bound", type=float, default=None,
                    dest="linear_bound",
                    help="amplitude, as a fraction of r0, beyond which the "
                         "record is treated as nonlinear and excluded")
    ap.add_argument("--surface-tension", type=float, default=1.0)
    ap.add_argument("--rho", type=float, default=1.0)
    ap.add_argument("--mu", type=float, default=2.0e-2,
                    help="dynamic viscosity of the liquid, for the exact "
                         "viscous prediction")
    args = ap.parse_args(argv)

    files = args.records or sorted(glob.glob("runs/records/*.dat"))
    if not files:
        print("no records found under runs/records/", file=sys.stderr)
        return 1

    bound = (args.linear_bound if args.linear_bound is not None
             else transient.LINEAR_BOUND)

    rows = []
    warnings = []
    for path in files:
        label = os.path.basename(path)[:-4]
        # Only the failures that a damaged or unfinished record actually produces are caught: a missing or unreadable file (OSError) and text that is not a numeric table (ValueError). Anything else is a bug in this script and should surface as one rather than as a row saying "unreadable".
        try:
            d = np.loadtxt(path, comments="#", ndmin=2)
        except (OSError, ValueError) as exc:
            rows.append([label, "-", "unreadable", "-", "-", "-", "-", "-",
                         "see below"])
            warnings.append(f"{label}: could not be read as a numeric table "
                            f"({exc})")
            continue
        if d.shape[1] < N_COLUMNS:
            rows.append([label, "-", "bad columns", "-", "-", "-", "-", "-",
                         "see below"])
            warnings.append(f"{label}: expected {N_COLUMNS} columns (time, "
                            f"amplitude, mean radius) as written by "
                            f"bin/dump_amplitude, found {d.shape[1]}")
            continue
        if d.shape[0] < 3:
            rows.append([label, "-", f"{d.shape[0]} rows", "-", "-", "-",
                         "-", "-", "too short"])
            continue

        t, a, r = d[:, 0], d[:, 1], d[:, 2]
        kr0 = kr0_of(label, args.kr0)

        # The window is bounded on both sides. Its upper limit is where the amplitude leaves the linear range, because a fit that straddles break-up describes neither regime; its lower limit is placed inside the linear part rather than as a fraction of the whole record, since a run carried well past linearity has a large fraction of its record in the wreckage. An explicit --t-min still overrides the lower bound. The rule and its fallbacks live in transient.choose_window so that every script applies the same one.
        t_min, t_hi, nonlinear, window_note = transient.choose_window(
            t, a, r0=args.r0, bound=bound, keep=args.frac, t_min=args.t_min)
        if window_note:
            warnings.append(f"{label}: {window_note}")

        window = (t >= t_min) & (t <= t_hi)
        sigma, r2, n = fit(t[window], a[window], t_min)

        # Corrected for the transient still present at the end of the linear window. Where the amplitude leaves the linear range before the transient decays, no window is clean and this extrapolation is the only defensible estimate; it is reported alongside the direct fit rather than in place of it.
        verdict = transient.assess(t, a, r0=args.r0, bound=bound)
        sigma_corrected = verdict.get("sigma_inf", float("nan"))

        # bin/dump_amplitude writes -1 where no interface was found, so negative and non-finite mean radii are missing values, not radii, and are excluded before the drift is measured.
        r_ok = r[np.isfinite(r) & (r > 0.0)]
        drift = ((r_ok.max() - r_ok.min()) / r_ok.mean() if r_ok.size >= 2
                 else float("nan"))

        stable = False
        period = float("nan")
        vis = float("nan")
        if HAVE_DISPERSION:
            # dispersion takes the dimensional wavenumber k and forms k*r0 itself, so kr0 is divided by r0 here, as run_lsa._theory does. Passing kr0 directly applied r0 twice and was wrong for any r0 other than 1.
            k = kr0 / args.r0
            inv = float(dispersion.growth_rate(k, r0=args.r0,
                                               sigma=args.surface_tension,
                                               rho=args.rho))
            # The exact viscous root is the appropriate reference at the finite Ohnesorge number of the decks; the inviscid rate is only an upper bound. The quadratic approximation is used only if an older dispersion module lacks the exact solver.
            viscous_fn = getattr(dispersion, "growth_rate_viscous_exact",
                                 dispersion.growth_rate_viscous)
            vis = float(viscous_fn(k, r0=args.r0, sigma=args.surface_tension,
                                   rho=args.rho, mu=args.mu))
            # Above kr0 = 1 the squared growth rate is negative, so the mode does not decay exponentially but oscillates. A log-linear fit of an oscillation is meaningless, and reporting one invites the reader to compare a number against a theory that does not predict it, so the period is reported instead.
            s2q = float(dispersion.growth_rate_squared(
                k, r0=args.r0, sigma=args.surface_tension, rho=args.rho))
            if s2q < 0.0:
                stable = True
                omega = math.sqrt(-s2q)
                period = 2.0 * math.pi / omega if omega > 0 else float("nan")
        else:
            inv = float("nan")

        dev = (sigma - inv) / inv * 100.0 if inv and math.isfinite(inv) \
            and inv != 0 else float("nan")
        dev_vis = (sigma - vis) / vis * 100.0 if math.isfinite(vis) \
            and vis > 0 else float("nan")

        note = ""
        if stable:
            span = (t[-1] - t[0]) / period if period > 0 else float("nan")
            ratio = a[-1] / a[0] if a[0] != 0 else float("nan")
            note = (f"stable: oscillates, T={period:.2f}, "
                    f"record spans {span:.2f} T, a_end/a_0={ratio:.3f}")
            rows.append([label, f"{kr0:.2f}", "oscillatory", "-", "0 (stable)",
                         "-", "0 (stable)", "-", note])
            continue
        if label.endswith(".base"):
            # The unperturbed twin carries no imposed mode, so its amplitude is discretisation noise and a fitted rate is meaningless. It is reported because the amplitude itself is the useful quantity: it bounds the noise floor against which every perturbed case is measured.
            note = f"twin: noise floor a={a.max():.2e}"
        elif not math.isfinite(sigma):
            note = "no fit"
        elif math.isfinite(inv) and inv > 0 and sigma > inv * 1.02:
            note = "ABOVE INVISCID: check the deck"
        elif nonlinear and math.isfinite(sigma_corrected):
            # A record that reached the nonlinear range is not necessarily unusable, but its direct fit is biased low, so the corrected value is what should be quoted.
            note = (f"nonlinear past t={t_hi:.1f}; "
                    f"transient-corrected sigma={sigma_corrected:.4f}")
        elif math.isfinite(r2) and r2 < 0.99:
            note = "poor fit: widen the window"
        elif abs(drift) > 0.01:
            note = "mass drift > 1 per cent"

        rows.append([label, f"{kr0:.2f}", f"{sigma:.6f}", f"{r2:.5f}",
                     f"{inv:.6f}" if math.isfinite(inv) else "-",
                     f"{dev:+.2f} %" if math.isfinite(dev) else "-",
                     f"{vis:.6f}" if math.isfinite(vis) else "-",
                     f"{dev_vis:+.2f} %" if math.isfinite(dev_vis) else "-",
                     note or f"{t_min:.2f}<t<{t_hi:.2f}, n={n}"])

    headers = ["case", "k r0", "sigma", "r^2", "inviscid", "vs inviscid",
               f"viscous (mu={args.mu:g})", "vs viscous", "note"]
    if HAVE_REPORT:
        print(report.table(rows, headers,
                           colalign=["left", "right", "right", "right",
                                     "right", "right", "right", "right",
                                     "left"],
                           title="Growth rates fitted to the amplitude records"))
    else:
        print("  ".join(headers))
        for row in rows:
            print("  ".join(str(c) for c in row))
    # Printed after the table rather than inside it, so that the reasons are readable in full; a rate whose case appears here came from a fallback window or was not fitted at all.
    if warnings:
        print("\nwarnings:", file=sys.stderr)
        for line in warnings:
            print(f"  {line}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
