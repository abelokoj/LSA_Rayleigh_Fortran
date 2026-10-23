#!/usr/bin/env python3
"""
print_thesis_tables.py
======================
Recompute the numerical tables of Chapters 5 and 6 of the dissertation from the data in results/, and print them in the layout of the dissertation, so that each entry can be checked against the printed table.

    python3 plotting/thesis/print_thesis_tables.py                       # the stored results
    python3 plotting/thesis/print_thesis_tables.py --snapshot runs \
        --weighting runs/analysis/weighting_test_output.txt              # a new set of runs (REPRODUCE.md, Level 3)

Tables reproduced:

  5.3  growth rates of the 16-cell sweep (interface fit, level-set and velocity decompositions) against Chandrasekhar's relation
  5.4  mesh study at x = 0.7, eps = 5e-3, on 8, 16 and 32 cells per radius
  5.5  frequency and damping of the stable mode at x = 1.1
  5.6  velocity-decomposition rate against the sampling stride
  6.1  interface rates of the sweep against the single-fluid, two-fluid and wall-bounded two-fluid relations
  6.2  proportionality to the amplitude on 8, 16 and 32 cells per radius
  6.3  joint velocity and level-set decomposition against the weight lambda

Inputs (all in results/ by default): analysis/summary.csv, analysis/records/*.txt and records/stab110_16.dat of the follow-up snapshot (which contains every case of the October 2026 campaign as well as the follow-up runs), analysis/studies/stride_study.csv, the linearity reports linearity_08/16/32.txt, and the raw output of the weighting test. The wall-bounded relation is plotting/confined_dispersion.py (Appendix D of the dissertation). Needs numpy and scipy only. Values are printed with one more digit than the dissertation where rounding could otherwise hide a difference.
"""
import csv
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "plotting"))
import confined_dispersion as cd  # noqa: E402

SNAP = os.path.join(ROOT, "results", "followup_2026-10-22")
WEIGHTING = os.path.join(ROOT, "results", "estimators_2026-10-22", "weighting_test_output_2026-10-22.txt")
OH, RHO2, MU2, RW = 0.02, 1.225e-3, 2.6e-4, 4.0  # campaign fluids: Oh, rho2/rho1, mu2 (mu1 = Oh), wall radius

SWEEP = ["sw020_16", "sw030_16", "sw040_16", "sw050_16", "sw060_16", "sw070_16", "sw080_16", "sw090_16"]


def summary():
    rows = {}
    for r in csv.DictReader(open(os.path.join(SNAP, "analysis", "summary.csv"))):
        r = {k.strip(): (v or "").strip() for k, v in r.items()}
        rows[r["label"]] = r
    return rows


def num(r, key):
    return float(r[key]) if r.get(key) else float("nan")


def pct(a, b):
    return 100.0 * (a / b - 1.0)


def spread(label):
    """Start-time spread of the interface rate, in percent, from the record report written by run_lsa record."""
    txt = open(os.path.join(SNAP, "analysis", "records", f"{label}.txt")).read()
    m = re.search(r"growth rate s =\s+([\d.]+).*?spread\s+(\S+)\s+\(", txt, re.S)
    return 100.0 * float(m.group(2)) / float(m.group(1))


def line(cells, widths):
    return "  ".join(str(c).rjust(w) for c, w in zip(cells, widths))


def title(s):
    print("\n" + s + "\n" + "-" * len(s))


def table_5_3(R):
    title("Table 5.3  Growth rates of the 16-cell sweep, eps = 2.5e-3 (departures from Chandrasekhar's relation, %)")
    w = (5, 9, 9, 8, 7, 9, 9)
    print(line(["x", "Exact", "Interf.", "e-folds", "Spread", "Lev.set", "Velocity"], w))
    for lab in SWEEP:
        r = R[lab]; ex = num(r, "visc_exact")
        v = pct(num(r, "dmd_rankscan_median"), ex)
        vs = f"{v:+.2f}" if lab not in ("sw080_16", "sw090_16") else f"({v:+.2f})"
        print(line([f"{num(r, 'kr0'):.1f}", f"{ex:.6f}", f"{pct(num(r, 'rate_interface'), ex):+.2f}", f"{num(r, 'efolds'):.2f}",
                    f"{spread(lab):.3f}", f"{pct(num(r, 'dmd_levelset_median'), ex):+.2f}", vs], w))
    print("Velocity: median over truncation ranks (dmd_rankscan_median); in parentheses where the decomposition is not a measurement (Sec. 5.7).")


def table_5_4(R):
    title("Table 5.4  Mesh study at x = 0.7, eps = 5e-3 (departures from Chandrasekhar's relation, %)")
    w = (6, 9, 8, 7, 9, 9, 10)
    print(line(["Cells", "Interf.", "e-folds", "Spread", "Lev.set", "Velocity", "Vel. IQR"], w))
    for n, lab in ((8, "meas08"), (16, "meas16"), (32, "meas32")):
        r = R[lab]; ex = num(r, "visc_exact"); vel = num(r, "dmd_rankscan_median")
        print(line([n, f"{pct(num(r, 'rate_interface'), ex):+.2f}", f"{num(r, 'efolds'):.2f}", f"{spread(lab):.3f}",
                    f"{pct(num(r, 'dmd_levelset_median'), ex):+.2f}", f"{pct(vel, ex):+.2f}",
                    f"{100 * num(r, 'dmd_rankscan_iqr') / vel:.1f} %"], w))
    print("Vel. IQR: interquartile range of the velocity rate over truncation ranks, relative to the rate (flags the 16-cell value).")


def table_5_5():
    title("Table 5.5  Stable mode at x = 1.1 (damped-cosine fit to the signed interface coefficient, plotting/fit_stable_band.py)")
    rec = os.path.relpath(os.path.join(SNAP, "records", "stab110_16.dat"), ROOT)
    out = subprocess.run([sys.executable, os.path.join("plotting", "fit_stable_band.py"), rec, "--kr0", "1.1", "--t-min", "0"],
                         capture_output=True, text=True, check=True, cwd=ROOT).stdout
    print(out.rstrip())
    om = float(re.search(r"omega\s+([\d.]+)", out).group(1)); dl = float(re.search(r"delta\s+([\d.]+)", out).group(1))
    # exact complex root of Chandrasekhar's relation, s = -delta + i omega (plotting/dispersion.py)
    import dispersion
    from scipy.optimize import fsolve
    res = lambda v: (lambda z: [z.real, z.imag])(dispersion._viscous_residual(complex(v[0], v[1]), 1.1, OH))
    sr, si = fsolve(res, [-0.03, 0.33], xtol=1e-14)
    print(f"Exact root of Chandrasekhar's relation at x = 1.1, Oh = {OH}: s = {sr:.6f} + {si:.6f} i")
    print(f"Differences: omega {pct(om, si):+.2f} %, delta {pct(dl, -sr):+.1f} %")


def table_5_6():
    title("Table 5.6  Velocity-decomposition rate against the sampling stride (departure from Chandrasekhar's relation, %)")
    rows = list(csv.DictReader(open(os.path.join(SNAP, "analysis", "studies", "stride_study.csv"))))
    w = (5, 8, 9, 9, 9)
    print(line(["x", "start", "stride 1", "stride 2", "stride 4"], w))
    for lab in SWEEP[:6]:
        dep = {int(r["stride"]): float(r["visc_exact_dep_pct_median"]) for r in rows if r["label"] == lab and r["kind"] == "vel"}
        st = next(csv.DictReader(open(os.path.join(SNAP, "analysis", lab, "studies", "stride_vel.csv"))))["start"].strip()
        print(line([lab[2] + "." + lab[3], st, f"{dep[1]:+.2f}", f"{dep[2]:+.2f}", f"{dep[4]:+.2f}"], w))
    print("Median over truncation ranks 4, 8, 16 and full rank; window start in stored snapshots.")


def table_6_1(R):
    title("Table 6.1  Interface rates of the sweep against three relations (departure, %)")
    w = (5, 8, 13, 11, 16)
    print(line(["x", "e-folds", "single-fluid", "two-fluid", "two-fluid, wall"], w))
    for lab in SWEEP:
        r = R[lab]; k = num(r, "kr0"); s = num(r, "rate_interface")
        wall = cd.growth_rate(k, OH, RHO2, MU2, RW)
        print(line([f"{k:.1f}", f"{num(r, 'efolds'):.2f}", f"{pct(s, num(r, 'visc_exact')):+.2f}",
                    f"{pct(s, num(r, 'visc_two_fluid')):+.2f}", f"{pct(s, wall):+.2f}"], w))
    print("Two-fluid, wall: the relation with the no-slip wall at R_w = 4 R_0 (plotting/confined_dispersion.py).")


def table_6_2(R):
    title("Table 6.2  Proportionality at x = 0.7 on three grids (interface rate vs unbounded two-fluid relation, %, e-folds in brackets)")
    sets = {8: ("amp08_025", "amp08_050", "amp08_100"), 16: ("sw070_16", "lin050_16", "lin100_16"), 32: ("amp32_025", "amp32_050", "amp32_100")}
    w = (6, 16, 16, 16, 12)
    print(line(["Cells", "eps = 2.5e-3", "eps = 5e-3", "eps = 1e-2", "norm spread"], w))
    for n, labs in sets.items():
        cells = []
        for lab in labs:
            r = R[lab]
            cells.append(f"{pct(num(r, 'rate_interface'), num(r, 'visc_two_fluid')):+.2f} ({num(r, 'efolds'):.2f})")
        txt = open(os.path.join(SNAP, "analysis", f"linearity_{n:02d}.txt")).read()
        sp = re.search(r"scaled-norm spread\s+([\d.]+) %", txt).group(1)
        print(line([n] + cells + [f"{float(sp):.2f} %"], w))


def table_6_3():
    title("Table 6.3  Joint decomposition against the weight lambda (departure from the two-fluid relation with the wall, %)")
    th = {"sw040_16": (0.4, RW, MU2), "sw050_16": (0.5, RW, MU2), "sw060_16": (0.6, RW, MU2), "sw070_16": (0.7, RW, MU2),
          "meas32": (0.7, RW, MU2), "conf15_32": (0.7, 1.5, MU2), "visc100r8_16": (0.7, 8.0, OH),
          "sw080_16": (0.8, RW, MU2), "meas16": (0.7, RW, MU2)}
    txt = open(WEIGHTING).read()
    w = (13, 6, 8, 7, 7, 7, 7, 7, 8)
    print(line(["case", "ratio", "vel.only", "0.01", "0.1", "1", "10", "100", "l.s.only"], w))
    for lab, (k, rw, mu2) in th.items():
        blk = txt.split(f"== {lab}:")[1].split("==")[0]
        ratio = re.search(r"at lambda = 1: ([\d.]+)", blk).group(1)
        rates = [float(v) for v in re.findall(r"rate ([\d.]+)", blk)]
        s = cd.growth_rate(k, OH, RHO2, mu2, rw)
        print(line([lab, ratio] + [f"{pct(x, s):+.2f}" for x in rates], w))
    print("sw080_16 and meas16: velocity decomposition flagged. Rates are read from the raw output of plotting/weighting_test.py.")


def main():
    global SNAP, WEIGHTING
    import argparse
    ap = argparse.ArgumentParser(description="Recompute the tables of Chapters 5 and 6 of the dissertation.")
    ap.add_argument("--snapshot", default=SNAP, help="folder holding analysis/ and records/ (default: the stored follow-up snapshot; use runs for a new set of runs)")
    ap.add_argument("--weighting", default=WEIGHTING, help="raw output of plotting/weighting_test.py (default: the stored output)")
    a = ap.parse_args()
    SNAP, WEIGHTING = os.path.abspath(a.snapshot), os.path.abspath(a.weighting)
    print(f"snapshot: {os.path.relpath(SNAP, ROOT)}    weighting output: {os.path.relpath(WEIGHTING, ROOT)}")
    R = summary()
    table_5_3(R); table_5_4(R); table_5_5(); table_5_6(); table_6_1(R); table_6_2(R); table_6_3()
    print()


if __name__ == "__main__":
    main()
