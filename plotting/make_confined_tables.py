"""
make_confined_tables.py
=======================
Writes the confined two-fluid viscous theory tables used by plot_followup_ch6.py, after checking the implementation (confined_dispersion.py) against the tested Fortran tables.

    python make_confined_tables.py THEORY_DIR

THEORY_DIR must already hold theory_vs_k.csv and theory_vs_mu_ratio.csv from `run_lsa theory` (the unconfined two-fluid viscous relation of the Fortran code, campaign fluid properties). Checks, each of which stops the script on failure:

  1. without the wall, confined_dispersion reproduces every Fortran value in both tables to a relative 1e-10;
  2. with the wall, the rate at k r0 = 0.7 converges to the unconfined value as the wall moves out (|difference| < 1e-9 relative at Rw = 20).

Tables written (capillary units; Oh = 0.02, rho_2/rho_1 = 1.225e-3; mu_2 = 2.6e-4 unless stated):
  confined_vs_k_Rw1p5.csv   k r0 = 0.05 .. 0.995, wall at 1.5 R0
  confined_vs_Rw_k0p7.csv   k r0 = 0.7, wall at 1.2 .. 8.5 R0, for mu_2/mu_1 = 0.013 and 1
  confined_vs_mu_Rw4.csv    k r0 = 0.7, wall at 4 R0, mu_2/mu_1 = 0.01 .. 1.3
"""
import csv, os, sys
import numpy as np
from confined_dispersion import growth_rate

OH, RHO2, MU2 = 0.02, 1.225e-3, 2.6e-4

def main(tdir):
    tk = list(csv.DictReader(open(os.path.join(tdir, "theory_vs_k.csv"))))
    tm = list(csv.DictReader(open(os.path.join(tdir, "theory_vs_mu_ratio.csv"))))
    worst = 0.0
    for r in tk:
        k = float(r["kr0"]); f = float(r["two_fluid_viscous"])
        worst = max(worst, abs(growth_rate(k, OH, RHO2, MU2) / f - 1))
    for r in tm:
        m = float(r["mu_ratio"]); f = float(r["two_fluid_viscous"])
        worst = max(worst, abs(growth_rate(0.7, OH, RHO2, OH * m) / f - 1))
    print(f"check 1: largest relative difference from the Fortran tables, {len(tk) + len(tm)} values: {worst:.1e}")
    if worst > 1e-10:
        sys.exit("check 1 failed")
    su = growth_rate(0.7, OH, RHO2, MU2); s20 = growth_rate(0.7, OH, RHO2, MU2, 20.0)
    print(f"check 2: wall at 20 R0 differs from unconfined by {abs(s20 / su - 1):.1e}")
    if abs(s20 / su - 1) > 1e-9:
        sys.exit("check 2 failed")

    with open(os.path.join(tdir, "confined_vs_k_Rw1p5.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["kr0", "Rw", "unconfined", "confined"])
        for r in tk:
            k = float(r["kr0"]); w.writerow([f"{k:.4f}", 1.5, r["two_fluid_viscous"], f"{growth_rate(k, OH, RHO2, MU2, 1.5):.12e}"])
    with open(os.path.join(tdir, "confined_vs_Rw_k0p7.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["mu_ratio", "Rw", "unconfined", "confined"])
        for m in (0.013, 1.0):
            uc = growth_rate(0.7, OH, RHO2, OH * m)
            for Rw in np.round(np.arange(1.2, 8.501, 0.05), 3):
                w.writerow([m, f"{Rw:.3f}", f"{uc:.12e}", f"{growth_rate(0.7, OH, RHO2, OH * m, Rw):.12e}"])
    with open(os.path.join(tdir, "confined_vs_mu_Rw4.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["mu_ratio", "Rw", "unconfined", "confined"])
        for r in tm:
            m = float(r["mu_ratio"]); w.writerow([r["mu_ratio"], 4.0, r["two_fluid_viscous"], f"{growth_rate(0.7, OH, RHO2, OH * m, 4.0):.12e}"])
    print("tables written to", tdir)

if __name__ == "__main__":
    main(sys.argv[1])
