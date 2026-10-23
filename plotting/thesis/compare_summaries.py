#!/usr/bin/env python3
"""
compare_summaries.py
====================
Compare the growth rates of a new analysis (runs/analysis/summary.csv, written by the analysis job after a rerun) with the stored ones in results/, case by case.

    python3 plotting/thesis/compare_summaries.py runs/analysis/summary.csv
    python3 plotting/thesis/compare_summaries.py NEW.csv results/followup_2026-10-22/analysis/summary.csv

For each case present in both files it prints the interface rate and the two decomposition rates of each, and their relative difference. Cases present in only one file are listed at the end. A rerun on the same solver commit, rank count and deck reproduces the stored interface rates to round-off; a different rank count or compiler changes the round-off and moves the decomposition rates by a small fraction of a percent, which is why the difference is printed rather than tested for equality.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "..", "..", "results", "followup_2026-10-22", "analysis", "summary.csv")
COLS = (("rate_interface", "interface"), ("dmd_levelset_median", "level set"), ("dmd_rankscan_median", "velocity"))


def load(path):
    rows = {}
    for r in csv.DictReader(open(path)):
        r = {k.strip(): (v or "").strip() for k, v in r.items()}
        rows[r["label"]] = r
    return rows


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    new, ref = load(sys.argv[1]), load(sys.argv[2] if len(sys.argv) > 2 else REF)
    common = [l for l in ref if l in new]
    print(f"{'case':16s}" + "".join(f"{name:>34s}" for _, name in COLS))
    print(f"{'':16s}" + "".join(f"{'stored':>12s}{'new':>12s}{'diff %':>10s}" for _ in COLS))
    worst = 0.0
    for lab in common:
        cells = []
        for key, _ in COLS:
            a, b = ref[lab].get(key, ""), new[lab].get(key, "")
            if a and b:
                d = 100.0 * (float(b) / float(a) - 1.0)
                if key == "rate_interface":
                    worst = max(worst, abs(d))
                cells.append(f"{float(a):12.6f}{float(b):12.6f}{d:+10.3f}")
            else:
                cells.append(f"{a[:10] or '-':>12s}{b[:10] or '-':>12s}{'':>10s}")
        print(f"{lab:16s}" + "".join(cells))
    print(f"\n{len(common)} cases compared; largest difference of the interface rate: {worst:.3f} %")
    only_ref = [l for l in ref if l not in new]
    only_new = [l for l in new if l not in ref]
    if only_ref:
        print("not in the new analysis: " + " ".join(only_ref))
    if only_new:
        print("only in the new analysis: " + " ".join(only_new))


if __name__ == "__main__":
    main()
