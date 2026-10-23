# Chapter 6 figures (22 October 2026)

Made by `plotting/plot_followup_ch6.py` from `results/followup_2026-10-22/analysis/summary.csv`; each in PNG (600 dpi), PDF and SVG, house style `pub_style.py`. Redraw with (from the top of the repository):

    mkdir -p figures && cp -r results/ch6_figures_2026-10-22/theory figures/
    cd plotting && python3 make_confined_tables.py ../figures/theory && cd ..
    python3 plotting/plot_followup_ch6.py results/followup_2026-10-22/analysis/summary.csv --theory figures/theory --outdir figures

- `ch6_confinement`: (a) growth rate against kR_0: sweep with the wall at 4 R_0 (circles) and the wall at 1.5 R_0 (squares; all at 16 cells, plus the single 32-cell check at kR_0 = 0.7, cross), with the two-fluid viscous relation unconfined (solid) and with the wall at 1.5 R_0 (dashed); (b) at kR_0 = 0.7, departure from the unconfined relation against the wall radius: measured interface fit and level-set DMD, and the confined relation as a curve.
- `ch6_viscous_ambient`: (a) growth rate at kR_0 = 0.7 against μ_2/μ_1: Chandrasekhar's relation (no gas), the two-fluid viscous relation unconfined (solid) and with the wall at 4 R_0 (dashed), the runs with the wall at 4 R_0 and, at μ_2/μ_1 = 1, at 6 and 8 R_0; (b) at μ_2/μ_1 = 1, departure from the unconfined relation against the wall radius, with the confined relation as a curve.

**Theory curves.** Both figures draw the exact two-fluid viscous relation (Tomotika 1935), the same relation used for every percentage. The curves are read from `theory/theory_vs_k.csv` and `theory/theory_vs_mu_ratio.csv`, written on the cluster on 22 October by the tested Fortran code with the campaign's fluid properties passed explicitly:

    run_lsa theory --kr0 K --r0 1 --sigma 1 --rho-l 1 --mu-l 0.02 --rho-g 1.225e-3 --mu-g MU_G --csv FILE

with K from 0.05 to 0.995 in steps of 0.005 at MU_G = 2.6e-4 (the campaign gas, μ_2/μ_1 = 0.013), and K = 0.7 at MU_G = 0.02 × (μ_2/μ_1) for 44 ratios from 0.01 to 1.3, including the four that were run. One value per call, because `run_lsa` accepts at most 64 wavenumbers per call. The tables agree with the `visc_two_fluid` column of `summary.csv` to within 5e-10 at every sampled point; the script checks this and stops if they disagree.

**Confined theory.** The runs have a no-slip wall at r = R_w (`ns.hi_bc = 5`). `plotting/confined_dispersion.py` solves the two-fluid viscous problem with that wall (6 × 6 determinant: both Bessel families in the ambient, no slip at R_w). `make_confined_tables.py` first checks it against the Fortran unconfined tables in `theory/` (largest relative difference 4.8e-13 over 234 values) and checks convergence to them as the wall moves out, then writes `theory/confined_vs_k_Rw1p5.csv`, `confined_vs_Rw_k0p7.csv` and `confined_vs_mu_Rw4.csv`. `confined_comparison.csv` (from `compare_confined.py`) compares every measured rate with both relations. The derivation is Appendix D of the dissertation, and the comparison is in Chapter 6.
