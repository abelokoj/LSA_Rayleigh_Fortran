# Python scripts: figures, tables and independent checks

The Fortran tools in `src/` do the analysis and write plain text. The Python scripts here draw the figures, recompute the tables of the dissertation from the stored results, and provide the independent checks of the theory (the wall-bounded dispersion relation is implemented here and checked against the Fortran library). They need `numpy`, `scipy` and `matplotlib` (`pip install -r plotting/requirements.txt`); none of them needs the solver or its plotfiles. `REPRODUCE.md` at the top of the repository gives the command for every table and figure.

## Dissertation tables and figures (`thesis/`)

| Script | What it does |
|---|---|
| `thesis/print_thesis_tables.py` | recomputes Tables 5.3–5.6 and 6.1–6.3 from `results/` and prints them in the layout of the dissertation |
| `thesis/make_extra_figures.py` | figures `ch5_window`, `ch5_proportionality`, `ch5_profiles`, `ch5_mesh`, `ch5_stride`, `ch6_weighting` |
| `thesis/make_dmd_evidence_figures.py` | figures `ch5_dmd_evidence` and `ch5_mode_shapes` |
| `thesis/compare_summaries.py` | compares a new `runs/analysis/summary.csv` with the stored one, case by case |
| `plot_campaign_results.py` | figures `campaign_dispersion` and `campaign_stable_band` |
| `plot_campaign_gallery.py` | the full figure set of the campaign, including `gallery_field_meas16` (the perturbation-field figure of Chapter 5) |
| `plot_followup_ch6.py` | figures `ch6_confinement` and `ch6_viscous_ambient` |

## Theory and checks

| Script | What it does |
|---|---|
| `dispersion.py` | inviscid, Chandrasekhar (exact viscous) and two-fluid inviscid relations; Python counterpart of `src/dispersion_mod.f90` |
| `confined_dispersion.py` | the two-fluid viscous relation with the no-slip outer wall (Appendix D of the dissertation), and its eigenfunction |
| `make_confined_tables.py` | checks `confined_dispersion.py` against the Fortran tables and writes the confined theory tables of Chapter 6 |
| `compare_confined.py` | every Chapter 6 rate against the unbounded and the wall-bounded relation (`results/ch6_figures_2026-10-22/confined_comparison.csv`) |
| `check_inviscid_limit.py` | the wall-bounded viscous relation tends to the closed-form inviscid confined relation as the viscosity vanishes |
| `check_mode_shapes_confined.py` | leading DMD velocity mode against the wall-bounded eigenfunction (rms departures quoted in Chapter 5) |
| `fit_stable_band.py` | damped-cosine fit of the stable mode at k R₀ = 1.1 (Table 5.5) |
| `weighting_test.py` | joint decomposition of velocity and level set against the weight λ (Table 6.3); needs the lossless exports, so it runs where the analysis job ran |
| `plot_mode_shape_archive.py` | reads a DMD mode archive (`studies/modes/`) and compares it with the unbounded eigenfunction |

## Supporting modules

`pub_style.py` (house style: Computer Modern, inward ticks, PNG/PDF/SVG at 600 dpi), `bundle.py` (the plain-text interface-profile format written by `bin/dump_profiles`), `animate.py` (interface and field animations), `plot_pinchoff.py`, `plot_ohnesorge.py`, `plot_schematic.py` (the remaining gallery figures), and `plot_results.py`, `fit_records.py`, `transient.py`, `nsmfp.py`, `dmd.py`, `interface_mode.py`, `report.py`, which are the Python counterparts of the Fortran tools used by the figure scripts above. `yt` is needed only by `plot_results.py --rundir`, which reads plotfiles directly.
