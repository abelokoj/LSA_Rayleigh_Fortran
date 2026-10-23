# LSA_Rayleigh_Fortran

This is the analysis code and data behind the dissertation *Global Linear Stability Analysis of Capillary-Driven Interfacial Instability: A Jacobian-Free Approach* (Abel Okojunu, Department of Mathematics, Florida State University; advisor Prof. Mark Sussman). It is provided so that the dissertation committee can inspect the code, run its tests, and reproduce every reported number from the data it came from.

**To reproduce the results, follow [`REPRODUCE.md`](REPRODUCE.md).** It has three levels: check the code (5 minutes), recompute every table and figure from the stored data (15 minutes, no solver needed), and rerun the simulations on a cluster.

## The method in one paragraph

The growth rate of the Rayleigh–Plateau instability of a liquid column is measured with an unmodified two-phase Navier–Stokes solver: the coupled level set and continuous moment-of-fluid (CLS-CMOF) method of Ye et al. (2024), in `amrex_implicit_interfaces`. Two runs are made from the same initial state, one with a small sinusoidal displacement of the interface and one without (the *twin*). Their difference, sampled in time, evolves under the linearized operator of the discrete solver (the NS-MFP idea of Ranjan, Unnikrishnan and Gaitonde, 2020). The growth rate is extracted in three independent ways:

- a least-squares fit of the interface mode amplitude;
- dynamic mode decomposition (DMD) of the velocity difference;
- DMD of the level-set difference in a narrow band around the interface.

The results are compared with exact linear theory: Chandrasekhar's viscous relation, Tomotika's two-fluid viscous relation, and that relation with the no-slip outer wall of the computational domain, derived in Appendix D of the dissertation.

## Quick check (any machine with gfortran and LAPACK)

```bash
make            # builds bin/run_lsa, bin/make_inputs, bin/dump_amplitude, bin/dump_profiles, bin/svd_and_diagonalize
make test       # the automated test suite; must end with ALL TEST PROGRAMS PASSED
./bin/run_lsa theory --kr0 0.2,0.4,0.6,0.7,0.8,0.9 --r-wall 4      # every theoretical relation used in the dissertation
```

`run_lsa theory` prints the inviscid, Weber, Chandrasekhar, two-fluid inviscid, two-fluid viscous (Tomotika) and, with `--r-wall`, wall-bounded two-fluid viscous growth rates for the campaign fluids (Oh = 0.02, ρ₂/ρ₁ = 1.225 × 10⁻³, μ₂/μ₁ = 0.013).

## Where each result comes from

| Dissertation result | Data in this repository | Reproduced by (`REPRODUCE.md`) |
|---|---|---|
| Dispersion relation, k R₀ = 0.2–0.95 (Table 5.3, Ch. 5) | `results/followup_2026-10-22/analysis/summary.csv`, `records/` | `scripts/refit_records.sh`, `plotting/thesis/print_thesis_tables.py`, `plotting/plot_campaign_results.py` |
| Mesh study, 8, 16, 32 cells per R₀ (Table 5.4) | `summary.csv` rows `meas08/16/32` | `print_thesis_tables.py`, `plotting/thesis/make_extra_figures.py` |
| Stable mode, k R₀ = 1.1 (Table 5.5) | `records/stab110_16.dat` | `print_thesis_tables.py` (through `plotting/fit_stable_band.py`), `plot_campaign_results.py` |
| Sampling stride (Table 5.6) | `analysis/studies/stride_study.csv` | `print_thesis_tables.py` |
| Decomposition evidence and mode shapes (Ch. 5) | `analysis/sw070_16/` (`dmd`, `dmd_ls`, `converge.csv`), `analysis/<case>/studies/modes/` | `plotting/thesis/make_dmd_evidence_figures.py`, `plotting/check_mode_shapes_confined.py` |
| Comparison with three relations (Table 6.1) | `summary.csv` | `print_thesis_tables.py` |
| Linearity at 8, 16, 32 cells per R₀ (Table 6.2) | `analysis/linearity_08.txt`, `_16.txt`, `_32.txt`, `summary.csv` | `print_thesis_tables.py` |
| Confined column and viscous ambient (Ch. 6) | `results/ch6_figures_2026-10-22/` (figures, `confined_comparison.csv`, theory tables) | `plotting/compare_confined.py`, `plotting/plot_followup_ch6.py`, `plotting/make_confined_tables.py` |
| Weighting of velocity and level set (Table 6.3) | `results/estimators_2026-10-22/weighting_test_output_2026-10-22.txt` | `print_thesis_tables.py`; the test itself is `plotting/weighting_test.py` |
| Choice of estimator and the quoting rule (Ch. 4) | `results/estimators_2026-10-22/table.md`, `all_estimators.csv` | — |
| Base-state drift of the twin (Ch. 5) | `analysis/studies/drift_sw070_16.csv` | `scripts/run_studies.sh drift` |
| Restarting (Krylov) variant at k R₀ = 0.7 (Ch. 5) | `results/restarting_members_2026-10-23/` | `krylov_test/` |
| Figure gallery and break-up | `results/gallery_2026-10-16/` | `plotting/plot_campaign_gallery.py` |
| The September 2026 runs, discussed as history | `archived_records/` | `run_lsa record` |

Every run has a run record, `results/*/run_records/<case>/manifest.txt`: the solver commit `dddb827c` (whose solver source is identical to commit `0b8df42` of Prof. Sussman's repository), the executable's MD5 checksum, the rank count, the case line and the modules.

## Layout

| Path | Contents |
|---|---|
| `REPRODUCE.md` | the commands that reproduce every result |
| `src/`, `tests/`, `Makefile` | the Fortran 2008 analysis code and its tests |
| `plotting/` | Python scripts for the figures, the tables (`plotting/thesis/`) and the independent checks of the theory |
| `scripts/` | Slurm and shell scripts: running cases, the analysis job, the studies, comparison with the stored results |
| `inputs/` | the solver input template |
| `campaign_*.txt`, `rerun_meas.txt`, `twins_meas*.txt`, `full_analysis.txt` | the case files of every run |
| `krylov_test/` | the solver's restarting Krylov driver: run scripts, analysis of its members, and the tested fixes proposed for its eigenvalue routine (`upstream_fixes/`) |
| `results/` | the measured data (records, analysis tables, decomposition outputs, run records) and figures; no plotfiles |
| `archived_records/` | the amplitude records of the September 2026 runs |
| `docs/` | `code_overview.md` (what each tool does, and its verification) and `solver_notes.md` (observations about the solver that bear on the method and on running it) |
