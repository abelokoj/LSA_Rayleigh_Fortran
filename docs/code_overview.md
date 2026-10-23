# The analysis code

This is a reference for what each tool does. `REPRODUCE.md` at the top of the repository gives the commands that reproduce the results of the dissertation; this file explains the tools those commands call.

The method is the Navier–Stokes-based mean-flow-perturbation (NS-MFP) approach of Ranjan, Unnikrishnan and Gaitonde (*J. Comput. Phys.* **403**, 109076, 2020), adapted to the two-phase solver `amrex_implicit_interfaces` (Chapter 4 of the dissertation). The solver is not modified. Each measurement is a pair of runs from the same initial state: one with a small sinusoidal displacement of the interface, and its unperturbed *twin*. Their difference, sampled in time, evolves under the linearized operator of the discrete solver. The tools below generate the decks of such a pair, read the solver's plotfiles, form the difference, and extract the growth rate in three ways: a least-squares fit of the interface mode amplitude, and dynamic mode decomposition (DMD) of the velocity difference and of the level-set difference in a narrow band around the interface.

## Programs

The analysis is written in Fortran 2008 and needs only a Fortran compiler with LAPACK and BLAS: `make` builds the programs into `bin/`, and `make test` builds and runs the tests. Python is used only for figures, tables and independent checks (`plotting/README.md`).

| Program | What it does |
|---|---|
| `bin/make_inputs` | writes the solver decks of a perturbed case and its twin from `inputs/inputs_rayleigh_template.txt` (`single`, `sweep-eps`, `sweep-k`). The two decks differ in exactly one functional line, `ns.radblob` (the amplitude), which a test enforces, because the twin subtraction is valid only under that condition. Options set the wavenumber, amplitude (directly, or as a multiple of the cell size with `--eps-cells`), resolution, wall radius, fluid properties, run length in e-folding times and number of snapshots. |
| `bin/run_lsa` | the analysis driver, with subcommands: `analyse` (twin difference and DMD from the plotfiles or from a saved export), `export` (save the matched twin-difference sequence losslessly), `record` (admissible-window growth rate of an interface-amplitude record, with its start-time sensitivity and the departure from each theoretical relation), `converge` (DMD over grids of start, length, stride and rank), `linearity` (the proportionality test: the twin difference scaled by the amplitude must be independent of the amplitude), `theory` (every dispersion relation used; `--r-wall R` adds the relation with the no-slip wall at r = R), `drift` (drift of the unperturbed twin from its initial state) and `dispersion`. `run_lsa` without arguments lists the options. |
| `bin/dump_amplitude` | the interface-amplitude record of a run: time, Fourier amplitude of the interface radius at the imposed wavenumber, mean radius, the signed cosine and sine coefficients, and the mean of r² |
| `bin/dump_profiles` | the interface profile r(z) at every snapshot, and optionally a field slice, as a plain-text bundle for the figures |
| `bin/svd_and_diagonalize` | the companion-matrix (Krylov) form of the decomposition on an exported sequence, used for the restarting variant (`krylov_test/`) |

## Source modules (`src/`)

| Module | Contents |
|---|---|
| `dispersion_mod` | the inviscid (Rayleigh), Weber, exact viscous (Chandrasekhar), two-fluid inviscid, two-fluid viscous (Tomotika) and wall-bounded two-fluid viscous relations, with exponentially scaled modified Bessel functions |
| `plotfile_mod` | a native reader of AMReX plotfiles (no `yt` needed), for both byte orders and multi-box grids |
| `nsmfp_mod` | pairing of the perturbed and twin snapshots, the twin difference, the narrow band, and the reader of the solver's restarting-driver members |
| `dmd_mod` | SVD-based DMD through LAPACK, with rank truncation, rank scans, mode selection and the reduced-map residual |
| `companion_mod` | the companion-matrix form of the decomposition (weighted least squares, Gram matrix) |
| `interface_mode_mod` | the interface radius from the zero level set and its Fourier coefficients |
| `record_mod` | the admissible window of an amplitude record: after the start-up transient has settled, before the amplitude reaches a/R₀ = 0.1 |
| `proportionality_mod` | the linearity test |
| `archive_mod` | the lossless export of a twin-difference sequence and the DMD archive (singular values, eigenvalues, modes) |
| `bundle_mod` | the plain-text interface-profile format |
| `make_inputs_mod`, `table_mod`, `util_mod` | deck generation, tabular output, utilities |

## Scripts (`scripts/`)

| Script | What it does |
|---|---|
| `run_all.sh` | runs every case of a case file: generates the decks, advances each run to completion (restarting from its newest checkpoint as often as needed), extracts the amplitude record and the profile bundle, and writes a run record (solver commit and checksum, rank count, case line, modules, job IDs) |
| `submit_hpc_array.sh`, `_array_task.sh` | one Slurm array task per case of a case file |
| `submit_campaign.sh` | every case file of the dissertation with its rank count and wall-clock limit, plus the analysis job |
| `_post_task.sh`, `merge_summaries.sh`, `analyse_campaign.sh` | the analysis job: record fits, exports, decompositions, sensitivity grids, the 16-cell linearity test, `runs/analysis/summary.csv` |
| `run_studies.sh` | the band, stride, base-drift and mode studies, from the exports |
| `linearity_08_32.sh` | the linearity test at 8 and 32 cells per radius |
| `refit_records.sh` | recomputes every interface growth rate of a results snapshot and compares it with the stored value |
| `estimate_cost.sh` | the cost and wall-clock time of each case of a case file |
| `install_into_solver_tree.sh` | copies this repository to `NS_fluids_lib/LSA_Rayleigh_Fortran` in the solver tree, where the run scripts expect it |
| `run_nsmfp.sh` | one twin pair on a workstation, from decks to growth rate |

## Case files

Each line of a case file is one run: `label kr0 eps cells_per_r0 rmax n_periods snapshots [make_inputs options]`, where `n_periods` is the run length in e-folding times of the viscous growth rate. A label ending in `.base` is the unperturbed twin. The case files of the dissertation are listed in `REPRODUCE.md`, Section 3.4; `full_analysis.txt` holds the September 2026 runs whose records are in `archived_records/`.

## Verification

`make test` runs eight test programs. They write their scratch files to a private directory under `$TMPDIR` (or `/tmp`) whose name contains a space, so that two runs on a shared login node cannot collide and every test that reads plotfiles also proves that paths with spaces work.

| Test program | What it checks |
|---|---|
| `test_dispersion` | the modified Bessel functions against quadrature of their integral representation and against scipy values; every dispersion relation against its limits (inviscid, long-wave, wall at infinity, vanishing viscosity, the stable band) and against independent solutions; the wall-bounded relation against the values of `plotting/confined_dispersion.py` |
| `test_dmd` | recovery of planted growth rates and frequencies, rank truncation, degenerate subspaces, input validation |
| `test_proportionality` | the linearity test, including nonlinear contamination against round-off, mismatched snapshot times and the twin invariant of generated decks |
| `test_end_to_end` | synthetic plotfiles in the genuine on-disk format read back through the production parser, differenced and decomposed, recovering the exact growth rate; a planted base-state drift two orders of magnitude above the perturbation makes the test fail if the twin subtraction is omitted or the snapshots are mispaired |
| `test_interface` | the interface radius and amplitude against planted interfaces, multi-box and big-endian plotfiles, decoy names, holes in the box layout, single precision |
| `test_bundle` | the profile-bundle format, written and read back |
| `test_campaign` | the admissible-window analysis against a known start-up, the DMD residual, the lossless export and archive, the volume-weighted inner product, the signed coefficients, the deck options, the restarting-driver reader |
| `test_companion` | the companion-matrix solver on Krylov sequences built from known eigenvalues |

The plotfile reader carries the most risk, since there is no Fortran equivalent of `yt`. During development a plotfile written by the Fortran writer was read correctly by `yt`, and the writer and reader round-trip in both byte orders.

Some code comments refer to the "September 2026 evidence report". It is the report of the first (September 2026) study, `NS_fluids_lib/analysis_report.typ` in the solver fork `abelokoj/amrex_implicit_interfaces` at commit `dddb827c`; the dissertation supersedes it.
