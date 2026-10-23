# Reproducing the results of the dissertation

This file gives the commands that reproduce every table and data figure of Chapters 5 and 6 of the dissertation, at three levels of effort. Each level is self-contained, and each lists the result to expect. All commands are run from the top of this repository unless a step says otherwise.

| Level | What it reproduces | Needs | Time |
|---|---|---|---|
| 1. Check the code | the analysis code builds and passes its tests; the theoretical growth rates | gfortran, LAPACK | 5 minutes |
| 2. Recompute from the stored data | every growth rate, table and data figure, from the amplitude records and decomposition outputs in `results/` | Level 1, Python 3 with numpy, scipy, matplotlib | 15 minutes |
| 3. Rerun the simulations | the stored data itself, from the solver | a Slurm cluster, the solver `amrex_implicit_interfaces` | several days of queue and run time |

Level 2 is the one most readers will want: it starts from the raw interface-amplitude records and the decomposition outputs of the runs and recomputes every number the dissertation quotes, without the solver.

## Level 1. Check the code

On any Linux machine, or WSL on Windows (on macOS, add `LAPACK_LIBS="-framework Accelerate"` to `make`):

```bash
sudo apt-get install -y gfortran make liblapack-dev libblas-dev
git clone https://github.com/abelokoj/LSA_Rayleigh_Fortran.git
cd LSA_Rayleigh_Fortran
make            # bin/run_lsa, bin/make_inputs, bin/dump_amplitude, bin/dump_profiles, bin/svd_and_diagonalize
make test       # must end with: ALL TEST PROGRAMS PASSED
```

`make test` checks the dispersion relations against limits, published values and independent scipy solutions (including the wall-bounded relation of Appendix D), the decomposition on sequences with known eigenvalues, the plotfile reader and writer, the interface-amplitude fit, the linearity test, the campaign tools, and the companion-matrix solver.

The theoretical growth rates used throughout, for the campaign fluids (Oh = 0.02, ρ₂/ρ₁ = 1.225 × 10⁻³, μ₂/μ₁ = 0.013, the defaults of `run_lsa`):

```bash
./bin/run_lsa theory --kr0 0.7 --r-wall 4
```

Expected at k R₀ = 0.7: inviscid 0.343327, Chandrasekhar (Eq. 3.10) 0.329367, two-fluid viscous (Tomotika) 0.328689, two-fluid viscous with the no-slip wall at 4 R₀ 0.328617. With `--r-wall 1.5` the last value is 0.318305 (a 3.16 % reduction).

## Level 2. Recompute every table and figure from the stored data

```bash
python3 -m pip install -r plotting/requirements.txt
mkdir -p figures          # every output of this level goes here (ignored by git)
```

### 2.1 Every interface growth rate, from the amplitude records

```bash
scripts/refit_records.sh
```

This fits every amplitude record in `results/followup_2026-10-22/records/` with `bin/run_lsa record`, exactly as the analysis job did, and prints the stored and recomputed rate of each of the 55 cases side by side. Expected: a largest difference of 0 %. The full report of each fit (window, e-folding times, start-time spread, departure from each relation) is written to `runs/refit/followup_2026-10-22/<case>.txt`.

### 2.2 The tables of Chapters 5 and 6

```bash
python3 plotting/thesis/print_thesis_tables.py
```

This prints Tables 5.3 (the 16-cell sweep), 5.4 (mesh study), 5.5 (stable mode), 5.6 (sampling stride), 6.1 (the three reference relations), 6.2 (linearity on three grids) and 6.3 (weighting of velocity and level set), computed from `results/` and from `plotting/confined_dispersion.py` for the wall-bounded relation. Every entry should match the dissertation. The start-time spreads are printed with three decimals; the dissertation rounds them to two.

The Chapter 6 comparison of every confined and viscous-ambient run with both relations:

```bash
python3 plotting/compare_confined.py results/followup_2026-10-22/analysis/summary.csv figures/confined_comparison.csv
```

Expected: the table of `results/ch6_figures_2026-10-22/confined_comparison.csv`, with every rate whose window exceeds one e-folding time within ±0.2 % of the wall-bounded relation (`conf15k080_16`, a 0.7 e-fold window, is +0.42 %).

### 2.3 The independent checks of the theory

```bash
cp -r results/ch6_figures_2026-10-22/theory figures/          # a working copy: the next command rewrites the confined tables in the folder it is given
cd plotting
python3 make_confined_tables.py ../figures/theory   # Appendix D: agreement with the Fortran tables, and the wall at 20 R0
python3 check_inviscid_limit.py                     # Appendix D: vanishing-viscosity limit
cd ..
python3 plotting/check_mode_shapes_confined.py results/followup_2026-10-22/analysis   # Chapter 5: mode shapes
```

Expected: a largest relative difference of 4.8 × 10⁻¹³ from the Fortran tables over 234 values, and 3.7 × 10⁻¹⁴ between the relation with the wall at 20 R₀ and the unbounded one, with confined tables in `figures/theory/` whose values equal the stored ones; the gap to the inviscid confined relation shrinking in proportion to √Oh down to about 10⁻⁵ at Oh = 10⁻¹⁰; and rms departures of the leading velocity mode from the wall-bounded eigenfunction of 0.81, 0.86, 0.89 and 0.92 % at k R₀ = 0.4, 0.5, 0.6 and 0.7 (5.5–9.8 % against the unbounded one).

The unbounded theory tables (`theory_vs_k.csv`, `theory_vs_mu_ratio.csv`) were written by the Fortran library; to regenerate the first with the Level 1 build (one wavenumber per call, since `run_lsa` accepts at most 64 values per call):

```bash
mkdir -p figures/theory
echo "kr0,inviscid,weber,chandrasekhar,two_fluid_inviscid,two_fluid_viscous" > figures/theory/theory_vs_k.csv
for k in $(python3 -c "print(' '.join(f'{0.05+0.005*i:.3f}' for i in range(190)))"); do
  ./bin/run_lsa theory --kr0 $k --csv figures/theory/tmp.csv > /dev/null && tail -1 figures/theory/tmp.csv | tr -d ' ' >> figures/theory/theory_vs_k.csv
done
rm figures/theory/tmp.csv
```

Its values equal those of `results/ch6_figures_2026-10-22/theory/theory_vs_k.csv`.

### 2.4 The figures

Each command writes PNG, PDF (and for most, SVG) files into `figures/`, under the name given in the first column.

| Figure in the dissertation | File written | Command |
|---|---|---|
| dispersion relation of the 16-cell sweep (Ch. 5) | `campaign_dispersion` | `python3 plotting/plot_campaign_results.py results/campaign_2026-10-09 --outdir figures` |
| stable mode at k R₀ = 1.1 (Ch. 5) | `campaign_stable_band` | same command |
| interface amplitude, local rate and fitting windows (Ch. 5) | `ch5_window` | `python3 plotting/thesis/make_extra_figures.py results/campaign_2026-10-09 results/estimators_2026-10-22/weighting_test_output_2026-10-22.txt plotting figures` |
| proportionality to the amplitude (Ch. 5) | `ch5_proportionality` | same command |
| interface profiles (Ch. 5) | `ch5_profiles` | same command |
| mesh study (Ch. 5) | `ch5_mesh` | same command |
| sampling stride (Ch. 5) | `ch5_stride` | same command |
| weighting of velocity and level set (Ch. 6) | `ch6_weighting` | same command |
| singular values, eigenvalues, rank and residual (Ch. 5) | `ch5_dmd_evidence` | `python3 plotting/thesis/make_dmd_evidence_figures.py results/followup_2026-10-22/analysis plotting figures` |
| leading mode against the eigenfunction (Ch. 5) | `ch5_mode_shapes` | same command |
| perturbation field with the interface (Ch. 5) | `gallery/gallery_field_meas16` | `python3 plotting/plot_campaign_gallery.py results/campaign_2026-10-09 --two-fluid results/two_fluid_2026-10-16/records --export results/exports/meas16 --export-case meas16 --no-gif --outdir figures/gallery` |
| confinement (Ch. 6) | `ch6_confinement` | `python3 plotting/plot_followup_ch6.py results/followup_2026-10-22/analysis/summary.csv --theory results/ch6_figures_2026-10-22/theory --outdir figures` |
| viscous ambient (Ch. 6) | `ch6_viscous_ambient` | same command |

The gallery command also draws the rest of the campaign's figure set (growth curves, local rates, the sweep, break-up, the Ohnesorge plane) and combines them into `figures/gallery/all_figures.pdf`; without `--no-gif` it also writes the animations, which takes about six minutes. The panels of `ch5_mesh` and `ch5_stride` plot the values of Tables 5.4 and 5.6 as typed into the script; those values are recomputed from the data by `print_thesis_tables.py` (Section 2.2). The schematic figures of Chapters 3 and 4 (the column, the MAC grid, MOF and CMOF, the NS-MFP variants) are drawings, not data, and are not part of this repository.

### 2.5 What Level 2 does not recompute

The decompositions themselves (singular values, eigenvalues, modes) are read from their stored outputs in `results/*/analysis/<case>/dmd`, `dmd_ls` and `studies/`. Repeating them needs the lossless twin-difference exports (`perturbation.bin`, up to 70 MB per case, about 1 GB in all), which stayed on the cluster and are listed in each snapshot's `SKIPPED.txt`; Level 3 regenerates them. The one export included is `results/exports/meas16/`, for the perturbation-field figure. Likewise the weighting test (Table 6.3) and the field-direction linearity test (Section 5.5, `results/linearity_fields_2026-10-23/`) are read from their stored outputs and repeated in Level 3.

## Level 3. Rerun the simulations

The runs used the solver `amrex_implicit_interfaces` at commit `dddb827c` (recorded, with the MD5 checksum `eeeee45161e1267ea7b36012ec1531d8` of the executable, in every `results/*/run_records/<case>/manifest.txt`), built with `gnu/13` and `openmpi/4.1.6` on the FSU Research Computing Center cluster and run under Slurm. The scripts assume Slurm; on another scheduler, `scripts/run_all.sh` runs a case file sequentially with `mpirun`.

### 3.1 Obtain and build the solver

Commit `dddb827c` is commit `0b8df42` of Prof. Sussman's repository with the September analysis package added; the solver source (every `.F90`, `.cpp`, `.H` file and the makefiles) is identical in the two. Either gives the same executable:

```bash
export AII=$HOME/amrex_implicit_interfaces
git clone https://github.com/msussman42/amrex_implicit_interfaces.git $AII
cd $AII && git checkout 0b8df42
# or, the exact commit of the runs:
#   git clone https://github.com/abelokoj/amrex_implicit_interfaces.git $AII && cd $AII && git checkout dddb827c
#   rm -rf NS_fluids_lib/LSA_Rayleigh_Fortran NS_fluids_lib/LSA_Rayleigh_Python   # the September package it carries; Section 3.2 installs the current one

module purge; module load gnu/13; module load openmpi      # the toolchain of the runs; adapt to your cluster
cd $AII/NS_fluids_lib && make -j16                          # amr2d.gnu.MPI.ex (not a DEBUG build)
```

Use `0b8df42` and not the current head: at some later commits (for example `b014c1f`) every fluid-only calculation stops at start-up with `expecting num_FSI_outer_sweeps>=2` (`docs/solver_notes.md`, Section 2). AMReX is vendored in the repository; no separate clone is needed.

### 3.2 Install this repository into the solver tree, build and test

The run scripts expect to sit in `NS_fluids_lib/LSA_Rayleigh_Fortran`, with the solver one level up:

```bash
cd $HOME && git clone https://github.com/abelokoj/LSA_Rayleigh_Fortran.git
cd LSA_Rayleigh_Fortran && scripts/install_into_solver_tree.sh
export LSA=$AII/NS_fluids_lib/LSA_Rayleigh_Fortran
cd $LSA && make -j8 && make test
```

### 3.3 Smoke test (ten minutes)

```bash
cd $LSA
printf 'smoke_16 0.7 2.5e-3 16 4.0 4.3 160 --eps-cells 0.04\n' > smoke_case.txt
salloc -t 00:10:00 -n 8 -A <account> bash -c 'module purge; module load gnu/13; module load openmpi; LAUNCHER=srun timeout 300 ./scripts/run_all.sh smoke_case.txt 8'
grep -m3 '^STEP' runs/smoke_16/solver.out && rm -rf runs/smoke_16 runs/records/smoke_16.dat* smoke_case.txt
```

Lines beginning `STEP` mean the decks, the solver and the launcher work together. Anything else is explained in `runs/smoke_16/solver.err`; `docs/solver_notes.md` lists the failures seen in practice.

### 3.4 Submit every run

```bash
cd $LSA
ACCOUNT=<account> scripts/submit_campaign.sh
```

This submits each case file as a Slurm job array, with the rank count and wall-clock limit it was run with, and an analysis job that starts when all arrays have ended. The case files and what they produce:

| Case file | Runs | Ranks | Results in the dissertation |
|---|---|---|---|
| `campaign_16_wide.txt` | sweep k R₀ = 0.2, 0.3, 0.4 with twins | 32 | Tables 5.3, 5.6, 6.1 |
| `campaign_16.txt` | sweep k R₀ = 0.5–0.9 with twins; ε = 5 × 10⁻³ and 10⁻² at k R₀ = 0.7; time step halved; wall at 6 R₀; gas density and gas viscosity ÷ 10; stable mode k R₀ = 1.1 | 16 | Tables 5.3–5.6, 6.1, 6.2 |
| `rerun_meas.txt` | mesh study at k R₀ = 0.7, ε = 5 × 10⁻³, on 8, 16 and 32 cells, with twins | 16 | Table 5.4, mode shapes |
| `campaign_next_16.txt` | ε = 10⁻³ near the cut-off; wall at 3, 2 and 1.5 R₀; μ₂/μ₁ = 0.1, 0.3, 1 | 16 | Chapter 6 |
| `campaign_amp08.txt`, `campaign_amp32.txt` | three amplitudes at 8 and at 32 cells, with one twin each | 8, 16 | Table 6.2 |
| `campaign_check.txt` | wall at 1.5 R₀ on 32 cells; μ₂/μ₁ = 1 with the wall at 6 and 8 R₀ | 16 | Chapter 6, Table 6.3 |

`scripts/estimate_cost.sh <case file> <ranks>` prints the expected cost of each file. The 16-cell cases take from a few hours to about half a day each; the 32-cell cases about a day each on 16 ranks. A task that reaches its wall-clock limit is not lost: submitting the same command again skips finished cases and resumes the others from their newest checkpoint. Progress: `squeue --me`, and `tail runs/array/run_<job>_<task>.log`.

The stored `meas16.base` and `meas32.base` were made in October as twins of September runs, with `twins_meas16.txt` and `twins_meas32.txt` (`--twin-of`, which copies the partner's own deck). `rerun_meas.txt` states the same configuration as ordinary twin pairs, so that it can be run on an empty account; see its header.

### 3.5 The analysis

The analysis job (`scripts/analyse_campaign.sh`, through `scripts/_post_task.sh`) fits every record, exports the matched twin differences, decomposes them, runs the sensitivity grids and the 16-cell linearity test, and writes `runs/analysis/summary.csv`. Its log is `runs/array/post-<job>.out` and must end with `status 0`. Then, on a compute node or as batch jobs:

```bash
cd $LSA
sbatch -A <account> -n 1 -c 4 --mem=16G -t 04:00:00 -o runs/array/studies-%j.out scripts/run_studies.sh all       # band, stride, drift and mode studies
sbatch -A <account> -n 1 -c 2 --mem=8G  -t 02:00:00 -o runs/array/linearity-%j.out scripts/linearity_08_32.sh 08 32   # Table 6.2, 8 and 32 cells
sbatch -A <account> -n 1 -c 4 --mem=32G -t 04:00:00 -o runs/array/linfields-%j.out scripts/linearity_fields.sh        # Section 5.5: scaled fields compared in direction (velocity, level set whole and band 2), 8, 16, 32 cells
python3 plotting/weighting_test.py runs/analysis sw040_16 sw050_16 sw060_16 sw070_16 sw080_16 meas16 meas32 conf15_32 visc100r8_16 > runs/analysis/weighting_test_output.txt   # Table 6.3
```

### 3.6 Compare with the stored results

```bash
python3 plotting/thesis/compare_summaries.py runs/analysis/summary.csv
python3 plotting/thesis/print_thesis_tables.py --snapshot runs --weighting runs/analysis/weighting_test_output.txt
scripts/refit_records.sh runs
```

The first prints the stored and the new interface, level-set and velocity rates of every case side by side; the second prints the tables of the dissertation from the new runs. The figure commands of Section 2.4 accept the new data in place of the stored snapshots: `runs` in place of `results/campaign_2026-10-09` and `runs/analysis` in place of `results/followup_2026-10-22/analysis`.

Bitwise reproducibility of the runs has not been tested. A different compiler, MPI library or rank count changes the round-off; expect the interface rates to agree closely with the stored ones, and the decomposition rates, which are computed from the difference of two runs, to differ more, most in the cases whose velocity decomposition is already flagged (`sw080_16`, `sw090_16`, `meas16`). Within a twin pair the rank count must be the same, which the case files and `scripts/submit_campaign.sh` ensure. The September runs (`full_analysis.txt`: 8 cells, ε = 2 × 10⁻²) are not needed for any table; they are kept in `archived_records/` and discussed in the dissertation as history.

### 3.7 The restarting variant

The comparison with the solver's built-in restarting Krylov driver (Chapter 5) uses a patched build in its own worktree; `krylov_test/README.md` gives the commands. Its per-member results are in `results/restarting_members_2026-10-23/`.

## Where things are

| Path | Contents |
|---|---|
| `results/followup_2026-10-22/` | the complete analysis of every run: `analysis/summary.csv` (one line per case), `analysis/records/` (the fit of every record), `analysis/<case>/` (decomposition outputs, sensitivity grids, studies), `records/` (the amplitude records), `run_records/` (solver commit, checksum, ranks and case line of every run) |
| `results/campaign_2026-10-09/` | the same for the first campaign alone, with the interface-profile bundles used by the figures |
| `results/two_fluid_2026-10-16/`, `ch6_figures_2026-10-22/`, `estimators_2026-10-22/`, `gallery_2026-10-16/`, `restarting_members_2026-10-23/` | the two-fluid comparison, the Chapter 6 figures and theory tables, every estimator for every case with the quoting rule and the weighting test, the figure gallery, and the restarting-variant members |
| `results/exports/meas16/` | the one lossless twin-difference export included (for the perturbation-field figure) |
| `archived_records/` | the September 2026 amplitude records |
| `docs/code_overview.md` | what each tool does, and how the code is verified |
| `docs/solver_notes.md` | observations about the solver that bear on the method and on running it |
