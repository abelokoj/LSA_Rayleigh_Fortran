# krylov_test: the solver's built-in restarting Krylov driver (the restarting variant of Chapter 5)

This folder holds the trial of the restarting procedure of Prof. Sussman's unpublished draft description of the method (steps 1–12), which the solver implements as its built-in LSA driver (`amr.LSA_activate`, `amr.LSA_nsteps_krylov_subspace_method`). Every member restarts from the base state, so the perturbation never grows over a long run. Everything the trial writes stays here, so it can be removed without touching the campaign.

| File | Purpose |
|---|---|
| `krylov_run.sh` | Slurm job: base state, Krylov cycle, then DMD of the members (narrow-band level set and velocity) |
| `lsa_levelset_scale.patch` | solver change making the first member's level-set amplitude an input (`ns.LSA_levelset_scale`, default 2.0 = the former hard-coded value) |
| `runs/<label>/` | solver output of the trial (created on the cluster) |
| `analysis/<label>/` | DMD archives, rank scans and the lossless export of the members (created on the cluster) |
| `slurm-lsa_krylov-<id>.out` | job logs (created on the cluster) |

**Separation from the main work:**

- **Solver:** the patched solver is built in its own git worktree, `~/aii_lsa`. The campaign's executable in `$AII/NS_fluids_lib` is never rebuilt.
- **Output:** nothing is written to `$LSA/runs`, so `analyse_campaign.sh` and `summary.csv` never see the trial.
- **Main code:** the only addition is the reader `run_lsa --krylov`, with `find_lsa_members` and `build_krylov_matrix`. It is used only when that option is given, and is covered by `test_campaign`.

**To run** (from `$LSA`):

```bash
cd $AII && git worktree add ~/aii_lsa $(git rev-parse HEAD)
cd ~/aii_lsa && git apply $LSA/krylov_test/lsa_levelset_scale.patch
cd NS_fluids_lib && make -j16
cd $LSA
SOLVER_EXE=~/aii_lsa/NS_fluids_lib/amr2d.gnu.MPI.ex sbatch --export=ALL krylov_test/krylov_run.sh kry070_16 0.7 16 4.0 60 0.08 0.04
```

**To remove the trial completely:**

```bash
rm -rf $LSA/krylov_test                         # trial output and scripts on the cluster
cd $AII && git worktree remove --force ~/aii_lsa   # the patched solver and its build
```

The `--krylov` reader in `run_lsa` is inert unless that option is given.

## Upstream chain (LSA_activate = 1 then 2)

Since 2026-10-15 the upstream solver completes the chain itself. `LSA_activate=2` reads the members, prints their inner products, calls `fort_process_LSA` (PROB_CPP_PARMS.F90: eq. (67) by SVD with LAPACK, the companion matrix, `ZGGEV`, the most dangerous mode as the eigenvalue of largest real part), and writes the mode as a checkpoint.

- `krylov_upstream.sh` runs that chain on our deck in a separate worktree of the upstream solver.
- `lsa_log_spectrum.py` rebuilds the Gram matrix from the solver log and repeats the decomposition two ways:
  - with `fort_process_LSA`'s rule: no singular-value cutoff, largest real part;
  - with the rule of `svd_and_diagonalize_sussman.f90`: cutoff `rcond`, growing near-stationary mode of largest amplitude.

**Build the upstream worktree** (`USE_LIBLINEAR = TRUE` is the upstream default and needs the bundled LAPACK built first):

```bash
git -C $AII fetch origin && git -C $AII worktree add ~/aii_up origin/master
cd ~/aii_up/NS_fluids_lib/BLASLAPACK/LINEAR && make && cp liblinear.a ../../
cd ~/aii_up/NS_fluids_lib && make -j16 USE_LIBLINEAR=TRUE LIBRARY_LOCATIONS=$PWD
```

**Run** (from `$LSA`):

```bash
SOLVER_EXE=$(ls -t ~/aii_up/NS_fluids_lib/amr2d*.ex | head -1) sbatch --export=ALL krylov_test/krylov_upstream.sh kup070_16 0.7 16 4.0 20 0.25
```

**Points to confirm with Prof. Sussman before reading the numbers:**

1. **Members are rescaled field by field.** In the version of 2026-10-15 each member is rescaled before it is re-added: the level set to a maximum of 2 dx in the band, and the velocity to `LSA_velocity_scale`. The companion relation A V = V S holds only for the unscaled sequence. The printed Gram matrix therefore needs the scale factors (`lsa_log_spectrum.py --scales`), and with different factors for the level set and the velocity no single correction is exact.
2. **The cutoff and the selection rule.** `fort_process_LSA` keeps every positive singular value and selects the largest real part. On a noisy sequence both can pick spurious modes, which `lsa_log_spectrum.py` shows side by side.

**To remove:** `git -C $AII worktree remove --force ~/aii_up`.
