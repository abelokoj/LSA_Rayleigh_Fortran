# Notes on the solver and on running it

Each observation here comes from building and running `amrex_implicit_interfaces`, not from reading its source or from synthetic tests. Each is stated with the symptom by which it is recognised, so that it can be identified again rather than rediscovered. `REPRODUCE.md` gives the commands; this file explains the choices behind some of them.

## 1. Building the solver

AMReX is vendored in the repository at `amrex-master`, at the same commit as the solver, so no separate clone is needed despite the suggestion in the `GNUmakefile` header. The build needs MPI (`libopenmpi-dev openmpi-bin` on Debian and Ubuntu; on the FSU RCC cluster the runs used the modules `gnu/13` and `openmpi/4.1.6`). From `NS_fluids_lib`, `make -j16` produces `amr2d.gnu.MPI.ex`, of about 109 MB. A report of `amr2d.gnu.DEBUG.ex is up to date` means `DEBUG = TRUE` is set in `GNUmakefile`: the MPI executable is then neither built nor refreshed, and a debug build is several times slower.

If a build is interrupted, `make` may leave an object file of zero length, which it then treats as current. The symptom appears much later, at the link stage, as an undefined reference to a symbol whose source apparently compiled, for example `amrex::FillPatchTower`. The remedy is

```bash
find tmp_build_dir -name '*.o' -size 0 -delete && make -j16
```

## 2. A defect that stops every fluid-only calculation, and the solver commit

At commit `cbcee7c`, every calculation without an ice, rigid or elastic material stopped during initialisation with `expecting num_FSI_outer_sweeps>=2`, the decks distributed with the solver (`run2d/inputs.growthrate.LSA`, `run2d/inputs.capillary`) included. The quantity `num_FSI_outer_sweeps` is derived, not an input: it starts at one and is incremented for each ice or FSI material, so for fluid-only problems it is always one, and an assertion that it be at least two rejects the whole class. Of the assertions introduced with the elastic-material work, the one reached by fluid-only problems lies in the fall-through branch.

It was corrected upstream at commit `0b8df42`, where the assertion is guarded by a test on the primary material; both a distributed deck and a generated deck then advance normally. At commit `b014c1f` it was present again, for distributed and generated decks alike. The runs of the dissertation therefore used `0b8df42`, as commit `dddb827c` of the fork `abelokoj/amrex_implicit_interfaces`, which adds the September analysis package to `0b8df42` and leaves the solver source unchanged; every run record in `results/` gives `dddb827c` and the same executable checksum. Whether later upstream commits carry the defect has not been checked, so a later commit should be tested with a distributed deck before use:

```bash
cd $AII/NS_fluids_lib/run2d && mpirun -n 1 ../amr2d.gnu.MPI.ex inputs.growthrate.LSA     # lines beginning STEP within seconds
```

Running a distributed deck is the most direct way to tell a solver failure from a fault in a generated deck.

## 3. Output naming

The solver does not honour `amr.plot_file`. With `ns.visual_nddata_format = 1` it writes directories named `nddataPLT00000023`, together with `MOF_PLT*` directories and Tecplot files, rather than `plt00023`. A tool that expects the documented name finds no output and reports an empty run rather than a naming mismatch. `nddataPLT` is the default prefix of the tools here (`--prefix` overrides it). The plotfiles are in the standard AMReX format with 46 components, among them `x_velocity`, `y_velocity`, `PRES_MG`, the volume fractions `F01` and `F02`, and the level sets `L0101` and `L0202`.

## 4. The quiescent column is an exact discrete equilibrium

The Fortran plotfile reader and `yt` agree to the last digit on genuine solver output. On the unperturbed column the largest velocity after 23 steps is 6.0 × 10⁻¹⁵: the base state stays at rest to machine precision. This is the premise of replacing the constraining body force of NS-MFP by a twin run (Chapter 4 of the dissertation); the twin then removes the remaining drift (`scripts/run_studies.sh drift`, Chapter 5).

## 5. Choice of observable

**The velocity difference is zero at the first snapshot.** The perturbation is imposed on the interface, so both runs of a pair start from rest. A decomposition that ranks modes by their projection on the first snapshot then returns every amplitude as zero; the first snapshots must be discarded (`--skip`), and in practice the decomposition window starts after the start-up transient, at the same time as the interface fit.

**The level set over the whole domain dilutes the signal.** The level set is a signed distance, so the twin difference is about ε sin(kz) throughout the domain, while only the neighbourhood of the interface takes part in the instability; in one run its norm stayed constant to 0.6 % while the interface amplitude grew by 37 %. Away from the interface the re-distanced level set also carries non-growing components. The level set is therefore decomposed in a narrow band, |r − R₀| < 2Δr.

**The interface amplitude is the most direct estimate.** The Fourier amplitude of the interface radius at the imposed wavenumber is the quantity the dispersion relation predicts; for r = R₀ + a(t) cos kz, linear theory gives a(t) = a(0) e^{st}, so a log-linear fit returns the growth rate. It is used alongside the decompositions, which supply the spectrum and the mode shapes; Chapter 4 states which estimate is quoted.

## 6. Amplitude against the mesh

An amplitude of 10⁻³ on a mesh of spacing 0.125 is under one per cent of a cell, below the scale of the interface reconstruction, and the run returned noise. The amplitude must be scaled against the mesh rather than against the state (the ε₀ = 10⁻⁵ of Ranjan et al. is a fraction of a free-stream state variable, not an interface displacement). The October 2026 campaign uses ε = 0.04 Δr (`--eps-cells 0.04`, 2.5 × 10⁻³ at 16 cells per radius), and the linearity test of Chapters 5 and 6 confirms the response is proportional over 2.5 × 10⁻³ to 10⁻².

## 7. Restarts and the order of snapshots

Checkpoint restarts work, and AMReX accepts the settings as command-line overrides, so the deck need not be edited (`scripts/run_all.sh` does this). The checkpoint interval must be smaller than the number of steps one allocation completes, or no checkpoint is written and each restart repeats from the beginning while the growing number of plotfiles gives every appearance of progress. Restarts write an extra plotfile, so the number of output directories exceeds `max_step / plot_int`; pairing the twins is by step number, which is unaffected. A time series must be ordered by time, not by step number or by the order a directory listing returns.

## 8. What an amplitude record does not carry

The amplitude record of a case is enough to fit a growth rate and to draw the growth curve and the local rate, but not to draw the interface shape or a field slice. The sweep scripts therefore write a profile bundle (`bin/dump_profiles`) for every case while the plotfiles are present: the interface profile at every snapshot in plain text, a few hundred kilobytes against tens of gigabytes of plotfiles. The interface figures of the dissertation are drawn from these bundles.

A Fortran detail found in making the format portable: written with the default descriptor `es24.16`, a real whose exponent needs three digits loses its exponent letter (4.66e-310 appears as `4.6598606362008991-310`). Fortran reads that back correctly, but no other language does. The writers use `es25.16e3`, and the reader also accepts the two-digit form.

## 9. The start-up transient and the fitting window

The September 2026 runs (ε = 2 × 10⁻², 8 cells per radius, `archived_records/`) gave growth rates about four per cent below the viscous relation. That departure was first read as discretisation error, but it came from the fitting window: the initial condition perturbs the interface and not the velocity, so the state starts as a mixture of the growing mode and decaying ones, and a fit that includes the start-up measures the mixture and comes out low. At ε = 2 × 10⁻² the amplitude reaches a/R₀ = 0.1, beyond which the response is no longer linear, after at most 0.24 e-folding times of settled growth, so those runs have no clean window at all (`archived_records/README.md`).

The October 2026 campaign was designed around this. The amplitude is eight times smaller (ε = 2.5 × 10⁻³), so that for k R₀ ≤ 0.7 more than one e-folding time of settled growth fits below a/R₀ = 0.1, and every rate is fitted over an admissible window that opens when the local growth rate has settled and closes at a/R₀ = 0.1 (`run_lsa record`), with the spread of refits from later start times as its uncertainty. For k R₀ = 0.2 to 0.7 the remaining departure of the interface rates from the single-fluid viscous relation, −0.29 to −0.45 %, is mostly the surrounding gas and the outer wall, which that relation omits: against the two-fluid relation with the wall it falls to −0.01 to −0.19 % (Table 6.1, Chapter 7). The two wavenumbers nearest the cut-off, whose windows are shorter than one e-folding time, are discussed separately there.
