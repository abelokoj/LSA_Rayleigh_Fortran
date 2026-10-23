#!/bin/bash
# =====================================================================
# krylov_run.sh -- run the solver's built-in restarting Krylov LSA driver for one case, then decompose the members.
# =====================================================================
# This is the restarting procedure of Prof. Sussman's unpublished draft description of the method (steps 1-12), which the solver implements as its built-in LSA driver: from a base state S0 the solver advances M steps to S_alt = N^M(S0); every Krylov member k then restarts from S0 with the previous difference added, advances the same M steps, and the new difference dS_k = N^M(S0 + dS_{k-1}) - S_alt is the next member. The perturbation is never allowed to grow over a long run, so the sequence can be as long as wanted without leaving the linear range, unlike the continuous twin (whose clean window is 0.6-0.8 e-folds). The solver implements it already (amr.LSA_activate, amr.LSA_nsteps_krylov_subspace_method; NavierStokes3.cpp and main.cpp); this script drives it and hands the members to run_lsa.
#
# Everything this test writes stays inside $LSA/krylov_test/ (runs/, analysis/, Slurm logs), so it never mixes with the campaign in $LSA/runs and can be removed with  rm -rf $LSA/krylov_test ~/aii_lsa  (see README.md here).
#
# Submit from $LSA (it is a Slurm batch script):
#
#     sbatch krylov_test/krylov_run.sh                               # k r0 = 0.7, 16 cells/r0, 60 members
#     sbatch krylov_test/krylov_run.sh kry070_16 0.7 16 4.0 60 0.08 0.04
#              arguments: LABEL KR0 CELLS RMAX NMEMBERS TSAMPLE LSSCALE
#
#   NMEMBERS  number of Krylov members (amr.LSA_nsteps_krylov_subspace_method)
#   TSAMPLE   sampling interval in capillary times; M = round(TSAMPLE / fixed_dt) steps per member
#   SKIP      (environment) members dropped before the DMD, default 0; the first trial showed the sequence needs about 10 t_c to settle.
#   LSSCALE   level-set amplitude of the first member, in cells (ns.LSA_levelset_scale). The solver hard-codes 2 cells,
#             i.e. 0.125 r0 at 16 cells/r0, far outside the linear range (a/r0 <= 0.1); 0.04 gives 0.0025 r0, the amplitude
#             of the October campaign. Needs krylov_test/lsa_levelset_scale.patch; without it the script stops unless
#             ALLOW_UNPATCHED=1 is set.
#
# Notes from reading the solver at commit dddb827 (October 2026):
#   - only the first member is normalised; later members are raw differences, cut to a narrow band |phi| < 3 dx in the
#     level set. At sigma ~ 0.33 and TSAMPLE = 0.08 the members grow by about 2.7 % each, so 60 members grow ~5x, which
#     stays linear from 0.0025 r0.
#   - the first member has zero velocity (the default eigenvector perturbs only the level set, LS = dx sin(2 pi z/yblob)),
#     so ns.yblob must be the wavelength, which make_inputs ensures.
#   - members are written as AMReX plotfiles nddataPLT<step>LSA<k> (amr.LSA_plot_index = 0; other indices write Tecplot
#     files and appear to fail in member 0).
# ---------------------------------------------------------------------
#SBATCH --job-name=lsa_krylov
#SBATCH --ntasks=16
#SBATCH --time=06:00:00
#SBATCH --account=genacc_q
#SBATCH --output=krylov_test/slurm-%x-%j.out
#SBATCH --error=krylov_test/slurm-%x-%j.err
#SBATCH --mail-type=END,FAIL
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
TEST="$PWD/krylov_test"          # all output of this test lives here
mkdir -p "$TEST/runs" "$TEST/analysis"
LABEL="${1:-kry070_16}"; KR0="${2:-0.7}"; CELLS="${3:-16}"; RMAX="${4:-4.0}"
NMEM="${5:-60}"; TSAMPLE="${6:-0.08}"; LSSCALE="${7:-0.04}"
EXE="${SOLVER_EXE:-../amr2d.gnu.MPI.ex}"
EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"
NPROCS="${SLURM_NTASKS:-16}"
SKIP="${SKIP:-0}"                # members dropped from the start of the decomposition (startup transient of the power iteration)
SOLVER_OPTS="${SOLVER_OPTS:-mac.mac_abs_tol=1.0e-10 mac.visc_abs_tol=1.0e-10}"

unset SLURM_STEP_ID SLURM_STEPID SLURM_STEP_NUM_TASKS SLURM_STEP_NODELIST || true
set +u
if type module >/dev/null 2>&1; then module purge; module load gnu/13; module load openmpi; fi
set -u

[[ -x "$EXE" ]] || { echo "error: solver $EXE not found" >&2; exit 1; }
[[ -x bin/make_inputs && -x bin/run_lsa ]] || { echo "error: run make first" >&2; exit 1; }
if ! grep -q "LSA_levelset_scale" "$EXE" 2>/dev/null; then
  if [[ -z "${ALLOW_UNPATCHED:-}" ]]; then
    echo "error: the solver does not know ns.LSA_levelset_scale; apply krylov_test/lsa_levelset_scale.patch and rebuild," >&2
    echo "       or set ALLOW_UNPATCHED=1 to run with the hard-coded first-member amplitude of 2 cells (nonlinear)." >&2
    exit 1
  fi
  echo "WARNING: unpatched solver: first member amplitude is 2 cells"
fi

DIR="$TEST/runs/$LABEL"
mkdir -p "$DIR/decks"
# The base deck is the unperturbed deck of make_inputs (radblob = 0, yblob = lambda, fixed dt); n-periods and snapshots do not matter here.
./bin/make_inputs single --kr0 "$KR0" --eps 1e-3 --cells-per-r0 "$CELLS" --rmax "$RMAX" \
    --n-periods 4 --snapshots 100 --outdir "$DIR/decks" > "$DIR/decks/plan.txt"
deck=$(grep -E '/inputs\.k[^ ]*\.base$' "$DIR/decks/plan.txt" | head -1)
[[ -f "$deck" ]] || { echo "error: base deck not generated" >&2; exit 1; }
cp "$deck" "$DIR/inputs"
dt=$(awk '$1=="ns.fixed_dt"{print $3; exit}' "$DIR/inputs")
M=$(awk -v t="$TSAMPLE" -v d="$dt" 'BEGIN{m=int(t/d+0.5); if (m<1) m=1; print m}')
echo "case $LABEL: k r0 = $KR0, $CELLS cells/r0, rmax $RMAX; fixed_dt $dt, M = $M steps per member ($TSAMPLE t_c), $NMEM members, first-member level set ${LSSCALE} cells"

LAUNCH=(srun --cpu-bind=none -n "$NPROCS")
COMMON=($SOLVER_OPTS amr.max_grid_size=16 amr.LSA_nsteps_krylov_subspace_method=$NMEM amr.LSA_plot_index=0 \
        ns.LSA_levelset_scale=$LSSCALE)

{
  echo "--- $(date '+%Y-%m-%d %H:%M:%S') krylov_run.sh"
  echo "label=$LABEL kr0=$KR0 cells=$CELLS rmax=$RMAX members=$NMEM M=$M dt=$dt tsample=$TSAMPLE lsscale=$LSSCALE"
  echo "solver=$EXE md5=$(md5sum "$EXE" | cut -d' ' -f1) commit=$(git -C "$(dirname "$EXE")" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "ranks=$NPROCS slurm_job=${SLURM_JOB_ID:-none}"
} >> "$DIR/manifest.txt"

cd "$DIR"
if [[ ! -f chk00000LSA00000/Header ]]; then
  echo "[$(date +%T)] base state: writes chk00000LSA00000"
  "${LAUNCH[@]}" "$EXE" inputs "${COMMON[@]}" amr.LSA_activate=0 max_step=1 amr.check_int=1 amr.plot_int=1 \
      < /dev/null > base.out 2> base.err
fi
[[ -f chk00000LSA00000/Header ]] || { echo "error: chk00000LSA00000 was not written; see $DIR/base.err" >&2; exit 1; }

echo "[$(date +%T)] Krylov cycle: $NMEM members of $M steps"
"${LAUNCH[@]}" "$EXE" inputs "${COMMON[@]}" amr.restart=chk00000LSA00000 amr.LSA_activate=1 \
    amr.LSA_initial_levelSteps=0 max_step=$M stop_time=1.0e6 amr.plot_int=$M amr.check_int=$M \
    < /dev/null > lsa.out 2> lsa.err || echo "WARNING: the LSA run ended with an error; see $DIR/lsa.err"
cd - > /dev/null

n=$(ls -d "$DIR"/nddataPLT*LSA* 2>/dev/null | wc -l)
echo "[$(date +%T)] $n member plotfiles written"

# Decomposition of the members: narrow-band level set (the observable that converged on the continuous-twin data) and velocity for comparison.
out="$TEST/analysis/$LABEL"; mkdir -p "$out"
./bin/run_lsa analyse --krylov "$DIR" --fields L0101 --band 2 --weight volume --kr0 "$KR0" --skip "$SKIP" \
    --rank-scan 2:16 --archive "$out/dmd_krylov_ls" > "$out/dmd_krylov_ls.log" 2>&1 || true
./bin/run_lsa analyse --krylov "$DIR" --fields x_velocity,y_velocity --weight volume --kr0 "$KR0" --skip "$SKIP" \
    --rank-scan 2:16 --archive "$out/dmd_krylov_vel" > "$out/dmd_krylov_vel.log" 2>&1 || true
./bin/run_lsa export --krylov "$DIR" --fields L0101,x_velocity,y_velocity --export "$out/export_krylov" \
    > "$out/export_krylov.log" 2>&1 || true
echo "level set: $(grep -E 'median rate|median vs' "$out/dmd_krylov_ls.log" | tr '\n' ' ')"
echo "velocity:  $(grep -E 'median rate|median vs' "$out/dmd_krylov_vel.log" | tr '\n' ' ')"
echo "[$(date +%T)] done; results in $out"
