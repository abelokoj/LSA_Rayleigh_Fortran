#!/bin/bash
# krylov_upstream.sh -- run the solver's own restarting Krylov LSA chain with the current upstream code (msussman42/amrex_implicit_interfaces, LSA_activate = 1 then 2), and re-analyse the Gram matrix it prints.
#
# Upstream (from 2026-10-15) completes the chain inside the solver: LSA_activate=1 writes the members chk<M>LSA<k>, k = 1..N; LSA_activate=2 reads them, prints the inner products "i,j,dot_product i j value", calls fort_process_LSA (PROB_CPP_PARMS.F90: eq. 67 by SVD with LAPACK, companion matrix, eigenvalues by ZGGEV, most dangerous mode = largest real eigenvalue) and writes the mode as chk<M>LSA<N+1>. This script runs that chain on the deck of our campaign and then runs krylov_test/lsa_log_spectrum.py on the log, which repeats the decomposition from the printed Gram matrix with fort_process_LSA's rule and with the rule of svd_and_diagonalize_sussman.f90 (cutoff and largest amplitude).
#
# Prerequisite: a separate worktree of the upstream solver built with USE_LIBLINEAR = TRUE (see krylov_test/README.md, "Upstream chain"), so the campaign executable is never touched:
#   git -C $AII fetch origin && git -C $AII worktree add ~/aii_up origin/master
#   cd ~/aii_up/NS_fluids_lib/BLASLAPACK/LINEAR && make && cp liblinear.a ../../
#   cd ~/aii_up/NS_fluids_lib && make -j16 USE_LIBLINEAR=TRUE LIBRARY_LOCATIONS=$PWD
#
# Usage (from $LSA):
#   SOLVER_EXE=~/aii_up/NS_fluids_lib/amr2d.gnu.MPI.ex sbatch --export=ALL krylov_test/krylov_upstream.sh [LABEL KR0 CELLS RMAX NMEMBERS TSAMPLE]
#   defaults: kup070_16 0.7 16 4.0 20 0.25
#
# Caution, checked against the upstream source of 2026-10-15 and to be confirmed with Prof. Sussman: (1) upstream rescales every member before it is re-added (level set to a maximum of 2 dx in the band, velocity to LSA_velocity_scale), so the printed Gram matrix is that of rescaled members and its companion eigenvalues need the scale factors (lsa_log_spectrum.py --scales) to be eigenvalues of the Jacobian; (2) the restart name for LSA_activate=2 and the member checkpoint names follow AMReX_AmrCore.cpp of that version and may change.
#SBATCH --job-name=lsa_krylov_up
#SBATCH --ntasks=16
#SBATCH --time=12:00:00
#SBATCH --account=genacc_q
#SBATCH --output=krylov_test/slurm-%x-%j.out
#SBATCH --error=krylov_test/slurm-%x-%j.err
#SBATCH --mail-type=END,FAIL
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
TEST="$PWD/krylov_test"
LABEL="${1:-kup070_16}"; KR0="${2:-0.7}"; CELLS="${3:-16}"; RMAX="${4:-4.0}"
NMEM="${5:-20}"; TSAMPLE="${6:-0.25}"
EXE="${SOLVER_EXE:?set SOLVER_EXE to the upstream executable built with USE_LIBLINEAR=TRUE}"
EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"
NPROCS="${SLURM_NTASKS:-16}"
SOLVER_OPTS="${SOLVER_OPTS:-mac.mac_abs_tol=1.0e-10 mac.visc_abs_tol=1.0e-10}"
unset SLURM_STEP_ID SLURM_STEPID SLURM_STEP_NUM_TASKS SLURM_STEP_NODELIST || true
set +u
if type module >/dev/null 2>&1; then module purge; module load gnu/13; module load openmpi; fi
set -u
[[ -x "$EXE" ]] || { echo "error: solver $EXE not found" >&2; exit 1; }
grep -q "fort_process_LSA" "$EXE" || { echo "error: $EXE has no fort_process_LSA (not the upstream LSA_activate=2 code)" >&2; exit 1; }

DIR="$TEST/runs/$LABEL"
mkdir -p "$DIR/decks" "$TEST/analysis/$LABEL"
./bin/make_inputs single --kr0 "$KR0" --eps 1e-3 --cells-per-r0 "$CELLS" --rmax "$RMAX" \
    --n-periods 4 --snapshots 100 --outdir "$DIR/decks" > "$DIR/plan.txt"
base=$(grep -E '/inputs\.k.*\.base$' "$DIR/plan.txt" | awk '{print $NF}' | tail -1)
cp "$base" "$DIR/inputs"
dt=$(awk -F= '/^ns.fixed_dt/{print $2+0}' "$DIR/inputs")
M=$(awk -v t="$TSAMPLE" -v d="$dt" 'BEGIN{m=int(t/d+0.5); if (m<1) m=1; print m}')
TAU=$(awk -v m="$M" -v d="$dt" 'BEGIN{printf "%.10g", m*d}')
LAUNCH=(srun --cpu-bind=none -n "$NPROCS")
COMMON=($SOLVER_OPTS amr.max_grid_size=16 amr.LSA_nsteps_krylov_subspace_method=$NMEM amr.LSA_plot_index=0)
{
  echo "--- $(date '+%Y-%m-%d %H:%M:%S') krylov_upstream.sh"
  echo "label = $LABEL  kr0 = $KR0  cells = $CELLS  rmax = $RMAX  members = $NMEM  M = $M  tau = $TAU"
  echo "solver = $EXE  md5 = $(md5sum "$EXE" | cut -c1-32)  commit = $(git -C "$(dirname "$EXE")" rev-parse --short HEAD 2>/dev/null || echo ?)"
} >> "$DIR/manifest.txt"
echo "case $LABEL: M = $M steps per member, tau = $TAU, $NMEM members"
cd "$DIR"

# 1. base state: one step without LSA, writes chk00000LSA00000
if [[ ! -f chk00000LSA00000/Header ]]; then
  "${LAUNCH[@]}" "$EXE" inputs "${COMMON[@]}" amr.LSA_activate=0 max_step=1 amr.check_int=1 amr.plot_int=1 \
      < /dev/null > base.out 2> base.err
fi
[[ -f chk00000LSA00000/Header ]] || { echo "error: chk00000LSA00000 not written; see $DIR/base.err" >&2; exit 1; }

# 2. the Krylov members (LSA_activate=1), each M steps, checkpoints chk<M>LSA<k>
last=$(printf "chk%05dLSA%05d" "$M" "$NMEM")
if [[ ! -f "$last/Header" ]]; then
  echo "[$(date +%T)] members (LSA_activate=1)"
  "${LAUNCH[@]}" "$EXE" inputs "${COMMON[@]}" amr.restart=chk00000LSA00000 amr.LSA_activate=1 \
      amr.LSA_initial_levelSteps=0 max_step=$M stop_time=1.0e6 amr.plot_int=$M amr.check_int=$M \
      < /dev/null > lsa1.out 2> lsa1.err || echo "WARNING: LSA_activate=1 ended with an error; see $DIR/lsa1.err"
fi
ls -d chk*LSA* 2>/dev/null | sed 's/^/  /'

# 3. inner products, fort_process_LSA, most dangerous mode (LSA_activate=2)
echo "[$(date +%T)] Gram matrix and most dangerous mode (LSA_activate=2)"
"${LAUNCH[@]}" "$EXE" inputs "${COMMON[@]}" amr.restart=chk00000LSA00000 amr.LSA_activate=2 \
    amr.LSA_initial_levelSteps=0 max_step=$M stop_time=1.0e6 amr.plot_int=$M amr.check_int=$M \
    < /dev/null > lsa2.out 2> lsa2.err || echo "WARNING: LSA_activate=2 ended with an error; see $DIR/lsa2.err"
cd - > /dev/null

# 4. re-analysis of the printed Gram matrix
out="$TEST/analysis/$LABEL"
cp "$DIR/lsa2.out" "$out/" 2>/dev/null || true
grep -E "i,j,dot_product|eigenvalue|max_index|X_dangerous" "$DIR/lsa2.out" > "$out/lsa2_summary.txt" 2>/dev/null || true
if command -v python3 >/dev/null 2>&1; then
  python3 krylov_test/lsa_log_spectrum.py "$DIR/lsa2.out" --tau "$TAU" --kr0 "$KR0" > "$out/spectrum.txt" 2>&1 || true
  cat "$out/spectrum.txt"
fi
echo "[$(date +%T)] done; results in $out"
