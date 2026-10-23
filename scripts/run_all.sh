#!/usr/bin/env bash
# =====================================================================
# run_all.sh -- advance every calculation listed in a sweep file, then extract
# the interface mode amplitude from each.
#
# Run from amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Fortran:
#
#     ./scripts/run_all.sh                          # full_analysis.txt, 1 rank
#     ./scripts/run_all.sh full_analysis.txt 16   # 16 ranks
#     ./scripts/run_all.sh my_subset.txt 8 ../amr2d.gnu.MPI.ex
#

# The launcher is selected by the LAUNCHER environment variable, which defaults to mpirun. Under Slurm it should be set to srun, which the FSU Research Computing Center documentation recommends in preference to mpirun because srun integrates with the scheduler and inherits the allocation, so the rank count and the placement come from the job rather than from a command-line argument that may disagree with it. Setting it also avoids the singleton failure mode, in which a launcher that cannot reach a process manager starts independent single-rank copies instead of one job:
#
#     LAUNCHER=srun ./scripts/run_all.sh full_analysis.txt 16
#
# The script is designed to be started once and left alone. It is also safe to
# interrupt and restart: a case whose amplitude record already exists is
# skipped, and a case that is partly advanced resumes from its most recent
# checkpoint. Progress is appended to runs/run_all.log, and a machine
# readable summary is written to runs/summary.tsv as each case finishes.
#
# Cases are attempted in file order and a failure does not stop the sweep, so
# that one unstable configuration cannot cost the whole run. Failures are
# recorded in the log and marked in the summary.
# =====================================================================
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

CASEFILE="${1:-full_analysis.txt}"
NPROCS="${2:-1}"
# The solver: the third argument, else $SOLVER_EXE, else ../amr2d.gnu.MPI.ex (this directory placed inside NS_fluids_lib). SOLVER_EXE lets a clone of the analysis repository live anywhere, e.g. export SOLVER_EXE=$AII/NS_fluids_lib/amr2d.gnu.MPI.ex.
EXE="${3:-${SOLVER_EXE:-../amr2d.gnu.MPI.ex}}"

# Launcher. Set LAUNCHER=srun under Slurm; see the header of this script.
LAUNCHER="${LAUNCHER:-mpirun}"
if [[ "$LAUNCHER" == mpirun* ]]; then
  LAUNCHER_OPTS="${LAUNCHER_OPTS:---allow-run-as-root --oversubscribe}"
else
  # --cpu-bind=none is required, not cosmetic. Slurm satisfies a request for N tasks from whatever cores are free, which on a busy cluster means a few cores on each of several nodes, in fragmented masks such as 0x140000 or 0x000001. srun's default binding cannot be mapped onto those and every launch is refused with "Unable to satisfy cpu bind request", so the solver never starts. Leaving placement to the operating system costs a little locality and avoids the failure entirely.
  LAUNCHER_OPTS="${LAUNCHER_OPTS:---cpu-bind=none}"
fi

RUNROOT="runs"
# The deck directory may be redirected. Generated deck names encode the wavenumber and amplitude but not the resolution or the domain size, so two cases that differ only in those would be written to the same filename; concurrent jobs must therefore each be given their own directory. See scripts/_array_task.sh.
DECKDIR="${DECKDIR:-$RUNROOT/decks}"
RECDIR="$RUNROOT/records"
# Figure data. Each case gets its own bundle of interface profiles and field slices, written straight after its amplitude record, so that every figure of the study can be redrawn later from plain text. The plotfiles usually stay on the machine that ran the sweep and amount to tens of gigabytes; a bundle is a few hundred kilobytes and travels in one copy. Set FIELDS empty to skip the field slices, which are the larger part.
BUNDLEDIR="${BUNDLEDIR:-$RUNROOT/bundles}"
FIELDS="${FIELDS:-y_velocity}"
# The log and the summary may be redirected so that concurrent jobs, such as the tasks of a Slurm array, each write their own and are merged afterwards rather than appending to one file at once.
LOG="${LOG:-$RUNROOT/run_all.log}"
SUMMARY="${SUMMARY:-$RUNROOT/summary.tsv}"

# Solver tolerances. These must sit well below the smallest perturbation
# amplitude in the sweep, or the smallest case measures solver residual rather
# than physics.
# Overridable, in the manner of CHECK_INT, so that a tolerance or a stop time may be varied for one sweep without editing this file.
SOLVER_OPTS="${SOLVER_OPTS:-mac.mac_abs_tol=1.0e-10 mac.visc_abs_tol=1.0e-10}"

# Checkpoint interval. This must be smaller than the number of steps a single
# scheduler allocation can complete, or no checkpoint is ever written and each
# restart silently repeats from the beginning.
CHECK_INT="${CHECK_INT:-50}"

# Checkpoints kept per case. AMReX never deletes an old checkpoint, and at CHECK_INT=50 a long run leaves dozens of full copies of its state; the October campaign filled the home quota this way, after which the solver could no longer write checkpoints or plotfiles and every new run aborted. Only the newest KEEP_CHK complete checkpoints are kept while a case runs. Once its record is complete, the newest KEEP_CHK_DONE are kept as a reference: a finished case is never resumed by the sweep, but they allow it to be restarted by hand to inspect or extend the end of the run. KEEP_CHK_DONE=0 removes them all.
KEEP_CHK="${KEEP_CHK:-2}"
KEEP_CHK_DONE="${KEEP_CHK_DONE:-2}"
PRUNE_EVERY="${PRUNE_EVERY:-120}"     # seconds between prunes during a run

mkdir -p "$RUNROOT" "$DECKDIR" "$RECDIR" "$BUNDLEDIR"

log () { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# Deletes all but the newest KEEP_CHK complete checkpoints of one case. A checkpoint counts as complete only if its Header exists (AMReX writes it last), so the one being written at that moment is never touched.
prune_checkpoints () {
  local dir=$1 keep=${2:-$KEEP_CHK}
  ( cd "$dir" 2>/dev/null || exit 0
    for c in chk[0-9]*; do
      [[ "$c" =~ ^chk[0-9]+$ && -d "$c" && -f "$c/Header" ]] && printf '%s\n' "$c"
    done | sort -V | head -n -"$keep" | while read -r c; do rm -rf -- "$c"; done )
}

if [[ ! -f "$CASEFILE" ]]; then
  log "ERROR: case file '$CASEFILE' not found"
  exit 1
fi
if [[ ! -x "$EXE" ]]; then
  log "ERROR: solver '$EXE' is not executable; build it with"
  log "       cd $(cd .. && pwd) && make -j16"
  exit 1
fi
# Resolved to an absolute path because each case is advanced from its own run
# directory, where a path relative to this one would no longer refer to it.
EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"
if [[ ! -x ./bin/make_inputs || ! -x ./bin/dump_amplitude || ! -x ./bin/dump_profiles ]]; then
  log "building the analysis tools"
  make >/dev/null || { log "ERROR: make failed"; exit 1; }
fi

[[ -f "$SUMMARY" ]] || printf 'label\tkr0\teps\tcells\trmax\tstatus\tsteps\tt_final\trecord\n' > "$SUMMARY"
# The summary is append-only and survives restarts, so the failure count at the end must consider only the rows this invocation adds. Counting the whole file kept a failure from an earlier attempt alive after the case had been rerun successfully, and every later job then exited non-zero.
SUMMARY_START=$(wc -l < "$SUMMARY")

if (( NPROCS > 1 )); then
  if [[ "$LAUNCHER" == mpirun* ]] && \
     mpirun -n 2 hostname 2>&1 | grep -q 'singletons will be started'; then
    log "ERROR: mpirun cannot reach a process manager and would start $NPROCS independent single-rank copies of the solver, all writing into the same directory. Re-run with 1 rank or repair the launcher; see docs/solver_notes.md and REPRODUCE.md, Section 3.3."
    exit 1
  fi
fi

log "sweep '$CASEFILE' on $NPROCS rank(s), solver '$EXE'"

# ---------------------------------------------------------------------
# Advance one case to completion, restarting as often as necessary.
#
# The solver is invoked repeatedly rather than once because a single
# invocation may be cut short by a scheduler limit; each invocation resumes
# from the latest checkpoint. The loop terminates when the step count stops
# advancing, which covers both normal completion and a stalled calculation.
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# Domain decomposition. AMReX aborts with "cannot have a dormant/idle processor on any level" when a rank is given no box, so the number of boxes must be at least the number of ranks. The box count follows from amr.n_cell and amr.max_grid_size, both of which are known once the deck exists, so the largest max_grid_size that supplies enough boxes is chosen here and passed on the command line. Where even the smallest sensible value cannot supply enough, the rank count is reduced instead, since running a case on fewer ranks is always preferable to not running it.
#
# Larger boxes are preferred because each carries a ghost region whose cost grows as the surface-to-volume ratio worsens, so the search runs from 32 downwards rather than simply taking the smallest.
# ---------------------------------------------------------------------
decompose () {
  local deck="$1" want="$2" nx nz mgs boxes cgmin smallest smallest_boxes
  nx=$(awk '/^amr.n_cell/{print $3; exit}' "$deck")
  nz=$(awk '/^amr.n_cell/{print $4; exit}' "$deck")
  if [[ -z "$nx" || -z "$nz" ]]; then echo "32 $want"; return; fi

  # The solver requires amr.max_grid_size >= cg.min_max_grid_size and aborts otherwise, so the floor is read from the deck rather than assumed. Lowering cg.min_max_grid_size to admit smaller boxes was rejected: it alters the bottom of the multigrid solve, and every case of the sweep should be advanced under the same solver settings so that the measured rates remain comparable. Reducing the rank count instead changes only the decomposition.
  cgmin=$(awk '/^cg.min_max_grid_size/{print $3; exit}' "$deck")
  [[ -z "$cgmin" ]] && cgmin=16

  smallest=""; smallest_boxes=0
  for mgs in 32 16 8; do
    (( mgs < cgmin )) && continue
    boxes=$(( ( (nx + mgs - 1) / mgs ) * ( (nz + mgs - 1) / mgs ) ))
    smallest=$mgs; smallest_boxes=$boxes
    if (( boxes >= want )); then echo "$mgs $want"; return; fi
  done

  # Not enough boxes at any permitted size: run on as many ranks as there are boxes.
  [[ -z "$smallest" ]] && { echo "$cgmin 1"; return; }
  (( smallest_boxes < 1 )) && smallest_boxes=1
  echo "$smallest $smallest_boxes"
}

# ---------------------------------------------------------------------
# Run record. The September 2026 evidence report could not say which solver build, rank count or deck produced each archived record (its section 11.4). Every invocation now appends a block to runs/<label>/manifest.txt recording exactly that, so that a result can be traced and a twin can be checked against its perturbed partner.
# ---------------------------------------------------------------------
write_manifest () {
  local dir="$1" label="$2" mgs="$3" ranks="$4" commit="unknown" exe_sum="unknown" deck_sum="unknown"
  commit=$(git -C "$(dirname "$EXE")" rev-parse HEAD 2>/dev/null || echo unknown)
  exe_sum=$(md5sum "$EXE" 2>/dev/null | awk '{print $1}')
  deck_sum=$(md5sum "$dir/inputs" 2>/dev/null | awk '{print $1}')
  {
    echo "--- $(date '+%Y-%m-%d %H:%M:%S')"
    echo "label           = $label"
    echo "case_line       = $CASE_LINE"
    echo "host            = $(hostname 2>/dev/null)"
    echo "slurm_job       = ${SLURM_JOB_ID:-none} array=${SLURM_ARRAY_JOB_ID:-none}_${SLURM_ARRAY_TASK_ID:-none}"
    echo "launcher        = $LAUNCHER $LAUNCHER_OPTS"
    echo "mpi_ranks       = $ranks"
    echo "max_grid_size   = $mgs"
    echo "solver          = $EXE"
    echo "solver_md5      = ${exe_sum:-unknown}"
    echo "solver_commit   = $commit"
    echo "deck_md5        = ${deck_sum:-unknown}"
    echo "solver_opts     = $SOLVER_OPTS amr.check_int=$CHECK_INT"
    if type module >/dev/null 2>&1; then
      echo "modules         = $(module -t list 2>&1 | tr '\n' ' ')"
    else
      echo "modules         = (module command not available)"
    fi
  } >> "$dir/manifest.txt" 2>/dev/null
}

advance () {
  local dir="$1" deck="$2" label="$3"
  mkdir -p "$dir"
  [[ -f "$dir/inputs" ]] || cp "$deck" "$dir/inputs"

  local mgs ranks
  read -r mgs ranks <<< "$(decompose "$dir/inputs" "$NPROCS")"
  if (( ranks < NPROCS )); then
    log "  $label: only $ranks box(es) available, running on $ranks rank(s) rather than $NPROCS" >&2
  fi

  write_manifest "$dir" "$label" "$mgs" "$ranks"

  local prev=-1 now=0 attempt=0
  while (( attempt < 200 )); do
    attempt=$(( attempt + 1 ))
    local last restart=""
    # The newest complete checkpoint. sort -V orders the step numbers numerically, so chk100000 follows chk99950 once the step count gains a digit, where a plain sort would put it first and restart from an older state. A checkpoint is accepted only if its Header exists, because AMReX writes the Header last: a job stopped at its wall-clock limit or pre-empted on backfill2 while writing a checkpoint leaves a directory without one, and restarting from it would abort or read a half-written state.
    # Only names that are chk followed by digits are considered, so that a renamed or temporary checkpoint directory (chk00100.old, chk00100.temp) can never be chosen.
    last=$(cd "$dir" && for c in chk[0-9]*; do
             [[ "$c" =~ ^chk[0-9]+$ && -d "$c" && -f "$c/Header" ]] \
               && printf '%s\n' "$c"
           done | sort -V | tail -1)
    [[ -n "$last" ]] && restart="amr.restart=$last"

    # stdin is redirected because mpirun would otherwise consume the case
    # file that the main loop is reading, so that only the first case runs.
    # Detached from stdout: this function's stdout is captured by the caller, and a sleeping pruner holding it open would delay every case by up to PRUNE_EVERY seconds.
    ( while sleep "$PRUNE_EVERY"; do prune_checkpoints "$dir"; done ) </dev/null >/dev/null 2>&1 &
    local pruner=$!
    ( cd "$dir" && $LAUNCHER $LAUNCHER_OPTS -n "$ranks" \
        "$EXE" inputs $SOLVER_OPTS amr.max_grid_size=$mgs \
        amr.check_int=$CHECK_INT $restart \
        < /dev/null >> solver.out 2>> solver.err )
    kill "$pruner" 2>/dev/null; wait "$pruner" 2>/dev/null
    prune_checkpoints "$dir"

    now=$(grep -h '^STEP' "$dir/solver.out" 2>/dev/null | tail -1 \
          | awk '{print $3}')
    now=${now:-0}
    if [[ "$now" == "$prev" ]]; then
      # Zero steps on the first attempt means the solver never started rather than that it finished. Reporting the reason here saves reading two log files to discover that the launcher, and not the analysis, is at fault.
      if (( attempt <= 2 )) && [[ "$now" == "0" ]]; then
        log "  $label: the solver produced no time steps; see $dir/solver.err" >&2
        if [[ -s "$dir/solver.err" ]]; then
          # Skipped: the X11 and launcher chatter that heads the file and says nothing about the failure.
          log "  $label: $(grep -vE 'No protocol specified|^ *$' "$dir/solver.err" \
                 | grep -iE 'abort|error|fail|expecting|not found' \
                 | head -3 | tr '\n' ' ')" >&2
        fi
      fi
      break                     # no progress: finished, or stalled
    fi
    prev="$now"
    # To stderr: this function's stdout is captured by the caller.
    log "  $label: step $now" >&2
  done
  echo "$now"
}

# ---------------------------------------------------------------------
# Main loop over the case file.
# ---------------------------------------------------------------------
while read -r label kr0 eps cells rmax nper snaps rest; do
  CASE_LINE="$label $kr0 $eps $cells $rmax $nper $snaps ${rest:-}"
  # Skip blank lines and comments. The test is on the first field, so an
  # indented comment is handled as well as one in the first column.
  [[ -z "${label:-}" ]] && continue
  [[ "${label:0:1}" == "#" ]] && continue

  record="$RECDIR/${label}.dat"
  if [[ -s "$record" ]]; then
    log "$label: already complete, skipping"
    continue
  fi

  log "$label: kr0=$kr0 eps=$eps cells=$cells rmax=$rmax periods=$nper ${rest:-}"

  # Extra columns after the seventh are options. --twin-of LABEL is handled here; everything else is passed to make_inputs unchanged (for example --safety 0.125, --rho-g 1.225e-4, --mu-g 2.6e-5, --n-waves 8, --eps-cells 0.08), so that a case can vary any generator setting without a new column.
  extra=""; twin_of=""
  set -- ${rest:-}
  while (( $# > 0 )); do
    case "$1" in
      --twin-of) twin_of="${2:-}"; shift 2 ;;
      \#*) break ;;                       # trailing comment
      *) extra="$extra $1"; shift ;;
    esac
  done

  # A label ending in .base denotes the unperturbed twin, for which the
  # amplitude field is meaningless and is forced to zero.
  eps_use="$eps"
  case "$label" in
    *.base) eps_use="0.0" ;;
  esac

  if [[ -n "$twin_of" ]]; then
    # The twin of an existing perturbed run is made from that run's own deck, with ns.radblob set to zero and nothing else changed. This is the only way to pair a new twin with a perturbed run made earlier, possibly by an older generator: regenerating the deck could differ in a setting that the case line does not record. The twin must also run on the same number of ranks as its partner; the partner's manifest is checked when it has one.
    src="$RUNROOT/$twin_of/inputs"
    if [[ ! -f "$src" ]]; then
      log "$label: FAILED, --twin-of $twin_of but $src does not exist"
      printf '%s\t%s\t%s\t%s\t%s\tdeck_missing\t-\t-\t-\n' \
        "$label" "$kr0" "$eps" "$cells" "$rmax" >> "$SUMMARY"
      continue
    fi
    mkdir -p "$RUNROOT/$label"
    deck="$RUNROOT/$label/inputs.twin"
    sed -e 's/^\(ns\.radblob[[:space:]]*=[[:space:]]*\)[^#]*/\10.0  /' "$src" > "$deck"
    ndiff=$(diff <(grep -v '^[[:space:]]*#' "$src" | sed 's/[[:space:]]*#.*//') \
                 <(grep -v '^[[:space:]]*#' "$deck" | sed 's/[[:space:]]*#.*//') | grep -c '^<' || true)
    if [[ "$ndiff" != "1" ]]; then
      log "$label: FAILED, twin deck differs from $src in $ndiff functional lines (expected exactly 1, ns.radblob)"
      continue
    fi
    pr=$(awk -F'= *' '/^mpi_ranks/{r=$2} END{print r}' "$RUNROOT/$twin_of/manifest.txt" 2>/dev/null || true)
    if [[ -n "$pr" ]]; then
      want=$(decompose "$deck" "$NPROCS" | awk '{print $2}')
      if [[ "$pr" != "$want" ]]; then
        log "$label: WARNING, partner $twin_of ran on $pr rank(s) but this job would use $want; resubmit with $pr ranks or the twin subtraction is not exact"
      fi
    else
      log "$label: note, $twin_of has no manifest; make sure this job uses the same rank count as it did (the rank count is in the partner's manifest or job log)"
    fi
  else
    if ! mi_out=$(./bin/make_inputs single --kr0 "$kr0" --eps "$eps_use" \
          --cells-per-r0 "$cells" --rmax "$rmax" --n-periods "$nper" \
          --snapshots "$snaps" --outdir "$DECKDIR" $extra 2>&1); then
      log "$label: FAILED to generate deck: $(printf '%s' "$mi_out" | grep -i error | head -2 | tr '\n' ' ')"
      printf '%s\t%s\t%s\t%s\t%s\tdeck_failed\t-\t-\t-\n' \
        "$label" "$kr0" "$eps" "$cells" "$rmax" >> "$SUMMARY"
      continue
    fi
    # The deck paths are read from make_inputs' own output (it prints the twin deck, then the perturbed one), rather than reconstructed here from the case line, which goes wrong as soon as the amplitude is set by --eps-cells.
    case "$label" in
      *.base) deck=$(printf '%s\n' "$mi_out" | grep -E '/inputs\.k[^ ]*\.base$' | head -1) ;;
      *)      deck=$(printf '%s\n' "$mi_out" | grep -E '/inputs\.k[^ ]*\.eps[^ ]*$' | head -1) ;;
    esac
  fi
  if [[ -z "$deck" || ! -f "$deck" ]]; then
    log "$label: FAILED, expected deck '$deck' was not generated"
    printf '%s\t%s\t%s\t%s\t%s\tdeck_missing\t-\t-\t-\n' \
      "$label" "$kr0" "$eps" "$cells" "$rmax" >> "$SUMMARY"
    continue
  fi

  steps=$(advance "$RUNROOT/$label" "$deck" "$label")

  # The wavelength is one axial box length and is read back from the deck so
  # that it cannot drift from the value the solver actually used.
  lambda=$(awk '/^ns.yblob/{print $3; exit}' "$RUNROOT/$label/inputs")

  rows=0
  # A broadband run (ns.radblob3 = M > 0 harmonics of the box, deck ns.yblob = box length) also gets one record per harmonic, in $RECDIR/<label>_harmonics/hNN.dat.
  nharm=$(awk '/^ns.radblob3/{print int($3); exit}' "$RUNROOT/$label/inputs")
  if ./bin/dump_amplitude "$RUNROOT/$label" "$lambda" "$record" ${nharm:+$( (( nharm > 1 )) && echo "$nharm")} >/dev/null 2>&1
  then
    # grep -c prints 0 and also exits non-zero when nothing matches, so the count is taken first and defaulted afterwards rather than through '|| echo 0', which would append a second line and make the arithmetic below malformed.
    rows=$(grep -cve '^[[:space:]]*#' -e '^[[:space:]]*$' "$record" 2>/dev/null || true)
    rows=${rows:-0}
  fi

  # A record with fewer than two rows cannot yield a growth rate, and the extractor exits zero when it writes a header and nothing else. Judging the case on the number of rows rather than on that exit status is what stops a run in which the solver never advanced from being recorded as a success.
  # The advance loop stops when the step count stops changing, which is true of a finished run and equally of one that aborts at the same step on every restart. Completion is therefore judged against the deck: the run is complete when it reached max_step, or when its last snapshot lies within one plot interval of stop_time. Without this, a crashed case with a few rows of record was marked ok and, since its record existed, skipped on every later restart.
  complete=0
  if (( rows >= 2 )); then
    tfin=$(tail -1 "$record" | awk '{print $1}')
    complete=$(awk -v f="$RUNROOT/$label/inputs" -v n="$steps" -v t="$tfin" '
      BEGIN {
        while ((getline line < f) > 0) {
          split(line, w, /[ \t=]+/)
          if (w[1] == "max_step")      ms = w[2]
          if (w[1] == "stop_time")     st = w[2]
          if (w[1] == "amr.plot_int")  pi = w[2]
          if (w[1] == "ns.fixed_dt")   dt = w[2]
        }
        ok = 0
        if (ms != "" && n + 0 >= ms + 0) ok = 1
        if (st != "" && t + 1.01 * (pi + 0) * (dt + 0) >= st + 0) ok = 1
        print ok
      }')
  fi

  if (( rows >= 2 && complete == 0 )); then
    # Kept beside the record path for diagnosis, but not under the record name, so that a restart of the sweep resumes this case from its newest checkpoint instead of skipping it. The suffix follows .dat, so that the plotting scripts, which read records/*.dat, never take a partial record for a finished one.
    mv -f "$record" "${record}.partial"
    log "$label: INCOMPLETE, stopped at step $steps, t=$tfin, short of the deck's stop_time/max_step; partial record kept as ${record}.partial; see $RUNROOT/$label/solver.err"
    printf '%s\t%s\t%s\t%s\t%s\tincomplete\t%s\t%s\t-\n' \
      "$label" "$kr0" "$eps" "$cells" "$rmax" "$steps" "$tfin" >> "$SUMMARY"
  elif (( rows >= 2 )); then
    rm -f "${record}.partial"            # superseded by the finished record
    prune_checkpoints "$RUNROOT/$label" "$KEEP_CHK_DONE"   # kept for reference only
    log "$label: complete, $steps steps, $rows snapshots, t=$tfin, record $record"

    # The figure data is exported here, while the plotfiles are still present, because this is the only moment at which they are guaranteed to exist. A failure to export is logged and does not fail the case: the growth rate is already secured by the amplitude record, and losing the interface figure of one case is not a reason to mark a completed calculation as failed.
    casebundle="$BUNDLEDIR/$label"
    mkdir -p "$casebundle"
    cp "$record" "$casebundle/amplitude.dat"
    if ! ./bin/dump_profiles "$RUNROOT/$label" "$lambda" "$casebundle" ${FIELDS%% *} >/dev/null 2>&1; then
      log "$label: figure data export failed; amplitude record is unaffected"
    fi
    printf '%s\t%s\t%s\t%s\t%s\tok\t%s\t%s\t%s\n' \
      "$label" "$kr0" "$eps" "$cells" "$rmax" "$steps" "$tfin" "$record" \
      >> "$SUMMARY"
  else
    rm -f "$record"
    log "$label: no usable amplitude record ($rows row(s) from $steps step(s)); see $RUNROOT/$label/solver.err"
    printf '%s\t%s\t%s\t%s\t%s\tno_record\t%s\t-\t-\n' \
      "$label" "$kr0" "$eps" "$cells" "$rmax" "$steps" >> "$SUMMARY"
  fi
done < "$CASEFILE"

failed=$(awk -F'\t' -v s="$SUMMARY_START" \
           'NR>s && NR>1 && $6 != "ok" { n++ } END { print n+0 }' "$SUMMARY")
if (( failed > 0 )); then
  log "sweep finished with $failed failed case(s); summary in $SUMMARY"
else
  log "sweep finished; summary in $SUMMARY"
fi
echo
if command -v column >/dev/null 2>&1; then
  column -t -s $'\t' "$SUMMARY"
else
  awk -F'\t' '{printf "%-14s %-6s %-8s %-6s %-6s %-12s %-8s %-10s\n", \
    $1, $2, $3, $4, $5, $6, $7, $8}' "$SUMMARY"
fi

# A non-zero status when any case failed, so that a scheduler records the job as failed rather than as completed. A sweep that produced no usable record is not a success, and reporting it as one is how a broken run goes unnoticed.
exit $(( failed > 0 ? 1 : 0 ))
