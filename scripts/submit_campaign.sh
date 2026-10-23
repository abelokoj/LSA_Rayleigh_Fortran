#!/bin/bash
# =====================================================================
# submit_campaign.sh -- submit the runs of the dissertation in one command (Slurm).
# =====================================================================
# Submits each case file as its own Slurm array, with the rank count and wall-clock limit it was run with, then one small post-processing job that waits for all of them (--dependency=afterany) and runs scripts/analyse_campaign.sh. Run on a login node, from the solver-tree copy of this repository ($LSA):
#
#     scripts/submit_campaign.sh                              # every case file of the dissertation (table below)
#     scripts/submit_campaign.sh campaign_16.txt rerun_meas.txt   # only these files
#
#   case file               ranks  wall clock  dissertation
#   campaign_16_wide.txt      32   13 h        sweep k r0 = 0.2-0.4 (Tables 5.3, 5.6, 6.1)
#   campaign_16.txt           16    6 h        sweep k r0 = 0.5-0.9, linearity, time step, wall, gas, stable mode (Tables 5.3-5.6, 6.1, 6.2)
#   rerun_meas.txt            16   30 h        mesh study at 8, 16, 32 cells (Table 5.4)
#   campaign_next_16.txt      16   12 h        cut-off, confinement, viscous ambient (Chapter 6)
#   campaign_amp08.txt         8    4 h        linearity at 8 cells (Table 6.2)
#   campaign_amp32.txt        16   36 h        linearity at 32 cells (Table 6.2); packed on one node
#   campaign_check.txt        16   48 h        wall at 1.5 R0 on 32 cells; viscous ambient with the wall at 6 and 8 R0 (Chapter 6)
#
# The wall-clock limits are about 1.5 times the longest task of each file; scripts/estimate_cost.sh prints the estimate for any file. A time-out is not fatal: cases whose record exists in runs/records are skipped and unfinished cases resume from their newest checkpoint, so resubmitting the same command finishes the work. The September 2026 runs (full_analysis.txt, 8 cells, eps = 2e-2) are not needed for any table; to repeat them, pass full_analysis.txt explicitly (8 ranks, 12 h).
#
# Settings (environment): ACCOUNT (default genacc_q, the FSU RCC general-access account; set your own), THROTTLE (8 tasks at once per array), POST_TIME (04:00:00).
# ---------------------------------------------------------------------
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

if [[ -n "${SLURM_JOB_ID:-}" && -z "${LSA_ALLOW_NESTED:-}" ]]; then
  echo "error: run this on a login node, not inside a Slurm job (see submit_hpc_array.sh)" >&2
  exit 1
fi

export ACCOUNT="${ACCOUNT:-genacc_q}"
export THROTTLE="${THROTTLE:-8}"
ids=()

# ranks, wall clock and packing of each case file, as run for the dissertation
settings () {
  case "$1" in
    campaign_16_wide.txt) echo "32 13:00:00 0" ;;
    campaign_16.txt)      echo "16 06:00:00 0" ;;
    rerun_meas.txt)       echo "16 30:00:00 0" ;;
    campaign_next_16.txt) echo "16 12:00:00 0" ;;
    campaign_amp08.txt)   echo "8 04:00:00 0" ;;
    campaign_amp32.txt)   echo "16 36:00:00 1" ;;
    campaign_check.txt)   echo "16 48:00:00 0" ;;
    full_analysis.txt)    echo "8 12:00:00 0" ;;
    *)                    echo "16 12:00:00 0" ;;
  esac
}

submit () {   # casefile
  local f="$1" n w p out id
  [[ -f "$f" ]] || { echo "skip: $f not found"; return; }
  read -r n w p <<< "$(settings "$f")"
  echo ">>> $f on $n ranks, $w per task"
  if [[ "$p" == 1 ]]; then
    out=$(PACK=1 WALLTIME="$w" scripts/submit_hpc_array.sh "$f" "$n")
  else
    out=$(WALLTIME="$w" scripts/submit_hpc_array.sh "$f" "$n")
  fi
  echo "$out" | sed 's/^/    /'
  id=$(echo "$out" | awk '/^array job/{print $3}' | tail -1)
  [[ -n "$id" ]] && ids+=("$id")
}

files=("$@")
(( ${#files[@]} )) || files=(campaign_16_wide.txt campaign_16.txt rerun_meas.txt campaign_next_16.txt campaign_amp08.txt campaign_amp32.txt campaign_check.txt)
for f in "${files[@]}"; do submit "$f"; done

if (( ${#ids[@]} == 0 )); then
  echo "nothing submitted"; exit 1
fi
dep="afterany:$(IFS=:; echo "${ids[*]}")"
mkdir -p runs/array
post=$(sbatch --parsable --dependency="$dep" --job-name=lsa_post \
        --account="$ACCOUNT" --ntasks=1 --cpus-per-task=4 --mem=16G \
        --time="${POST_TIME:-04:00:00}" \
        --output="runs/array/post-%j.out" --error="runs/array/post-%j.err" \
        --mail-type=END,FAIL "$HERE/scripts/_post_task.sh")
echo ">>> post-processing job ${post%%;*}, waits for: ${ids[*]}"
echo "watch with: squeue --me ; results in runs/analysis/summary.csv"
