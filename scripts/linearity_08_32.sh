#!/bin/bash
# =====================================================================
# linearity_08_32.sh -- proportionality test at 8 and 32 cells per radius (Table 6.2 of the dissertation).
# =====================================================================
# The analysis job (scripts/analyse_campaign.sh, step 4b) writes the 16-cell report runs/analysis/linearity_16.txt. This script writes the 8- and 32-cell reports from the runs of campaign_amp08.txt and campaign_amp32.txt, with the same options: the full-state norm of the twin difference for eps = 2.5e-3, 5e-3 and 1e-2, compared up to the end of the admissible window of the largest amplitude (column t_hi of its record fit). It is the script that produced results/followup_2026-10-22/analysis/linearity_08.txt and linearity_32.txt.
#
# Run from the solver-tree copy of this repository ($LSA) after the analysis job has written runs/analysis/records/amp08_100.csv and amp32_100.csv:
#
#     sbatch -A <account> -n 1 -c 2 --mem=8G -t 02:00:00 -o runs/array/linearity-%j.out scripts/linearity_08_32.sh 08 32
#
# or, interactively on a compute node, scripts/linearity_08_32.sh 08 32.
# ---------------------------------------------------------------------
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v module >/dev/null 2>&1; then
  set +u; module purge; module load gnu/13; module load openmpi; set -u
fi
[[ -x bin/run_lsa ]] || { echo "bin/run_lsa missing; run make first" >&2; exit 1; }
(( $# )) || set -- 08 32
for g in "$@"; do
  rec="runs/analysis/records/amp${g}_100.csv"
  [[ -f "$rec" ]] || { echo "skip $g: $rec not found (run the analysis job first)"; continue; }
  tmax=$(awk -F, 'NR==2{print $4}' "$rec")
  ./bin/run_lsa linearity --base "runs/amp${g}.base" --kr0 0.7 --fields x_velocity,y_velocity --t-max "$tmax" \
    --case "2.5e-3=runs/amp${g}_025" --case "5e-3=runs/amp${g}_050" --case "1e-2=runs/amp${g}_100" \
    > "runs/analysis/linearity_${g}.txt"
  echo "== linearity, $g cells per radius"; grep -E "verdict|spread" "runs/analysis/linearity_${g}.txt"
done
