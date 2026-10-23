#!/bin/bash
# =====================================================================
# estimate_cost.sh -- grid, steps, core-hours and wall time of a case file.
# =====================================================================
# Plans every case with bin/make_inputs (nothing is run) and converts the step count and cell count into cost, from the one measured speed of this solver: about 2 s per step for the 32 x 72 reference grid on one core with tolerances of 1e-10 (FINDINGS_FROM_REAL_RUNS.md). Cost per step is taken to scale with the number of cells, and the parallel efficiency with the given rank count is assumed (EFF, default 0.6). Treat the numbers as +-2x; their use is to choose the rank count, the wall-clock request and the order of the cases, and to see which cases dominate.
#
#     scripts/estimate_cost.sh campaign_16.txt 16
#     EFF=0.5 SEC_PER_CELL_STEP=4.3e-4 scripts/estimate_cost.sh campaign_16.txt 32
#
# The wall-clock request printed is 1.5x the estimate, rounded up to the hour, which is what submit_campaign.sh passes to Slurm. A case line with --twin-of is costed like the case it twins.
# ---------------------------------------------------------------------
set -u
cd "$(dirname "${BASH_SOURCE[0]}")/.."
CASEFILE="${1:?usage: estimate_cost.sh CASEFILE NPROCS}"
NPROCS="${2:-16}"
EFF="${EFF:-0.6}"
SPC="${SEC_PER_CELL_STEP:-$(awk 'BEGIN{print 2.0/(32*72)}')}"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

printf '%-16s %6s %10s %8s %8s %10s %9s %8s %8s\n' label kr0 grid boxes ranks steps core_h wall_h request
tot=0; worst=0
declare -A line_of
while read -r label kr0 eps cells rmax nper snaps rest; do
  [[ -z "${label:-}" || "${label:0:1}" == "#" ]] && continue
  line_of[$label]="$kr0 $eps $cells $rmax $nper $snaps ${rest:-}"
  extra=""; twin=""
  set -- ${rest:-}
  while (( $# > 0 )); do
    case "$1" in --twin-of) twin="$2"; shift 2 ;; \#*) break ;; *) extra="$extra $1"; shift ;; esac
  done
  if [[ -n "$twin" && -n "${line_of[$twin]:-}" ]]; then
    read -r kr0 eps cells rmax nper snaps rest <<< "${line_of[$twin]}"
    extra=""; set -- ${rest:-}
    while (( $# > 0 )); do case "$1" in --twin-of) shift 2 ;; \#*) break ;; *) extra="$extra $1"; shift ;; esac; done
  fi
  e="$eps"; [[ "$label" == *.base ]] && e=0.0
  out=$(./bin/make_inputs single --kr0 "$kr0" --eps "$e" --cells-per-r0 "$cells" --rmax "$rmax" \
        --n-periods "$nper" --snapshots "$snaps" --outdir "$tmp" $extra 2>&1 | sed -n '2,4p')
  nr=$(grep -oE 'grid [0-9]+x[0-9]+' <<< "$out" | head -1 | sed 's/grid //; s/x.*//')
  nz=$(grep -oE 'grid [0-9]+x[0-9]+' <<< "$out" | head -1 | sed 's/.*x//')
  steps=$(grep -oE 'steps=[0-9]+' <<< "$out" | head -1 | cut -d= -f2)
  [[ -z "$steps" ]] && { echo "$label: could not plan ($out)"; continue; }
  boxes=$(( ((nr+15)/16) * ((nz+15)/16) ))
  ranks=$(( boxes < NPROCS ? boxes : NPROCS ))
  read -r coreh wallh req < <(awk -v s="$steps" -v c="$((nr*nz))" -v spc="$SPC" -v r="$ranks" -v eff="$EFF" \
     'BEGIN{ch=s*c*spc/3600; w=ch/(r*(r>1?eff:1)); q=int(1.5*w)+1; printf "%.1f %.2f %d\n", ch, w, q}')
  printf '%-16s %6s %10s %8s %8s %10s %9s %8s %7s:00\n' "$label" "$kr0" "${nr}x${nz}" "$boxes" "$ranks" "$steps" "$coreh" "$wallh" "$req"
  tot=$(awk -v a="$tot" -v b="$coreh" 'BEGIN{print a+b}')
  worst=$(awk -v a="$worst" -v b="$req" 'BEGIN{print (b>a)?b:a}')
done < "$CASEFILE"
echo
echo "total ${tot} core-hours; longest single case needs a ${worst}:00:00 request"
echo "(WALLTIME for submit_hpc_array.sh is per task, so the longest case sets it)"
