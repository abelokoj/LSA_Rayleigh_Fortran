#!/bin/bash
# run_studies.sh -- follow-up studies on the finished campaign, from the data the analysis job already left in runs/analysis (no solver runs).
#
#   band    narrow-band level-set DMD with bands of 1, 2, 3 and 4 cells around r = r0 (is the ~0.9 % level-set bias set by the band?)
#   stride  sensitivity of the DMD rate to the sampling interval: strides 1, 2, 4, 8, 16 of the stored snapshots, for the velocity and for the level-set band
#   drift   drift of the unperturbed twin from its initial state (run_lsa drift), for the base-state section of the verification chapter; needs the twin's plotfiles
#   modes   DMD archives with the leading modes kept (for the eigenfunction check, plotting/plot_mode_shape_archive.py)
#   all     all of the above
#
# Usage (from $LSA, on a compute node or as a batch job):
#   scripts/run_studies.sh band   [labels...]     default labels: sw050_16 sw070_16 meas32
#   scripts/run_studies.sh stride [labels...]     default: the sweep sw020_16 .. sw090_16
#   scripts/run_studies.sh drift  [labels...]     default: sw070_16 (reads runs/<label>.base)
#   scripts/run_studies.sh modes  [labels...]     default: sw040_16 sw050_16 sw060_16 sw070_16 meas32
#   sbatch -A genacc_q -n 1 -c 4 --mem=16G -t 04:00:00 -o runs/array/studies-%j.out scripts/run_studies.sh all
#
# Inputs: runs/analysis/<label>/export (velocity) and export_ls (level set), written by analyse_campaign.sh, and runs/analysis/records/<label>.csv (the interface window). Outputs: runs/analysis/studies/ (one table per study) and per-case logs under runs/analysis/<label>/studies/. Every step is skipped when its output exists; FORCE=1 redoes it.
# Settings: RUNROOT (default runs), ANALYSIS (default $RUNROOT/analysis).

set -u
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RUNROOT="${RUNROOT:-runs}"
ANALYSIS="${ANALYSIS:-$RUNROOT/analysis}"
FORCE="${FORCE:-0}"
STUDIES="$ANALYSIS/studies"
mkdir -p "$STUDIES"
[[ -x bin/run_lsa ]] || { echo "bin/run_lsa missing; run make first" >&2; exit 1; }
need () { [[ "$FORCE" == 1 || ! -e "$1" ]]; }

what="${1:-all}"; shift || true
labels=("$@")

kr0_of () {
  local k
  k=$(awk '/^#[ \t]*k r0[ \t]*=/{print $NF; exit}' "$RUNROOT/$1/inputs" 2>/dev/null)
  [[ -z "$k" ]] && k=$(awk -F, 'NR==2{print $2+0}' "$ANALYSIS/records/$1.csv" 2>/dev/null)
  echo "${k:-0.7}"
}

# Snapshots of an export inside the interface window [t_lo, t_hi] of the record fit: prints "skip maxs".
window_of () {
  local label="$1" exp="$2" tlo thi
  tlo=$(awk -F, 'NR==2{print $3}' "$ANALYSIS/records/$label.csv" 2>/dev/null)
  thi=$(awk -F, 'NR==2{print $4}' "$ANALYSIS/records/$label.csv" 2>/dev/null)
  if [[ -z "$tlo" || ! -f "$exp/times.txt" ]]; then echo "0 0"; return; fi
  awk -v lo="$tlo" -v hi="$thi" '{t=$1+0; if (t < lo-1e-9) s++; else if (t <= hi+1e-9) m++} END{print s+0, m+0}' "$exp/times.txt"
}

field () { awk -F'= *' -v k="$2" '$1 ~ "^"k {print $2; exit}' "$1" 2>/dev/null; }

# --- band study ------------------------------------------------------------
do_band () {
  local ls=("$@"); [[ ${#ls[@]} -eq 0 ]] && ls=(sw050_16 sw070_16 meas32)
  local tab="$STUDIES/band_study.csv"
  echo "label,kr0,band_cells,median_rate,iqr,ranks_with_rate,visc_exact" > "$tab"
  echo "== band study: ${ls[*]}"
  for label in "${ls[@]}"; do
    local exp="$ANALYSIS/$label/export_ls" out="$ANALYSIS/$label/studies" k sk mx b
    [[ -f "$exp/perturbation.bin" ]] || { echo "  $label: no level-set export, skipped"; continue; }
    mkdir -p "$out"; k=$(kr0_of "$label"); read -r sk mx <<< "$(window_of "$label" "$exp")"
    for b in 1 2 3 4; do
      if need "$out/band$b/rank_scan.txt"; then
        ./bin/run_lsa analyse --from-export "$exp" --kr0 "$k" --weight volume --band "$b" \
          --skip "$sk" --max-snapshots "$mx" --rank-scan 2:12 --archive "$out/band$b" \
          > "$out/band$b.log" 2>&1 || true
      fi
      echo "$label,$k,$b,$(field "$out/band$b/rank_scan.txt" median_rate),$(field "$out/band$b/rank_scan.txt" iqr),$(field "$out/band$b/rank_scan.txt" ranks_with_rate),$(field "$out/band$b/rank_scan.txt" visc_exact)" >> "$tab"
    done
    echo "  $label: done (window skip $sk, $mx snapshots)"
  done
  awk -F, 'NR==1{printf "%-10s %5s %5s %12s %10s %10s\n","label","kr0","band","median","iqr","vs_visc_%"; next}
           {printf "%-10s %5s %5s %12.6f %10.2e %+10.3f\n",$1,$2,$3,$4,$5,($7>0?100*($4/$7-1):0)}' "$tab"
}

# --- stride study ------------------------------------------------------------
do_stride () {
  local ls=("$@"); [[ ${#ls[@]} -eq 0 ]] && ls=(sw020_16 sw030_16 sw040_16 sw050_16 sw060_16 sw070_16 sw080_16 sw090_16)
  echo "== stride study: ${ls[*]}"
  for label in "${ls[@]}"; do
    local out="$ANALYSIS/$label/studies" k sk mx hi starts
    mkdir -p "$out"; k=$(kr0_of "$label")
    for kind in vel ls; do
      local exp opts
      if [[ $kind == vel ]]; then exp="$ANALYSIS/$label/export"; opts=(--weight volume)
      else exp="$ANALYSIS/$label/export_ls"; opts=(--weight volume --band 2); fi
      [[ -f "$exp/perturbation.bin" ]] || { echo "  $label/$kind: no export, skipped"; continue; }
      read -r sk mx <<< "$(window_of "$label" "$exp")"
      hi=$(( sk + mx )); (( mx == 0 )) && hi=0
      starts="$sk"
      if need "$out/stride_$kind.csv"; then
        ./bin/run_lsa converge --from-export "$exp" --kr0 "$k" "${opts[@]}" --max-snapshots "$hi" \
          --starts "$starts" --strides 1,2,4,8,16 --ranks 4,8,16,0 --csv "$out/stride_$kind.csv" \
          > "$out/stride_$kind.log" 2>&1 || true
      fi
    done
    echo "  $label: done"
  done
  # summary: median over ranks of the rate at each stride, using all snapshots of the window
  local tab="$STUDIES/stride_study.csv"
  echo "label,kind,stride,dt,median_rate,n_ranks,visc_exact_dep_pct_median" > "$tab"
  for label in "${ls[@]}"; do
    for kind in vel ls; do
      local f="$ANALYSIS/$label/studies/stride_$kind.csv"
      [[ -f "$f" ]] || continue
      awk -F, -v l="$label" -v kd="$kind" 'NR>1 && $7!="NaN" && $7!="" {key=$1; g[key]=g[key]" "$7; e[key]=e[key]" "$8; dt[key]=$2}
        END{for (s in g){n=split(g[s],a," "); m=split(e[s],b," "); asort(a); asort(b);
            med=(n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2); meb=(m%2?b[(m+1)/2]:(b[m/2]+b[m/2+1])/2);
            printf "%s,%s,%s,%s,%.8f,%d,%.4f\n", l, kd, s, dt[s], med, n, meb}}' "$f" | sort -t, -k3,3n >> "$tab"
    done
  done
  column -t -s, "$tab" 2>/dev/null || cat "$tab"
}

# --- base-state drift ------------------------------------------------------------
do_drift () {
  local ls=("$@"); [[ ${#ls[@]} -eq 0 ]] && ls=(sw070_16)
  echo "== base-state drift: ${ls[*]}"
  for label in "${ls[@]}"; do
    local base="$RUNROOT/$label.base" out="$STUDIES/drift_$label"
    [[ -d "$base" ]] || base="$RUNROOT/$label"
    if need "$out.csv"; then
      ./bin/run_lsa drift --base "$base" --fields L0101,x_velocity,y_velocity --csv "$out.csv" > "$out.txt" 2>&1 || true
    fi
    cat "$out.txt"
  done
}

# --- DMD archives with modes ------------------------------------------------------------
do_modes () {
  local ls=("$@"); [[ ${#ls[@]} -eq 0 ]] && ls=(sw040_16 sw050_16 sw060_16 sw070_16 meas32)
  echo "== DMD archives with modes: ${ls[*]}"
  for label in "${ls[@]}"; do
    local exp="$ANALYSIS/$label/export" out="$ANALYSIS/$label/studies/modes" k sk mx
    [[ -f "$exp/perturbation.bin" ]] || { echo "  $label: no velocity export, skipped"; continue; }
    k=$(kr0_of "$label"); read -r sk mx <<< "$(window_of "$label" "$exp")"
    if need "$out/modes.bin"; then
      mkdir -p "$(dirname "$out")"
      ./bin/run_lsa analyse --from-export "$exp" --kr0 "$k" --weight volume --skip "$sk" \
        --max-snapshots "$mx" --rank 8 --n-keep 4 --archive "$out" > "$out.log" 2>&1 || true
    fi
    echo "  $label: $(field "$out/summary.txt" leading_stationary_growth_rate)"
  done
}

case "$what" in
  band)   do_band "${labels[@]}" ;;
  stride) do_stride "${labels[@]}" ;;
  drift)  do_drift "${labels[@]}" ;;
  modes)  do_modes "${labels[@]}" ;;
  all)    do_band; do_stride; do_drift; do_modes ;;
  *) echo "usage: $0 band|stride|drift|modes|all [labels...]" >&2; exit 1 ;;
esac
