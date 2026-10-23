#!/usr/bin/env bash
# =====================================================================
# refit_records.sh -- recompute every interface growth rate of a results snapshot from its amplitude records, with the Fortran tool that produced them, and compare with the stored values.
# =====================================================================
# Run from the top of this repository after `make`:
#
#     scripts/refit_records.sh                                   # results/followup_2026-10-22 (every case of the dissertation)
#     scripts/refit_records.sh results/campaign_2026-10-09       # another snapshot
#     scripts/refit_records.sh runs                              # a new set of runs (records in runs/records, summary in runs/analysis)
#
# For every case in <snapshot>/analysis/summary.csv that has a record <snapshot>/records/<label>.dat, it runs
#
#     bin/run_lsa record --record <snapshot>/records/<label>.dat --kr0 <k> [gas options of the case] --csv <out>/<label>.csv
#
# exactly as the analysis job did (scripts/analyse_campaign.sh, step 1), and prints the stored and recomputed rate side by side. The gas options (--rho-g, --mu-g) are read from the case line of the run's run record (results/*/run_records/<case>/manifest.txt); cases without a manifest use the campaign fluids, which are run_lsa's defaults. The full report of each fit (window, e-folds, start-time spread, departure from each theory) is written to <out>/<label>.txt. Nothing in results/ is modified.
# ---------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
SNAP="${1:-results/followup_2026-10-22}"
OUT="${OUT:-runs/refit/$(basename "$SNAP")}"
[[ -x bin/run_lsa ]] || { echo "bin/run_lsa not found: run make first" >&2; exit 1; }
SUM="$SNAP/analysis/summary.csv"
[[ -f "$SUM" ]] || { echo "no $SUM" >&2; exit 1; }
mkdir -p "$OUT"

# gas options of a case, from the case line recorded in its manifest
opts_of () {
  local m
  m=$(ls results/*/run_records/"$1"/manifest.txt "$SNAP"/run_records/"$1"/manifest.txt runs/"$1"/manifest.txt 2>/dev/null | tail -1 || true)
  [[ -n "$m" ]] || return 0
  grep -m1 '^case_line' "$m" | tr ' ' '\n' | awk '/^--(rho-g|mu-g|mu-l)$/{o=$0; next} o!=""{printf "%s %s ", o, $0; o=""}'
}

printf '%-16s %6s %12s %12s %10s\n' case kr0 stored recomputed 'diff %'
n=0; worst=0
while IFS=, read -r label kr0 _ _ _ rate _; do
  [[ "$label" == label ]] && continue
  label=$(echo "$label" | tr -d ' '); kr0=$(echo "$kr0" | tr -d ' '); rate=$(echo "$rate" | tr -d ' ')
  rec="$SNAP/records/$label.dat"
  [[ -f "$rec" && -n "$rate" ]] || continue
  # shellcheck disable=SC2046
  ./bin/run_lsa record --record "$rec" --kr0 "$kr0" $(opts_of "$label") --csv "$OUT/$label.csv" > "$OUT/$label.txt" 2>&1 || { echo "$label: run_lsa record failed, see $OUT/$label.txt"; continue; }
  new=$(awk -F, 'NR==1{for(i=1;i<=NF;i++) if($i=="rate") c=i; next} NR==2{gsub(/ /,"",$c); print $c}' "$OUT/$label.csv")
  d=$(awk -v a="$rate" -v b="$new" 'BEGIN{printf "%+.4f", 100*(b/a-1)}')
  printf '%-16s %6s %12.6f %12.6f %10s\n' "$label" "$kr0" "$rate" "$new" "$d"
  worst=$(awk -v w="$worst" -v d="$d" 'BEGIN{d=(d<0)?-d:d; print (d>w)?d:w}')
  n=$((n + 1))
done < "$SUM"
echo
echo "$n records refitted; largest difference from the stored rate: $worst %  (reports in $OUT/)"
