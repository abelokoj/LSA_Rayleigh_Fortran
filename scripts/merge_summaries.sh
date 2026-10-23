#!/bin/bash
# =====================================================================
# Merge the per-task summaries written by a Slurm array run.
# =====================================================================
# Each array task writes its own summary so that concurrent tasks do not interleave lines in a shared file. This combines them into runs/summary.tsv, in the same format the sequential runner produces, so that the plotting and reporting code need not know which route produced the results.
#
# Run from amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Fortran:
#
#     ./scripts/merge_summaries.sh            # every array job, newest result per case
#     ./scripts/merge_summaries.sh 1234567    # one array job only
#
# A case resubmitted after a time-out or a pre-emption appears in the summaries of more than one array job. Only its newest row is kept: summaries are read in order of job ID, and a later row for a label replaces an earlier one. Deduplicating whole lines was not enough, because the stale no_record or incomplete row and the new ok row differ and both survived, so the plots read the case twice.
# ---------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

SRC="runs/array"
OUT="runs/summary.tsv"

shopt -s nullglob
JOBID="${1:-}"
if [[ -n "$JOBID" ]]; then
  parts=("$SRC"/summary_"${JOBID}"_*.tsv)
else
  parts=("$SRC"/summary_*.tsv)
fi
if (( ${#parts[@]} == 0 )); then
  echo "no per-task summaries under $SRC; nothing to merge" >&2
  exit 1
fi

# Header from the first part, then every data line from all of them, sorted by label so that the merged file is reproducible regardless of the order the tasks happened to finish in.
# Oldest job first (sort -V orders the numeric job and task IDs), so that the last row seen for a label is the newest.
mapfile -t parts < <(printf '%s\n' "${parts[@]}" | sort -V)
head -n 1 "${parts[0]}" > "$OUT"
for f in "${parts[@]}"; do tail -n +2 "$f"; done \
  | awk -F'\t' 'NF { row[$1] = $0 } END { for (k in row) print row[k] }' \
  | sort -t $'\t' -k1,1 >> "$OUT"

n=$(( $(wc -l < "$OUT") - 1 ))
echo "merged ${#parts[@]} task summaries into $OUT ($n case(s))"

if command -v column >/dev/null 2>&1; then
  column -t -s $'\t' "$OUT"
fi
