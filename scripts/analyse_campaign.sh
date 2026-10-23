#!/bin/bash
# =====================================================================
# analyse_campaign.sh -- every post-processing step, after the solver runs.
# =====================================================================
# Runs on the files the sweep left in runs/ and writes everything to runs/analysis/. Each step is skipped when its output already exists, so the script may be re-run after more cases finish and only the new work is done. It needs no Python; the figures are drawn afterwards from runs/analysis on a workstation.
#
# What it does, following the recommendations of the September 2026 evidence report:
#
#   1. records   run_lsa record on every interface-amplitude record: the admissible-window growth rate (after the startup transient, before a/r0 = 0.1), with the window and a start-time sensitivity table printed beside it -> runs/analysis/records/<label>.txt and .csv
#   2. export    for every perturbed case with a .base twin, the matched perturbation sequence Q_p - Q_b written once as a lossless binary (run_lsa export) -> runs/analysis/<label>/export. Every later DMD runs from this file, not from the plotfiles, which are then no longer needed for the analysis.
#   3. dmd       run_lsa analyse on the export over the SAME window the interface fit used, with the singular values, eigenvalues, residual and leading modes archived -> runs/analysis/<label>/dmd
#   4. converge  the DMD sensitivity study: start offset x snapshot count x sampling stride x rank -> runs/analysis/<label>/converge.csv
#   4b. linearity  full-state proportionality of the amplitude set at k r0 = 0.7 (run_lsa linearity) -> runs/analysis/linearity_16.txt
#   5. summary   one line per case comparing the interface rate, the DMD rate and theory -> runs/analysis/summary.csv
#
# Submit as a small job after the arrays (scripts/submit_campaign.sh does this with a dependency), or run directly in an interactive allocation:
#
#     sbatch -A genacc_q -n 1 -c 4 -t 04:00:00 --mem=16G scripts/analyse_campaign.sh
#
# Environment: FIELDS (default x_velocity,y_velocity), WEIGHT (default volume), RUNROOT (default runs), FORCE=1 to redo steps whose output exists.
# ---------------------------------------------------------------------
set -u
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

RUNROOT="${RUNROOT:-runs}"
OUT="$RUNROOT/analysis"
FIELDS="${FIELDS:-x_velocity,y_velocity}"
WEIGHT="${WEIGHT:-volume}"
FORCE="${FORCE:-0}"
mkdir -p "$OUT/records"

[[ -x bin/run_lsa ]] || { echo "bin/run_lsa missing; run make first" >&2; exit 1; }

need () { [[ "$FORCE" == 1 || ! -e "$1" ]]; }
# FORCE_DMD=1 redoes the DMD and sensitivity steps (from the existing exports) without re-reading the plotfiles.
need_dmd () { [[ "$FORCE" == 1 || "${FORCE_DMD:-0}" == 1 || ! -e "$1" ]]; }
# FORCE_RECORDS=1 redoes only the record fits (cheap; no plotfiles read), e.g. after a new theory column was added.
need_rec () { [[ "$FORCE" == 1 || "${FORCE_RECORDS:-0}" == 1 || ! -e "$1" ]]; }

# Fluid properties of a case from its deck (ns.denconst = rho_l rho_g, ns.viscconst = mu_l mu_g), as run_lsa options, so that the two-fluid theory uses the gas of that case.
props_of () {
  local deck="$RUNROOT/$1/inputs"
  [[ -f "$deck" ]] || return 0
  awk '/^ns\.denconst[ \t]*=/{sub(/#.*/,""); split($0,a,"="); split(a[2],d," "); rg=d[2]}
       /^ns\.viscconst[ \t]*=/{sub(/#.*/,""); split($0,a,"="); split(a[2],v," "); ml=v[1]; mg=v[2]}
       END{if (ml!="") printf "--mu-l %s ", ml; if (rg!="") printf "--rho-g %s ", rg; if (mg!="") printf "--mu-g %s", mg}' "$deck"
}

# k r0 of a case: from its deck header, else from the summary.
kr0_of () {
  local label="$1" k
  k=$(awk '/^#[ \t]*k r0[ \t]*=/{print $NF; exit}' "$RUNROOT/$label/inputs" 2>/dev/null)
  [[ -z "$k" ]] && k=$(awk -F'\t' -v l="$label" '$1==l{k=$2} END{print k}' "$RUNROOT"/summary.tsv "$RUNROOT"/array/summary_*.tsv 2>/dev/null)
  echo "${k:-0.7}"
}

# --- 1. interface records ---------------------------------------------
echo "== 1. interface-amplitude records"
for rec in "$RUNROOT"/records/*.dat; do  # harmonic records of broadband runs live in records/<label>_harmonics/ and are analysed by scripts/analyse_broadband.sh
  [[ -e "$rec" ]] || continue
  label=$(basename "$rec" .dat)
  case "$label" in *.base) continue ;; esac
  txt="$OUT/records/$label.txt"
  if need_rec "$txt"; then
    k=$(kr0_of "$label")
    # shellcheck disable=SC2046
    ./bin/run_lsa record --record "$rec" --kr0 "$k" $(props_of "$label") --csv "$OUT/records/$label.csv" \
        --starts-t "$(awk 'NR>1 && $1+0>0 {t[n++]=$1} END{for(i=1;i<=5;i++) printf "%s%s", (i>1?",":""), t[int(n*(0.45+0.08*i))]}' "$rec")" \
        > "$txt" 2>&1 || true
    echo "  $label: $(grep -E 'growth rate s' "$txt" | head -1)"
  fi
done

# --- 2-4. twin pairs ---------------------------------------------------
echo "== 2-4. twin pairs: export, DMD, sensitivity"
for bdir in "$RUNROOT"/*.base; do
  [[ -d "$bdir" ]] || continue
  label=$(basename "$bdir" .base)
  pdir="$RUNROOT/$label"
  [[ -d "$pdir" ]] || continue
  # Only pairs whose two runs are both finished: run_all.sh writes runs/records/<label>.dat only when a case completes. Exporting a pair that is still running would freeze an incomplete sequence, which later passes skip because the export exists.
  if [[ ! -s "$RUNROOT/records/$label.dat" || ! -s "$RUNROOT/records/$label.base.dat" ]]; then
    echo "  $label: pair not finished yet, skipped"
    continue
  fi
  dest="$OUT/$label"; mkdir -p "$dest"
  k=$(kr0_of "$label")

  if need "$dest/export/perturbation.bin"; then
    echo "  $label: exporting the matched perturbation sequence"
    ./bin/run_lsa export --perturbed "$pdir" --base "$bdir" --fields "$FIELDS" \
        --export "$dest/export" > "$dest/export.log" 2>&1 \
      || { echo "  $label: export FAILED, see $dest/export.log"; continue; }
  fi

  # The DMD window is the interface-fit window: snapshots with t_lo <= t <= t_hi, counted on the export's own times.
  skip=0; maxs=0
  csv="$OUT/records/$label.csv"
  if [[ -s "$csv" ]]; then
    read -r tlo thi < <(awk -F, 'NR==2{print $3, $4}' "$csv")
    skip=$(awk -v a="$tlo" '$1+0 < a-1e-9 {n++} END{print n+0}' "$dest/export/times.txt")
    maxs=$(awk -v a="$tlo" -v b="$thi" '$1+0 >= a-1e-9 && $1+0 <= b+1e-9 {n++} END{print n+0}' "$dest/export/times.txt")
  fi
  if need_dmd "$dest/dmd/summary.txt"; then
    ./bin/run_lsa analyse --from-export "$dest/export" --kr0 "$k" --weight "$WEIGHT" \
        --skip "$skip" --max-snapshots "$maxs" --rank-scan 4:16 --archive "$dest/dmd" --save "$dest/dmd.json" \
        > "$dest/dmd.log" 2>&1 || true
    echo "  $label: DMD (snapshots $skip+1 .. $((skip+maxs))): $(grep -E 'DMD: leading' "$dest/dmd.log" | head -1)"
  fi
  # Level-set observable in a narrow band (added 2026-10-08). The velocity difference gave rates that moved by several per cent with the rank; the level set L0101 restricted to |r - r0| < 2 dr gave rank-insensitive rates (IQR 0.1-0.5 % over ranks 2-12) on meas08, meas16 and sw080_16. Exported once (raw, all rows), then decomposed in the band over the interface window.
  if need "$dest/export_ls/perturbation.bin"; then
    ./bin/run_lsa export --perturbed "$pdir" --base "$bdir" --fields L0101 \
        --export "$dest/export_ls" > "$dest/export_ls.log" 2>&1 \
      || echo "  $label: level-set export FAILED, see $dest/export_ls.log"
  fi
  if [[ -s "$dest/export_ls/perturbation.bin" ]] && need_dmd "$dest/dmd_ls/rank_scan.txt"; then
    ./bin/run_lsa analyse --from-export "$dest/export_ls" --kr0 "$k" --weight volume --band 2 \
        --skip "$skip" --max-snapshots "$maxs" --rank-scan 2:12 --archive "$dest/dmd_ls" \
        > "$dest/dmd_ls.log" 2>&1 || true
    echo "  $label: level-set DMD $(grep -E 'median rate' "$dest/dmd_ls.log" | head -1)"
  fi
  if need_dmd "$dest/converge.csv"; then
    # Same upper end as the fit; start offsets at 1/2, 3/4 and 1 of the admissible start and a quarter of the window after it, so that the table shows how the rate moves as the startup transient is included or excluded (the start-time sensitivity of the report), crossed with sampling strides and ranks.
    hi=$(( skip + maxs )); (( maxs == 0 )) && hi=0
    starts="$((skip/2)),$((3*skip/4)),$skip,$((skip + maxs/4))"
    ./bin/run_lsa converge --from-export "$dest/export" --kr0 "$k" --weight "$WEIGHT" \
        --max-snapshots "$hi" --starts "$starts" --strides 1,2,4 --ranks 4,8,16,0 \
        --csv "$dest/converge.csv" > "$dest/converge.log" 2>&1 || true
    echo "  $label: sensitivity table $(grep -c . "$dest/converge.csv" 2>/dev/null) lines"
  fi
done

# --- 4b. proportionality over the full record ------------------------
# Full-state |Q'| proportionality across the amplitude set of campaign_16.txt, if all its runs exist: eps = 2.5e-3 (sw070_16), 5e-3 (lin050_16), 1e-2 (lin100_16), one twin sw070_16.base. Runs are matched by time, so they must share the plot cadence, which they do by construction.
# All four runs must be finished (their records exist); on 2026-10-08 this step ran while three of them were still going and compared 9 common times of incomplete runs.
lin_ready=1
for l in sw070_16 sw070_16.base lin050_16 lin100_16; do [[ -s "$RUNROOT/records/$l.dat" ]] || lin_ready=0; done
if (( lin_ready )) && need "$OUT/linearity_16.txt"; then
  echo "== 4b. proportionality (k r0 = 0.7, 16 cells)"
  eps_sw=$(awk '$1=="ns.radblob"{print $3; exit}' "$RUNROOT/sw070_16/inputs")
  # Compare only up to the end of the linear window of the LARGEST amplitude (lin100_16 reaches a/r0 = 0.1 first). On 2026-10-09 the whole records were compared, out to t = 13 where lin100_16 has a/r0 ~ 0.4, and the late-time divergence gave a spurious NONLINEAR verdict.
  lin_tmax=$(awk -F, 'NR==2{print $4}' "$OUT/records/lin100_16.csv" 2>/dev/null)
  ./bin/run_lsa linearity --base "$RUNROOT/sw070_16.base" --kr0 0.7 --fields "$FIELDS" ${lin_tmax:+--t-max $lin_tmax} \
      --case "${eps_sw}=$RUNROOT/sw070_16" --case "5e-3=$RUNROOT/lin050_16" \
      --case "1e-2=$RUNROOT/lin100_16" > "$OUT/linearity_16.txt" 2>&1 || true
  grep -E "verdict|spread" "$OUT/linearity_16.txt" | sed 's/^/  /'
fi

# --- 5. summary ----------------------------------------------------------
echo "== 5. summary"
sum="$OUT/summary.csv"
echo "label,kr0,t_lo,t_hi,n,rate_interface,efolds,settled,rate_dmd,dmd_residual,dmd_rankscan_median,dmd_rankscan_iqr,dmd_levelset_median,dmd_levelset_iqr,visc_exact,inviscid,visc_two_fluid" > "$sum"
for csv in "$OUT"/records/*.csv; do
  [[ -e "$csv" ]] || continue
  label=$(basename "$csv" .csv)
  row=$(awk -F, 'NR==2{print $2","$3","$4","$5","$6","$8","$9}' "$csv")
  rd=$(awk -F'= *' '/^leading_stationary_growth_rate/{print $2}' "$OUT/$label/dmd/summary.txt" 2>/dev/null)
  rr=$(awk -F'= *' '/^reduced_map_residual/{print $2}' "$OUT/$label/dmd/summary.txt" 2>/dev/null)
  rm_=$(awk -F'= *' '/^median_rate/{print $2}' "$OUT/$label/dmd/rank_scan.txt" 2>/dev/null)
  ri_=$(awk -F'= *' '/^iqr/{print $2}' "$OUT/$label/dmd/rank_scan.txt" 2>/dev/null)
  th=$(awk -F, 'NR==2{print $11","$12","$13}' "$csv")
  lm_=$(awk -F'= *' '/^median_rate/{print $2}' "$OUT/$label/dmd_ls/rank_scan.txt" 2>/dev/null)
  li_=$(awk -F'= *' '/^iqr/{print $2}' "$OUT/$label/dmd_ls/rank_scan.txt" 2>/dev/null)
  echo "$label,$row,${rd:-},${rr:-},${rm_:-},${ri_:-},${lm_:-},${li_:-},$th" >> "$sum"
done
if command -v column >/dev/null 2>&1; then column -t -s, "$sum"; else cat "$sum"; fi
echo "written to $OUT"
