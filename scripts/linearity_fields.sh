#!/bin/bash
# linearity_fields.sh
#
# Field-direction linearity test (thesis Section 5.5, review item 4) for the three amplitude sets at k r0 = 0.7:
#   16 cells: sw070_16 (eps 2.5e-3), lin050_16 (5e-3), lin100_16 (1e-2), twin sw070_16.base
#    8 cells: amp08_025, amp08_050, amp08_100, twin amp08.base
#   32 cells: amp32_025, amp32_050, amp32_100, twin amp32.base
#
# 1. Writes the lossless twin-difference exports (velocity and level set) that do not exist yet. Only sw070_16 was exported by the campaign analysis; the others are exported here from their plotfiles, against the shared twin of their set.
# 2. Runs plotting/linearity_fields.py on velocity and on level set, up to the time at which the largest amplitude of the set reaches a/r0 = 0.1 (column 4 of its record csv, as for the norm test).
#
# Run from $LSA (the LSA_Rayleigh_Fortran folder inside the solver tree), after git pull and install_into_solver_tree.sh:
#   sbatch -A genacc_q -n 1 -c 4 --mem=32G -t 04:00:00 -o runs/array/linfields-%j.out scripts/linearity_fields.sh
# Results: runs/analysis/linearity_fields/*.txt and *.csv (small; collect_results.sh picks them up, or copy the folder).
set -u
R=runs
A=runs/analysis
OUT=$A/linearity_fields
mkdir -p "$OUT"

export_one() {   # label base
  local c=$1 b=$2
  [[ -d "$R/$c" ]] || { echo "  $c: run folder missing (plotfiles deleted?), skipped"; return; }
  mkdir -p "$A/$c"
  [[ -s "$A/$c/export/perturbation.bin" ]] || ./bin/run_lsa export --perturbed "$R/$c" --base "$R/$b" \
      --fields x_velocity,y_velocity --export "$A/$c/export" > "$A/$c/export.log" 2>&1 || echo "  $c: velocity export FAILED"
  [[ -s "$A/$c/export_ls/perturbation.bin" ]] || ./bin/run_lsa export --perturbed "$R/$c" --base "$R/$b" \
      --fields L0101 --export "$A/$c/export_ls" > "$A/$c/export_ls.log" 2>&1 || echo "  $c: level-set export FAILED"
}

run_set() {   # name largest-label eps1=label1 eps2=label2 eps3=label3
  local name=$1 big=$2; shift 2
  local tmax; tmax=$(awk -F, 'NR==2{print $4}' "$A/records/$big.csv" 2>/dev/null)
  [[ -n "$tmax" ]] || { echo "  $name: no record csv for $big, skipped"; return; }
  for kind in export export_ls; do
    local args=()
    for item in "$@"; do args+=("${item%%=*}=$A/${item#*=}/$kind"); done
    python3 plotting/linearity_fields.py --tmax "$tmax" --csv "$OUT/${name}_${kind}.csv" "${args[@]}" > "$OUT/${name}_${kind}.txt" 2>&1
    echo "== $name $kind (t <= $tmax)"; grep "max rel_diff" "$OUT/${name}_${kind}.txt"
    if [[ $kind == export_ls ]]; then
      # the level set restricted to the narrow band |r - r0| < 2 cells, the observable the decomposition uses
      python3 plotting/linearity_fields.py --band 2 --tmax "$tmax" --csv "$OUT/${name}_${kind}_band2.csv" "${args[@]}" > "$OUT/${name}_${kind}_band2.txt" 2>&1
      echo "== $name ${kind}_band2 (t <= $tmax)"; grep "max rel_diff" "$OUT/${name}_${kind}_band2.txt"
    fi
  done
}

for c in sw070_16 lin050_16 lin100_16; do export_one "$c" sw070_16.base; done
for c in amp08_025 amp08_050 amp08_100; do export_one "$c" amp08.base; done
for c in amp32_025 amp32_050 amp32_100; do export_one "$c" amp32.base; done

run_set lin16 lin100_16 2.5e-3=sw070_16 5e-3=lin050_16 1e-2=lin100_16
run_set lin08 amp08_100 2.5e-3=amp08_025 5e-3=amp08_050 1e-2=amp08_100
run_set lin32 amp32_100 2.5e-3=amp32_025 5e-3=amp32_050 1e-2=amp32_100
