#!/bin/bash
# =====================================================================
# install_into_solver_tree.sh -- copy this code into the solver tree.
# =====================================================================
# Run on the cluster from the top of a clone of this repository. It copies the code (sources, tests, scripts, plotting, inputs, Makefile, case files, documents, results, archived records) into $AII/NS_fluids_lib/LSA_Rayleigh_Fortran, where the runs live and where the scripts expect the solver one level up. Nothing there is deleted: runs/, build/, bin/ and any file not in this folder are left as they are, and a dated copy of every file that is about to be replaced is kept first.
#
#     cd ~/LSA_Rayleigh_Fortran
#     scripts/install_into_solver_tree.sh                       # target $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
#     scripts/install_into_solver_tree.sh /some/other/target    # explicit target
#
# After a later git pull, run it again: only changed files are copied.
# ---------------------------------------------------------------------
set -euo pipefail
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${1:-${AII:?set AII=\$HOME/amrex_implicit_interfaces or give the target directory}/NS_fluids_lib/LSA_Rayleigh_Fortran}"

if [[ "$SRC" == "$(cd "$DEST" 2>/dev/null && pwd)" ]]; then
  echo "source and target are the same directory; nothing to do"; exit 0
fi
mkdir -p "$DEST"
stamp=$(date +%Y%m%d_%H%M%S)
backup="$DEST/.install_backup_$stamp"

# What is copied: everything in this folder except build products, run output and caches.
EXCLUDES=(--exclude=.git/ --exclude=runs/ --exclude=build/ --exclude=bin/ --exclude=figures/ --exclude=__pycache__/ --exclude='*.pyc' --exclude='.install_backup_*/')

echo "from: $SRC"
echo "to:   $DEST"
# rsync --backup keeps the previous version of every replaced file under $backup, with its relative path.
rsync -a --itemize-changes --backup --backup-dir="$backup" "${EXCLUDES[@]}" "$SRC"/ "$DEST"/ | grep -v '^\.d' || true
chmod +x "$DEST"/scripts/*.sh
if [[ -d "$backup" ]]; then
  echo "replaced files were saved in $backup"
else
  echo "no existing file was replaced"
fi
echo "next: cd $DEST && make -j8 && make test"
