#!/bin/bash
# =====================================================================
# Body of the post-processing job submitted by scripts/submit_campaign.sh; not intended to be run by hand.
# =====================================================================
# Merges the per-task summaries of the arrays, then runs every analysis step on what the arrays produced (scripts/analyse_campaign.sh). It is submitted with --dependency=afterany on the arrays, so it starts as soon as the last task ends, whether or not every task succeeded: the analysis skips what is missing, and a resubmission later only does the new work.
# ---------------------------------------------------------------------
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
set +u
module purge
module load gnu/13
module load openmpi
set -u
echo "post-processing job $SLURM_JOB_ID started $(date)"
./scripts/merge_summaries.sh || true
./scripts/analyse_campaign.sh
status=$?
echo "finished $(date), status $status"
exit $status
