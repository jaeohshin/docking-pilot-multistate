#!/bin/bash
#SBATCH -J rescore
#SBATCH -p normal
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 06:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/rescore/logs/%x_%A_%a.log
set -euo pipefail
D=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/rescore
mapfile -t JOBS < "$D/tasks_$TAG.txt"
IFS=' ' read -r S CASE ARM <<< "${JOBS[$((SLURM_ARRAY_TASK_ID + ${OFFSET:-0}))]}"
/gpfs/jaeohshin/work/envs/dock-vina/bin/python -u "$D/code/rescore.py" \
  --run "$RUN" --case "$CASE" --arm "$ARM" --scoring "$S" --out "$D/out/$TAG/$S"
