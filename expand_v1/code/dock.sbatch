#!/bin/bash
#SBATCH -J dock_expand
#SBATCH -p normal
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH -t 08:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1/logs/%x_%A_%a.log
# 2x2 docking pilot: one array task per (case, arm). Search budget is matched across arms (ANALYSIS_PLAN.md).
set -euo pipefail
D=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1
P=/gpfs/jaeohshin/work/envs/dock-vina/bin/python
mapfile -t JOBS < "$D/docked_tasks.txt"
IFS=' ' read -r CASE ARM <<< "${JOBS[$SLURM_ARRAY_TASK_ID]}"
$P -u "$D/code/run_docking.py" --case "$CASE" --arm "$ARM" --cpu 8
