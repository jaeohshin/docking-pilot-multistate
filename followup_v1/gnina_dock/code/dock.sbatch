#!/bin/bash
#SBATCH -J gnina_dock
#SBATCH -p a10,a40
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --gres=gpu:1
#SBATCH -t 04:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/gnina_dock/logs/%x_%A_%a.log
# 실험 C: gnina가 탐색과 선택을 모두 한다. arm 간 예산은 동일(PLAN.md).
set -uo pipefail
D=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/gnina_dock
P=/gpfs/jaeohshin/work/envs/dock-vina/bin/python
mapfile -t JOBS < "$D/docked_tasks.txt"
IFS=' ' read -r CASE ARM <<< "${JOBS[$((SLURM_ARRAY_TASK_ID + ${OFFSET:-0}))]}"
$P -u "$D/code/run_gnina_dock.py" --case "$CASE" --arm "$ARM" --cpu 8
