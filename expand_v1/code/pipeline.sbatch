#!/bin/bash
#SBATCH -J dp_expand_pipe
#SBATCH -p normal
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 04:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1/logs/%x_%j.log
# 확대 cohort: assembly 감사 -> qFit 상태 감사 -> apo pose 재정합 (CPU, 도킹 없음)
set -euo pipefail
cd /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1
/gpfs/jaeohshin/work/envs/multistate-analysis/bin/python -u code/run_pipeline.py
