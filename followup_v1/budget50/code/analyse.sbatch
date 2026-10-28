#!/bin/bash
#SBATCH -J dp_b50_analyse
#SBATCH -p normal
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/budget50/logs/%x_%j.log
set -euo pipefail
cd /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/budget50
/gpfs/jaeohshin/work/envs/dock-vina/bin/python -u code/analyse_budget50.py
