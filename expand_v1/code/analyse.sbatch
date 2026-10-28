#!/bin/bash
#SBATCH -J dp_expand_analyse
#SBATCH -p normal
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1/logs/%x_%j.log
set -euo pipefail
cd /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1
/gpfs/jaeohshin/work/envs/dock-vina/bin/python -u code/analyse_expand.py
