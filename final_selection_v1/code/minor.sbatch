#!/bin/bash
#SBATCH -J dp_minor_state
#SBATCH -p a10
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/final_selection_v1/logs/%x_%j.log
set -euo pipefail
cd /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/final_selection_v1
/gpfs/jaeohshin/work/envs/multistate-analysis/bin/python -u code/minor_state_residues.py
