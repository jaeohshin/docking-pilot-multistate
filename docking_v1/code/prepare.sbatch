#!/bin/bash
#SBATCH -J dock_prepare
#SBATCH -p a10
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/docking_v1/logs/%x_%j.log
# Ligand/receptor preparation for the 2x2 pilot (no docking).
set -euo pipefail
cd /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/docking_v1
/gpfs/jaeohshin/work/envs/dock-vina/bin/python -u code/prepare_inputs.py "$@"
