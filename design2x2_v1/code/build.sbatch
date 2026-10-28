#!/bin/bash
#SBATCH -J dp_build_receptors
#SBATCH -p a10
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/design2x2_v1/logs/%x_%j.log
# CPU-only receptor build from qFit models (no docking, inference or training).
set -euo pipefail
cd /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/design2x2_v1
/gpfs/jaeohshin/work/envs/multistate-analysis/bin/python -u code/build_receptors.py
