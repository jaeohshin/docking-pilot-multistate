#!/bin/bash
#SBATCH -J dp_expand_prep
#SBATCH -p normal
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 06:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1/logs/%x_%A_%a.log
set -euo pipefail
cd /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1
/gpfs/jaeohshin/work/envs/dock-vina/bin/python -u code/build_and_prep.py --shard $SLURM_ARRAY_TASK_ID --shards 20
