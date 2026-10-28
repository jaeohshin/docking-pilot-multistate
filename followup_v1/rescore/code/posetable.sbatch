#!/bin/bash
#SBATCH -J posetable
#SBATCH -p normal
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/rescore/logs/%x_%j.log
set -euo pipefail
D=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/rescore
/gpfs/jaeohshin/work/envs/dock-vina/bin/python -u "$D/code/pose_table.py" --run "$RUN" --out "$D/out/$TAG/pose_table.tsv.gz"
