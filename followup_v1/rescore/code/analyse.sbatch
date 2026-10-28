#!/bin/bash
#SBATCH -J rescore_analyse
#SBATCH -p normal
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 03:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/rescore/logs/%x_%j.log
set -euo pipefail
D=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/rescore
/gpfs/jaeohshin/work/envs/dock-vina/bin/python -u "$D/code/analyse_rescore.py" \
  --run "$RUN" --pose-table "$D/out/$TAG/pose_table.tsv.gz" \
  --scores vinardo="$D/out/$TAG/vinardo" gnina="$D/out/$TAG/gnina" --out "$D/out/$TAG"
