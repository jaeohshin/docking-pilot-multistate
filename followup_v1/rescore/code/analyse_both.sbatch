#!/bin/bash
#SBATCH -J rescore_both
#SBATCH -p normal
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 03:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/rescore/logs/%x_%j.log
# CNNscore 를 주 선택자로 쓴 재채점 재집계. CNNaffinity 결과(out/expand8)는 그대로 둔다.
set -uo pipefail
D=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/rescore
/gpfs/jaeohshin/work/envs/dock-vina/bin/python -u "$D/code/analyse_rescore.py" \
  --run /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1 --pose-table "$D/out/expand8/pose_table.tsv.gz" \
  --scores vinardo="$D/out/expand8/vinardo" gnina="$D/out/expand8g/gnina" \
  --out "$D/out/expand8_cnnscore"
