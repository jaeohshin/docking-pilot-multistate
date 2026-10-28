#!/bin/bash
#SBATCH -J chain_b50
#SBATCH -p normal
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH -t 00:20:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/rescore/logs/%x_%j.log
# budget50 도킹이 끝난 뒤에만 task 목록을 만들 수 있다(docking_log.json이 그때 생긴다).
# 그래서 목록 생성과 재채점 array 제출을 이 작업이 대신한다.
set -euo pipefail
D=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/rescore
RUN=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/budget50
P=/gpfs/jaeohshin/work/envs/dock-vina/bin/python
$P "$D/code/make_tasks.py" --run "$RUN" --tag b50
N=$(wc -l < "$D/tasks_b50.txt")
echo "b50 rescore tasks: $N"
PT=$(sbatch --parsable --export=ALL,RUN="$RUN",TAG=b50 "$D/code/posetable.sbatch")
DEPS="$PT"
OFF=0
while [ $OFF -lt $N ]; do
  LAST=$(( N - OFF - 1 )); [ $LAST -gt 999 ] && LAST=999
  J=$(sbatch --parsable --export=ALL,RUN="$RUN",TAG=b50,OFFSET=$OFF --array=0-$LAST%80 "$D/code/rescore.sbatch")
  echo "shard offset $OFF size $((LAST+1)) -> $J"
  DEPS="$DEPS:$J"
  OFF=$(( OFF + LAST + 1 ))
done
sbatch --dependency=afterany:$DEPS --export=ALL,RUN="$RUN",TAG=b50 "$D/code/analyse.sbatch"
