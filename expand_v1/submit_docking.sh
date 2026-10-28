#!/bin/bash
# 확대 cohort 도킹 제출. 입력 준비 array가 끝난 뒤 실행한다.
set -euo pipefail
D=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1
/gpfs/jaeohshin/work/envs/dock-vina/bin/python "$D/code/make_tasks_expand.py"
N=$(wc -l < "$D/docked_tasks.txt")
echo "tasks: $N"
sbatch --array=0-$((N-1))%60 "$D/code/dock.sbatch"
