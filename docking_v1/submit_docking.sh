#!/bin/bash
# 2x2 도킹 제출 (grammar, CPU). 입력 준비가 끝난 뒤 실행한다.
#   ssh grammar 'bash /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/docking_v1/submit_docking.sh'
set -euo pipefail
D=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/docking_v1
/gpfs/jaeohshin/work/envs/dock-vina/bin/python "$D/code/make_tasks.py"
N=$(wc -l < "$D/docked_tasks.txt")
echo "tasks: $N"
sbatch --array=0-$((N-1))%40 "$D/code/dock.sbatch"
