#!/bin/bash
#SBATCH -J dp_expand_audit
#SBATCH -p normal
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1/logs/%x_%j.log
# C 확대: packer용 상동성 차단을 푼 1,209 후보에 같은 기하 QC를 적용한다 (CPU, 도킹 없음).
set -euo pipefail
cd /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1
/gpfs/jaeohshin/work/envs/multistate-analysis/bin/python -u code/audit_geometry_expanded.py
