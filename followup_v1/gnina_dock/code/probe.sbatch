#!/bin/bash
#SBATCH -J gnina_probe
#SBATCH -p a10
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH -t 01:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/gnina_dock/logs/%x_%j.log
# gnina 탐색 비용만 측정한다. pose는 probe/ 에 쓰고 본 실행에서 쓰지 않는다(PLAN.md).
set -uo pipefail          # -e 는 쓰지 않는다: ls|head 의 SIGPIPE 로 조용히 죽는다
G=/gpfs/jaeohshin/work/envs/gnina-bin/gnina
B=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028
D=$B/expand_v1
O=$B/followup_v1/gnina_dock/probe
P=/gpfs/jaeohshin/work/envs/dock-vina/bin/python
nvidia-smi --query-gpu=name --format=csv,noheader
$G --version | head -2
mapfile -t CASES < <(ls "$D/docked")
for CASE in "${CASES[@]:0:3}"; do
  read -r CX CY CZ SX SY SZ REC <<< "$($P -c "
import json
m={c['case_id']:c for c in json.load(open('$D/inputs_manifest.json'))['cases']}['$CASE']['apo']
print(*m['box_center'], *m['box_size'], m['receptors']['single_0'])")"
  echo "case $CASE rec $REC box $SX"
  for EX in 8 32; do
    S=$(date +%s)
    $G --receptor "$D/inputs/$CASE/$REC" --ligand "$D/inputs/$CASE/ligand_dock.pdbqt" \
       --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
       --size_x "$SX" --size_y "$SY" --size_z "$SZ" \
       --exhaustiveness "$EX" --num_modes 9 --seed 20261028 --cpu 8 \
       --out "$O/${CASE}_ex${EX}.sdf" > "$O/${CASE}_ex${EX}.log" 2>&1
    echo "TIME case=$CASE exh=$EX rc=$? seconds=$(( $(date +%s) - S ))"
    tail -2 "$O/${CASE}_ex${EX}.log" | sed 's/^/    /'
  done
done
echo PROBE_DONE
