#!/bin/bash
#SBATCH -J gnina_cpu_probe
#SBATCH -p normal
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/gnina_dock/logs/%x_%j.log
# grammar(CPU)에서 gnina 탐색 비용 측정. pose는 probe_cpu/ 에 쓰고 본 실행에 쓰지 않는다.
set -uo pipefail
G=/gpfs/jaeohshin/work/envs/gnina-bin/gnina
D=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/expand_v1
O=/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/gnina_dock/probe_cpu
P=/gpfs/jaeohshin/work/envs/dock-vina/bin/python
mkdir -p $O
$G --version|head -1
mapfile -t CASES < <(ls "$D/docked")
for CASE in "${CASES[@]:0:2}"; do
  read -r CX CY CZ SX SY SZ REC <<< "$($P -c "
import json
m={c['case_id']:c for c in json.load(open('$D/inputs_manifest.json'))['cases']}['$CASE']['apo']
print(*m['box_center'], *m['box_size'], m['receptors']['single_0'])")"
  $P -c "
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
pm=PDBQTMolecule.from_file('$D/inputs/$CASE/ligand_dock.pdbqt',skip_typing=True)
rd=RDKitMolCreate.from_pdbqt_mol(pm,only_cluster_leads=False)[0]
w=Chem.SDWriter('$O/$CASE.sdf');w.write(rd,confId=0);w.close()"
  S=$(date +%s)
  $G --receptor "$D/inputs/$CASE/$REC" --ligand "$O/$CASE.sdf" \
     --center_x $CX --center_y $CY --center_z $CZ --size_x $SX --size_y $SY --size_z $SZ \
     --exhaustiveness 32 --num_modes 9 --seed 20261028 --cpu 8 \
     --out "$O/${CASE}_out.sdf" > "$O/$CASE.log" 2>&1
  echo "TIME_CPU case=$CASE rc=$? seconds=$(( $(date +%s) - S ))"
done
echo CPU_PROBE_DONE
