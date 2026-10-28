#!/bin/bash
#SBATCH -J dock_vina_env
#SBATCH -p a10
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/docking_v1/logs/%x_%j.log
# Docking environment for the 2x2 pilot (user-approved 2026-10-28): AutoDock Vina + Meeko + RDKit + openbabel, conda-forge.
set -euo pipefail
ENV=/gpfs/jaeohshin/work/envs/dock-vina
CONDA=/home/jaeohshin/miniforge3/bin/conda
curl -sS -o /dev/null -w "network %{http_code}\n" https://conda.anaconda.org/conda-forge/channeldata.json
[ -d "$ENV" ] && { echo "exists: $ENV"; exit 1; }
$CONDA create -y -p "$ENV" -c conda-forge --override-channels python=3.11 vina meeko rdkit openbabel gemmi numpy scipy pandas pdb-tools
$CONDA clean -a -y
"$ENV/bin/python" - <<'P'
import vina, meeko, rdkit, gemmi, numpy, scipy, pandas
from openbabel import openbabel
print('vina', vina.__version__, 'meeko', meeko.__version__, 'rdkit', rdkit.__version__, 'gemmi', gemmi.__version__, 'openbabel', openbabel.OBReleaseVersion())
P
ls "$ENV/bin" | grep -E "^(vina|obabel|mk_prepare_ligand.py|mk_prepare_receptor.py)$"
$CONDA list -p "$ENV" --export > /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/docking_v1/dock-vina_env.txt
