#!/bin/bash
#SBATCH -J dock_vina_pqr
#SBATCH -p a10
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:40:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/docking_v1/logs/%x_%j.log
# Add PDB2PQR + PROPKA (receptor protonation at pH 7.4) to dock-vina (user-approved 2026-10-28).
set -euo pipefail
ENV=/gpfs/jaeohshin/work/envs/dock-vina
CONDA=/home/jaeohshin/miniforge3/bin/conda
$CONDA install -y -p "$ENV" -c conda-forge --override-channels pdb2pqr propka || "$ENV/bin/pip" install --no-cache-dir pdb2pqr propka
$CONDA clean -a -y
"$ENV/bin/pdb2pqr30" --version || "$ENV/bin/pdb2pqr" --version
"$ENV/bin/python" -c "import propka, vina, meeko; print('propka', propka.__version__, 'vina', vina.__version__, 'meeko', meeko.__version__)"
$CONDA list -p "$ENV" --export > /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/docking_v1/dock-vina_env.txt
