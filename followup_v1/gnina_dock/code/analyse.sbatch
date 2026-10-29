#!/bin/bash
#SBATCH -J gnina_analyse
#SBATCH -p a10
#SBATCH -c 4
#SBATCH --mem=32G

#SBATCH -t 03:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/gnina_dock/logs/%x_%j.log
set -uo pipefail
cd /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/followup_v1/gnina_dock
/gpfs/jaeohshin/work/envs/dock-vina/bin/python -u code/analyse_gnina.py
