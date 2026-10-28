#!/bin/bash
#SBATCH -J dp_assembly_audit
#SBATCH -p a10
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/final_selection_v1/logs/%x_%j.log
# CPU-only structure audit of the 79 fixed candidates (no docking, inference or training).
set -euo pipefail
cd /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/final_selection_v1
/gpfs/jaeohshin/work/envs/multistate-analysis/bin/python -u code/audit_assembly_contacts.py
