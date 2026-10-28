#!/bin/bash
#SBATCH -J dp_qfit_states
#SBATCH -p a10
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/state_audit_v1/logs/%x_%j.log
# CPU-only qFit pocket-state audit (no docking, inference or training).
set -euo pipefail
cd /gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/state_audit_v1
/gpfs/jaeohshin/work/envs/multistate-analysis/bin/python -u code/audit_qfit_states.py
