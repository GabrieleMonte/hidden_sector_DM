#!/bin/bash
#SBATCH -J lam-crossing-survey-hP
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/lam_crossing_survey_hP_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/lam_crossing_survey_hP_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -q shared
#SBATCH -t 03:00:00
#
# Re-derive LAM_LIMS after the phase-entry fix, before rebuilding the grids.
#
#   sbatch GCE/scripts/lam_crossing_survey_hP.sbatch
#
# Shared QOS, 64 of a node's 256 cores: 30 nodes x 3 outputs = 90 sweeps, and
# the sweep stops one point past the crossing, so ~6 solves each at 4-9 s.
# Bills 64/256 of a node, about 0.2 node-hours against a ~215 node-hour
# balance.  PROBE_TIMEOUT caps any single solve at 300 s.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd /global/u1/g/gab97/secluded_DM

python GCE/scripts/lam_crossing_survey_hP.py 64
