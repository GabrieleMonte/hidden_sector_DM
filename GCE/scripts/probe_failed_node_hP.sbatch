#!/bin/bash
#SBATCH -J probe-failed-node-hP
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/probe_failed_node_hP_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/probe_failed_node_hP_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -q shared
#SBATCH -t 04:00:00
#
# The one SOLVE_FAILED node of the eps=1e-10 grid, re-solved with headroom.
#
#   sbatch GCE/scripts/probe_failed_node_hP.sbatch
#
# Shared QOS, 8 of a node's 256 cores: one worker per coupling point, so the
# wall time is the slowest single point, not their sum.  Worst case at
# PROBE_TIMEOUT=3600 is joint + retry + background = 3 h, hence the 4 h limit.
# Billing is 8/256 of a node, about 0.1 node-hours -- the allocation is down
# to ~215 node-hours and a whole node for this would not be worth it.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd /global/u1/g/gab97/secluded_DM

python GCE/scripts/probe_failed_node_hP.py 8
