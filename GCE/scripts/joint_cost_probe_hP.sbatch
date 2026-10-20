#!/bin/bash
#SBATCH -J joint-cost-probe-hP
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/joint_cost_probe_hP_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/joint_cost_probe_hP_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 96
#SBATCH -q shared
#SBATCH -t 01:15:00
#
# Attribute the eps=1e-10 joint cost before resubmitting the rebuild.
#
#   sbatch GCE/scripts/joint_cost_probe_hP.sbatch
#
# 4 nodes x 2 eps x 6 lam x 2 rtol = 96 half-split solves, one per core, each
# capped at 1800 s.  Shared QOS, 96 of 256 cores: about 0.5 node-hours.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd /global/u1/g/gab97/secluded_DM

python GCE/scripts/joint_cost_probe_hP.py 96
