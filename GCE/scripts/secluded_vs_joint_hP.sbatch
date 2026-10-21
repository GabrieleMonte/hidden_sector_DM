#!/bin/bash
#SBATCH -J secluded-vs-joint-hP
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/secluded_vs_joint_hP_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/secluded_vs_joint_hP_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -q shared
#SBATCH -t 02:00:00
#
# Does the secluded route reproduce the joint route at eps = 1e-9?
#
#   sbatch GCE/scripts/secluded_vs_joint_hP.sbatch
#
# 12 nodes x 2 routes = 24 sweeps, one per core.  Shared QOS, 32 of 256
# cores: about 0.25 node-hours.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd /global/u1/g/gab97/secluded_DM

python GCE/scripts/secluded_vs_joint_hP.py 32
