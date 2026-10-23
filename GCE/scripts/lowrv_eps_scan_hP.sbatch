#!/bin/bash
#SBATCH -J lowrv-eps-scan
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/%x_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/%x_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -q shared
#SBATCH -t 04:00:00
#
# Why only eps=1e-9 fails in the rv < 0.2 corner.  252 solves at a 1200 s cap,
# so the worst case is well inside the wall on 64 cores.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd /global/u1/g/gab97/secluded_DM

python GCE/scripts/lowrv_eps_scan_hP.py 64
