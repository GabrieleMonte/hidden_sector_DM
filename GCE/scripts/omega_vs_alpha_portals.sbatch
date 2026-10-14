#!/bin/bash
#SBATCH -J omega-alpha-portals
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/omega_vs_alpha_portals_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/omega_vs_alpha_portals_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH -q debug
#SBATCH -t 00:30:00
#
# Do BL / LiLj / B need rebuilding on the fixed crossing, and did they ever
# get the Phase-3 follow-up?
#
#   sbatch GCE/scripts/omega_vs_alpha_portals.sbatch
#
# Debug QOS: 3 portals x 16 nodes x 8 alphas x (1 secluded + 2 joint) = 1152
# solves on 128 workers, 9 waves.  A login-node smoke test measured 8-10 s per
# solve at rv=0.12; the stiff high-rv end is far slower, so PROBE_TIMEOUT caps
# a single solve at 600 s and the sampling is even rather than stiff-weighted.
#
# One node, not sixteen: the 1-node regular queue is ~17,700 deep, so debug is
# the only way a single-node job starts today.  The probe is small enough to
# fit, which is the point of sampling instead of re-running the grids.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd /global/u1/g/gab97/secluded_DM

mkdir -p GCE/output/diag

srun --ntasks=1 --cpus-per-task=256 --cpu-bind=none \
     python GCE/scripts/omega_vs_alpha_portals.py 128
