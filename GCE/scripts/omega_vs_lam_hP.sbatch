#!/bin/bash
#SBATCH -J omega-lam-hP
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/omega_vs_lam_hP_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/omega_vs_lam_hP_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH -q debug
#SBATCH -t 00:30:00
#
# Omega(lam) at the Higgs-portal secluded nodes that have no crossing.
#
#   sbatch GCE/scripts/omega_vs_lam_hP.sbatch
#
# Debug QOS: 20 nodes x 15 lam = 300 solves on 128 workers.  A login-node
# smoke test measured ~95 s per (node, lam) -- these are three background
# solves each, not one -- so ~3 waves of 128 is well inside 30 min.
# PROBE_TIMEOUT caps any single solve at 600 s so one stiff outlier cannot
# run the allocation out.
#
# Reads a snapshot, not the live grid: job 59713507 is still checkpointing
# relic_grid_hP_eps_secluded.npz with np.savez, which rewrites in place.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd /global/u1/g/gab97/secluded_DM

DIAG=GCE/output/diag
mkdir -p "$DIAG"
cp GCE/output/relic/relic_grid_hP_eps_secluded.npz \
   "$DIAG/relic_grid_hP_eps_secluded_snapshot.npz"

srun --ntasks=1 --cpus-per-task=256 --cpu-bind=none \
     python GCE/scripts/omega_vs_lam_hP.py 128
