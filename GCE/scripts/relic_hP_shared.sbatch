#!/bin/bash
#SBATCH -J relic-hP-shared
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_%x_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_%x_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH -q shared
#SBATCH -t 03:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#
# One Higgs-portal output, unsharded, on half a node of the shared QOS.
#
#   sbatch -J hP-e08 --export=ALL,SLICE=1e-8 GCE/scripts/relic_hP_shared.sbatch
#   sbatch -J hP-e07 --export=ALL,SLICE=1e-7 GCE/scripts/relic_hP_shared.sbatch
#
# Use this, not relic_hP_slice.sbatch, for an output whose nodes are cheap.
# Measured 2026-10-20 at eps=1e-8 over 24 grid nodes on the production axis:
# median 87.4 s, mean 124.8 s, max 1043.8 s at rv=0.950, so the whole grid is
# about 72 core-hours.  On 128 logical cores that is ~1 h wall, and the shared
# QOS bills 128/256 of a node, so ~0.5 node-hours.
#
# relic_hP_slice.sbatch asks for 16 regular nodes and 2048 workers against
# 2083 tasks, which is a SINGLE wave: the wall is set by the slowest node
# while 2047 workers idle, and the whole node-count is billed for it.  That
# cost eps=1e-9 19.8 node-hours for the same ~72 core-hours of work.  Keep it
# for outputs whose per-node cost is hours, not minutes.
#
# No sharding and no merge: one process owns the grid and checkpoints every
# node (CHECKPOINT_EVERY = 1), so a requeue resumes where it stopped.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd /global/u1/g/gab97/secluded_DM

: "${SLICE:?set SLICE=1e-7|1e-8|1e-9|1e-10|secluded via --export}"

python GCE/scripts/make_relic_grid_hP.py --only "$SLICE" -n 128
