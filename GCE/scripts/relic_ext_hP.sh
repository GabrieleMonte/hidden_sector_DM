#!/bin/bash
#SBATCH -J relic_ext_hP
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_ext_hP_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_ext_hP_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH --open-mode=append
#SBATCH --requeue
# -N, -q and -t are supplied on the sbatch command line.

module load python
source "$HOME/cosmo_env/bin/activate"

# The solves are single-threaded scipy ODE integrations; the parallelism is
# mp.Pool.  Pin any stray BLAS threading to 1 so the workers cannot
# oversubscribe the node.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd /global/u1/g/gab97/secluded_DM

# The three outputs write three different files, so they run as concurrent
# processes rather than one after another.  66 nodes per output caps the useful
# worker count at 66 each; 3 x 66 = 198 fits the node's 256 hardware threads.
# make_relic_grid checkpoints after every node, so a wall hit loses at most one.
NCORES=66
for ONLY in secluded 1e-7 1e-8; do
    python GCE/scripts/make_relic_grid_hP_ext.py --only "$ONLY" -n "$NCORES" &
done
wait
