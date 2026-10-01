#!/bin/bash
#SBATCH -J chanspec
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/chanspec_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/chanspec_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH --open-mode=append
#SBATCH --requeue
# -N, -q and -t are supplied on the sbatch command line.

module load python
source "$HOME/cosmo_env/bin/activate"

# Pythia is single-threaded per process; the parallelism is mp.Pool.  Pin any
# stray BLAS threading to 1 so 256 workers cannot oversubscribe the node.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Stage on $PSCRATCH: compute nodes reach /global/homes over DVS, and this run
# writes ~157k small .npz files.  Copied back to the repo when the grid is done.
export CHANNEL_SPECTRA_DIR=$PSCRATCH/channel_spectra
mkdir -p "$CHANNEL_SPECTRA_DIR"

cd /global/u1/g/gab97/secluded_DM

# One task per node, each forking a 256-way pool over the node's logical cores.
# Rank w of N takes grid columns w, w+N, ... -- the stride sharding the script
# already implements, so every rank gets a mix of cheap and expensive points.
srun --ntasks-per-node=1 -c 256 --cpu-bind=none bash -c \
  'python GCE/scripts/make_channel_spectra.py -n 256 \
       --worker $SLURM_PROCID --nworkers '"$SLURM_NTASKS"
