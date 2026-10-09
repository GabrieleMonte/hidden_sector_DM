#!/bin/bash
#SBATCH -J relic-p3
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_p3_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_p3_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH --open-mode=append
#SBATCH --requeue
#
# Phase-3 rerun campaign, Stage 3/4 (docs/phase3_rerun_campaign.md).
# One slice per job; -q and -t on the command line:
#
#   sbatch -q regular -t 24:00:00 GCE/scripts/relic_phase3.sbatch make_relic_grid_vP_xi_ini_2.5.py
#
# Single node by design: the relic grid has no --worker sharding; mp.Pool fans
# the (mX, rv) nodes over the node's 256 logical cores.  --requeue + the
# per-node checkpoint mask (CHECKPOINT_EVERY=8) make wall-clock kills free:
# a worst-case grid node is ~32 h (8 alphas x 1 retry x SOLVE_TIMEOUT=7200 s),
# so a 24 h slice may requeue once or twice near the end -- that is expected.
# Run phase3_preflight.py and phase3_invalidate.py BEFORE the first job.

module load python
source "$HOME/cosmo_env/bin/activate"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd /global/u1/g/gab97/secluded_DM

SCRIPT=${1:?usage: sbatch [-q ... -t ...] relic_phase3.sbatch <make_relic_grid_vP*.py>}

# The queue count printed here must match the campaign doc (782 for pass A,
# 460 for pass B) -- a surprise here means the invalidation step was skipped.
python GCE/scripts/"$SCRIPT" --dry-run

srun --ntasks-per-node=1 -c 256 --cpu-bind=none \
  python GCE/scripts/"$SCRIPT" -n 256
