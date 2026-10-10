#!/bin/bash
#SBATCH -J relic-portal
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_portal_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_portal_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH --open-mode=append
#SBATCH --requeue
#
# One (portal, xi_ini) slice of the shared relic grid.  -q and -t go on the
# command line; everything after the script name is passed to the driver:
#
#   sbatch -q regular -t 48:00:00 GCE/scripts/relic_portal.sbatch --portal BL
#   sbatch -q regular -t 48:00:00 GCE/scripts/relic_portal.sbatch --portal LiLj --xi-ini 2.5
#
# -t 48:00:00, not 24 h: a worst-case node is ~32 h (8 alphas x 1 retry x
# SOLVE_TIMEOUT=7200 s), --requeue does not cover a walltime kill, and the
# checkpoint records only COMPLETED nodes -- so a >24 h node would restart from
# zero in every 24 h slice.  48 h is the QOS maximum and costs nothing extra:
# the driver exits as soon as its queue drains.
#
# Single node by design: the grid has no --worker sharding, so mp.Pool fans the
# (mX, rv) nodes over one node's cores.  Different (portal, xi_ini) slices write
# different files, so run them as concurrent job steps in one allocation when
# you want several -- a multi-node request also clears the queue far faster than
# several 1-node ones (see 813de606bc).

module load python
source "$HOME/cosmo_env/bin/activate"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd /global/u1/g/gab97/secluded_DM

WORKERS=${WORKERS:-128}   # physical cores; 256 logical workers roughly double
                          # per-solve latency and push marginal nodes past
                          # SOLVE_TIMEOUT into a wasted 7200 s retry.

[ $# -gt 0 ] || { echo "usage: sbatch [-q .. -t ..] relic_portal.sbatch --portal <vP|BL|LiLj|B> [--xi-ini X] [...]" >&2; exit 1; }

# The dry run is a gate: it prints the portal line, what each output reuses and
# the queue.  load_state aborts here if the stored grid was built under the
# other Phase-3 setting, rather than mixing two physics into one file.
DRY=$(python GCE/scripts/make_relic_grid.py "$@" --dry-run 2>&1) || {
    echo "$DRY"; echo "ABORT: dry run refused -- see the message above." >&2; exit 1; }
echo "$DRY"

srun --ntasks-per-node=1 -c 256 --cpu-bind=none \
  python GCE/scripts/make_relic_grid.py "$@" -n "$WORKERS"
