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
#   sbatch -q regular -t 48:00:00 GCE/scripts/relic_phase3.sbatch make_relic_grid_vP_xi_ini_2.5.py
#
# Single node by design: the relic grid has no --worker sharding; mp.Pool fans
# the (mX, rv) nodes over the node's 256 logical cores.  --requeue + the
# per-node checkpoint mask (CHECKPOINT_EVERY=8) make wall-clock kills free:
# a worst-case grid node is ~32 h (8 alphas x 1 retry x SOLVE_TIMEOUT=7200 s).
# Submit with -t 48:00:00, NOT 24 h: --requeue covers node failure and preemption,
# but a job killed at its time limit ends TIMEOUT and is not requeued, and the
# checkpoint mask only records COMPLETED nodes -- so a >24 h node would restart
# from zero in every 24 h slice and never finish.  48 h is the QOS maximum and
# costs nothing extra: the driver exits as soon as its queue drains.
# Run phase3_preflight.py and phase3_invalidate.py BEFORE the first job.

module load python
source "$HOME/cosmo_env/bin/activate"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd /global/u1/g/gab97/secluded_DM

SCRIPT=${1:?usage: sbatch [-q ... -t ...] relic_phase3.sbatch <make_relic_grid_vP*.py> [workers]}
WORKERS=${2:-128}   # 128 = physical cores.  256 logical workers roughly
                    # doubles per-solve latency, which pushes marginal nodes
                    # past SOLVE_TIMEOUT into a wasted 7200 s retry.

# The dry run is a hard gate, not a note.  The queue count must match the campaign
# doc (782 pass A, 460 pass B, 662 stage 5), but that count alone cannot catch a
# skipped Stage 2: ONLY_RV_COLUMNS produces it either way.  "starting fresh" is the
# real tell -- it means phase3_invalidate.py never ran (or the stored alpha_X axis
# drifted), and continuing would drop the finished low-rv columns on the next _save.
DRY=$(python GCE/scripts/"$SCRIPT" --dry-run 2>&1)
echo "$DRY"
if grep -q "starting fresh" <<<"$DRY"; then
  echo "ABORT: $SCRIPT would start fresh -- run GCE/scripts/phase3_invalidate.py first," >&2
  echo "       and check ALPHAX_LIMS against the axis stored in the .npz." >&2
  exit 1
fi

srun --ntasks-per-node=1 -c 256 --cpu-bind=none \
  python GCE/scripts/"$SCRIPT" -n "$WORKERS"
