#!/bin/bash
#SBATCH -J relic-hP-slice
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_%x_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_%x_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 16
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH -q regular
#SBATCH -t 12:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#
# One Higgs-portal output, sharded 16 ways by rv column, one node per shard.
#
#   sbatch -J hP-secluded --export=ALL,SLICE=secluded GCE/scripts/relic_hP_slice.sbatch
#   sbatch -J hP-e09      --export=ALL,SLICE=1e-9     GCE/scripts/relic_hP_slice.sbatch
#   sbatch -J hP-e10      --export=ALL,SLICE=1e-10    GCE/scripts/relic_hP_slice.sbatch
#
# This replaces relic_hP.sbatch and relic_hP_secluded.sbatch.  The former ran
# secluded and eps=1e-9 unsharded on one node each: measured on 59713507, that
# left eps=1e-9 at 810 of 2083 nodes after nine hours, about 26 h of work, so
# it never reached the wall-clock limit and the merge never ran.  Every output
# needs the same sharding eps=1e-10 already had.
#
# A stride, not a block: the stiff nodes sit at high rv, so a block split would
# hand one shard most of the cost.  Each shard writes its own file; shards
# checkpoint per node, so a resubmit resumes where it stopped.
#
# The merge is NOT here.  It is a separate job submitted with
# --dependency=afterany, so it still runs when this one is cancelled or hits
# the limit -- which is what cost 59695521 and 59713507 their merges.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd /global/u1/g/gab97/secluded_DM

DRIVER=GCE/scripts/make_relic_grid_hP.py
LOGS=GCE/scripts/logs
WORKERS=128
NSHARDS=16
: "${SLICE:?set SLICE=secluded|1e-9|1e-10 via --export}"

run_one () {
    local k=$1
    local log=$LOGS/relic_hP_${SLICE}_${SLURM_JOB_ID}_s${k}.out
    local rc=0
    {
        echo "[$SLICE s$k] $(date '+%F %T') starting"
        srun --nodes=1 --ntasks=1 --cpus-per-task=256 --exclusive --cpu-bind=none \
             python "$DRIVER" --only "$SLICE" --shard "$k" --nshards "$NSHARDS" -n "$WORKERS"
        rc=$?
        echo "[$SLICE s$k] $(date '+%F %T') finished rc=$rc"
    } >"$log" 2>&1
    return $rc
}

echo "Higgs portal $SLICE, $NSHARDS shards on $SLURM_JOB_NUM_NODES nodes"
[ "${SLURM_JOB_NUM_NODES:-0}" -ge "$NSHARDS" ] || {
    echo "HALTED: need $NSHARDS nodes, have ${SLURM_JOB_NUM_NODES}." >&2
    exit 1; }

PIDS=""
for k in $(seq 0 $((NSHARDS - 1))); do
    run_one "$k" & PIDS="$PIDS $!:s$k"
done

fail=0
for entry in $PIDS; do
    pid=${entry%%:*}
    name=${entry#*:}
    wait "$pid"
    rc=$?
    printf '  %-6s rc=%s\n' "$name" "$rc"
    [ "$rc" -ne 0 ] && fail=1
done

echo
[ "$fail" -eq 0 ] && echo "$SLICE: all shards done" || echo "$SLICE: a shard failed"
exit "$fail"
