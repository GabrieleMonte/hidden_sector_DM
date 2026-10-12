#!/bin/bash
#SBATCH -J relic-hP
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 18
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH -q regular
#SBATCH -t 12:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#
# Higgs-portal relic grids at xi_ini = 1, Phase 3 on.
#
#   sbatch GCE/scripts/relic_hP.sbatch
#
# EIGHTEEN NODES, BECAUSE ONE SLICE DOMINATES.  Measured on the first run
# (59694081): secluded 2.36 s/it and eps=1e-9 1.30 s/it, both done inside an
# hour on one node each, but eps=1e-10 ran at 25 s/it -- about 87 minutes of
# CPU per node, 25x the others -- which alone would take ~14 h.
#
# So eps=1e-10 is split into NSHARDS concurrent shards, one node each, by rv
# column.  A stride, not a block: the stiff nodes sit at high rv, so a block
# split would hand one shard nearly all the cost.  Each shard writes its own
# file and the merge at the end folds them into the canonical grid, keeping
# whatever a previous run had already finished.
#
# The floor is one node's wall time -- a single stiff solve can run to
# SOLVE_TIMEOUT = 7200 s -- so more than ~16 shards buys nothing.  18 nodes
# also sits in the least contended queue bucket (the campaign measured 87
# pending jobs at 17-32 nodes against 421 at 9-16).
#
# Gamma_phi at epsX = 1e-9 is x0.002 the vector portal's, so phi is ~480x
# longer-lived at the same epsX: more entropy injection, and Phase 3 engaging
# far more often than in the vector campaigns.  That is why 1e-10 is so stiff.

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
SHARDED=1e-10                     # the slice that needs splitting

# run_one <label> <extra args...> -- one srun on one node.
# Ends in `return`, never `exit`: these are backgrounded and an exit would
# kill the caller's subshell.
run_one () {
    local label=$1
    shift
    local log=$LOGS/relic_hP_${SLURM_JOB_ID}_${label}.out
    local rc=0
    {
        echo "[$label] $(date '+%F %T') starting: $*"
        srun --nodes=1 --ntasks=1 --cpus-per-task=256 --exclusive --cpu-bind=none \
             python "$DRIVER" "$@" -n "$WORKERS"
        rc=$?
        echo "[$label] $(date '+%F %T') finished rc=$rc"
    } >"$log" 2>&1
    return $rc
}

echo "Higgs portal, $SLURM_JOB_NUM_NODES nodes: $SLURM_JOB_NODELIST"
echo "eps=$SHARDED split into $NSHARDS shards by rv column"
echo "per-slice logs: $LOGS/relic_hP_${SLURM_JOB_ID}_<label>.out"
[ "${SLURM_JOB_NUM_NODES:-0}" -ge $((NSHARDS + 2)) ] || {
    echo "HALTED: need $((NSHARDS + 2)) nodes, have ${SLURM_JOB_NUM_NODES}." >&2
    exit 1; }

PIDS=""
run_one secluded --only secluded & PIDS="$PIDS $!:secluded"
run_one 1e-9     --only 1e-9     & PIDS="$PIDS $!:1e-9"
for k in $(seq 0 $((NSHARDS - 1))); do
    run_one "${SHARDED}_s${k}" --only "$SHARDED" --shard "$k" --nshards "$NSHARDS" &
    PIDS="$PIDS $!:${SHARDED}_s${k}"
done

fail=0
for entry in $PIDS; do
    pid=${entry%%:*}
    name=${entry#*:}
    wait "$pid"
    rc=$?
    printf '  %-16s rc=%s\n' "$name" "$rc"
    [ "$rc" -ne 0 ] && fail=1
done

echo
echo "merging $SHARDED shards"
python "$DRIVER" --only "$SHARDED" --nshards "$NSHARDS" --merge
mrc=$?

echo
if [ "$fail" -eq 0 ] && [ "$mrc" -eq 0 ]; then
    echo "all slices complete"
else
    echo "a slice or the merge failed -- see its log; shards checkpoint per"
    echo "node, so a resubmit resumes where it stopped"
fi
exit $(( fail != 0 || mrc != 0 ))
