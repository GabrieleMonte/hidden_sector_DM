#!/bin/bash
#SBATCH -J relic-hP-sec
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_sec_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_sec_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 16
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH -q regular
#SBATCH -t 08:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#
# Higgs-portal SECLUDED only, on the 14-point lam axis.
#
#   sbatch GCE/scripts/relic_hP_secluded.sbatch
#
# Job 59713507 ran this slice on 8 lam points over [2e-5, 1] -- a factor 4.69
# per step -- and found 2114 of 6249 crossings, against 2856 on the earlier
# 8 points over [5e-3, 1], a factor 2.13 per step.  The 742 lost nodes are at
# the rv ~ 0.7 edge, where Omega peaks near the 0.12 target and resolution is
# what decides whether the crossing is bracketed.  N_LAM puts 14 points on the
# wide axis: a factor 2.30 per step, the density that produced 2856, without
# giving up the floor that supplies the 107 crossings below 5e-3.
#
# The 14-point axis does not match the stored 8-point one, so load_state
# starts this output fresh -- intended, and the 8-point grid is kept beside it
# as relic_grid_hP_eps_secluded_8pt_backup.npz.
#
# eps=1e-9 and eps=1e-10 keep 8 points and are untouched: their stored axes
# were checked bit-for-bit against this code while 59713507 was still running.
#
# Sharded by rv column like the main campaign, 16 ways: the slice took 9.4 h
# on one node at 8 points, so ~16 h at 14, but each shard holds ~130 nodes
# against 128 workers, so a shard is about one wave.  16 nodes also queues
# faster here than one does.

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

run_one () {
    local k=$1
    local log=$LOGS/relic_hP_sec_${SLURM_JOB_ID}_s${k}.out
    local rc=0
    {
        echo "[s$k] $(date '+%F %T') starting"
        srun --nodes=1 --ntasks=1 --cpus-per-task=256 --exclusive --cpu-bind=none \
             python "$DRIVER" --only secluded --shard "$k" --nshards "$NSHARDS" -n "$WORKERS"
        rc=$?
        echo "[s$k] $(date '+%F %T') finished rc=$rc"
    } >"$log" 2>&1
    return $rc
}

echo "Higgs portal secluded, 14 lam points, $SLURM_JOB_NUM_NODES nodes"
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
echo "merging secluded shards"
python "$DRIVER" --only secluded --nshards "$NSHARDS" --merge
mrc=$?

echo
if [ "$fail" -eq 0 ] && [ "$mrc" -eq 0 ]; then
    echo "secluded complete on the 14-point axis"
else
    echo "a shard or the merge failed -- shards checkpoint, so a resubmit resumes"
fi
exit $(( fail != 0 || mrc != 0 ))
