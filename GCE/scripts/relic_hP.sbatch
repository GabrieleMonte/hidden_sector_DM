#!/bin/bash
#SBATCH -J relic-hP
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 3
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH -q regular
#SBATCH -t 48:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#
# Higgs-portal relic grids at xi_ini = 1, Phase 3 on.
#
#   sbatch GCE/scripts/relic_hP.sbatch
#
# Three slices, one per node: secluded, eps=1e-9, eps=1e-10.  They write
# different files, so the wall clock is the slowest slice, not the sum.  These
# grids are new, so each starts fresh; make_relic_grid.py checkpoints per node,
# so a resubmit resumes where it stopped.
#
# Gamma_phi at epsX = 1e-9 is x0.002 the vector portal's, so phi is ~480x
# longer-lived at the same epsX: expect more entropy injection and Phase 3
# engaging more often than in the vector campaigns, hence the 48 h limit.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd /global/u1/g/gab97/secluded_DM

LOGS=GCE/scripts/logs
WORKERS=128

run_slice () {
    local only=$1
    local log=$LOGS/relic_hP_${SLURM_JOB_ID}_${only}.out
    local rc=0
    # Ends in `return`, never `exit`: these are backgrounded and an exit would
    # kill the caller's subshell.
    {
        echo "[$only] $(date '+%F %T') starting --only $only"
        srun --nodes=1 --ntasks=1 --cpus-per-task=256 --exclusive --cpu-bind=none \
             python GCE/scripts/make_relic_grid_hP.py --only "$only" -n "$WORKERS"
        rc=$?
        echo "[$only] $(date '+%F %T') finished rc=$rc"
    } >"$log" 2>&1
    return $rc
}

echo "Higgs portal, 3 nodes: $SLURM_JOB_NODELIST"
echo "per-slice logs: $LOGS/relic_hP_${SLURM_JOB_ID}_<output>.out"
python GCE/scripts/make_relic_grid_hP.py --dry-run

run_slice secluded & P_S=$!
run_slice 1e-9    & P_9=$!
run_slice 1e-10   & P_10=$!
wait $P_S;  RC_S=$?
wait $P_9;  RC_9=$?
wait $P_10; RC_10=$?

echo
echo "secluded rc=$RC_S   eps1e-9 rc=$RC_9   eps1e-10 rc=$RC_10"
if [ $RC_S -eq 0 ] && [ $RC_9 -eq 0 ] && [ $RC_10 -eq 0 ]; then
    echo "all slices complete"
else
    echo "a slice failed -- see its log; a resubmit resumes from the checkpoints"
fi
exit $(( RC_S != 0 || RC_9 != 0 || RC_10 != 0 ))
