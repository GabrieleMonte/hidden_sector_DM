#!/bin/bash
#SBATCH -J relic-bf-xi1
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_bf_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_bf_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 2
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH -q regular
#SBATCH -t 48:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#
# xi_ini = 1 backfill: the 1830 nodes left after the Phase-3 campaign
# (813de606bc).  Both slices run on the LEGACY path -- the drivers are in their
# restored configuration (aebd7ce518), PHASE3_KW = {} -- because every node here
# sits below rv column 57 and must match the trusted neighbours beside it.
#
#   secluded   1030 nodes  the midpoint columns the 34->67 embed left empty;
#                          brings xi=1 level with the xi=0.5 and 2.5 slices
#   eps=1e-10   800 nodes  the midpoint columns below col 57 (rv 0.052-0.582)
#
# Different output files, so they run concurrently, one node each.  Two nodes
# rather than one sequential pass: it halves the wall clock and a 1-node job
# draws from the most contended bucket on the machine.
#
#   sbatch GCE/scripts/relic_backfill_xi1.sbatch

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd /global/u1/g/gab97/secluded_DM

LOGS=GCE/scripts/logs
WORKERS=128

run_slice () {
    local label=$1 driver=$2
    local log=$LOGS/relic_bf_${SLURM_JOB_ID}_${label}.out
    local rc=0 dry
    {
        echo "[$label] $(date '+%F %T') starting $driver"
        dry=$(python "GCE/scripts/$driver" --dry-run 2>&1)
        echo "$dry"
        if grep -q "starting fresh" <<<"$dry"; then
            echo "[$label] ABORT: would start fresh -- the stored alpha_X axis"
            echo "         does not match ALPHAX_LIMS.  Refusing to overwrite."
            rc=1
        else
            srun --nodes=1 --ntasks=1 --cpus-per-task=256 --exclusive --cpu-bind=none \
                 python "GCE/scripts/$driver" -n "$WORKERS"
            rc=$?
            echo "[$label] $(date '+%F %T') finished rc=$rc"
        fi
    } >"$log" 2>&1
    return $rc
}

echo "xi_ini=1 backfill, 2 nodes: $SLURM_JOB_NODELIST"
echo "per-slice logs: $LOGS/relic_bf_${SLURM_JOB_ID}_<label>.out"

run_slice secluded make_relic_grid_vP_bf_sec.py & P_S=$!
run_slice eps1e-10 make_relic_grid_vP_bf_e10.py & P_E=$!
wait $P_S; RC_S=$?
wait $P_E; RC_E=$?

echo
echo "secluded rc=$RC_S   eps1e-10 rc=$RC_E"
[ $RC_S -eq 0 ] && [ $RC_E -eq 0 ] && echo "backfill complete" || echo "a slice failed -- see its log"
exit $(( RC_S != 0 || RC_E != 0 ))
