#!/bin/bash
#SBATCH -J relic-portals
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_portals_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_portals_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 9
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH -q regular
#SBATCH -t 48:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#
# BL, LiLj(mu-tau) and Baryon relic grids at xi_ini = 1, mediators at 10 TeV.
# Submit this ONE job and log out:
#
#   sbatch GCE/scripts/portal_pipeline.sbatch
#
# Rehearse the control flow without spending anything:
#
#   PORTAL_PIPELINE_DRYRUN=1 bash GCE/scripts/portal_pipeline.sbatch
#
# ONE QUEUE WAIT.  The preflight gates the slices, so running it as its own job
# would put a second full queue wait (measured 1.7-3.1 h) on the critical path,
# and sbatch is refused for -q interactive here.  So this job holds the whole
# allocation: it preflights on the head node, then runs every slice in place.
# The preflight leaves 8 nodes idle for ~2 h, ~16 node-hours out of ~2,100
# remaining, which is what buys the single wait.
#
# NINE SLICES, ONE PER NODE.  3 portals x 3 outputs (secluded, eps=1e-9,
# eps=1e-10).  All nine write different .npz files, so they are safe to run at
# once and the wall clock is the slowest slice, not the sum.  A 9-node request
# also sits in a far shorter queue than nine 1-node ones (the campaign measured
# 17,739 pending 1-node jobs against 421 in the 9-16 bucket, and this account's
# own 1-node waits range to 15.9 h against ~3 h for multi-node).
#
# Phase 3 is ON (the driver default): these are new grids with no legacy
# neighbours to match, and the 2026-10 campaign measured the conserved-S
# approximation biasing alpha_relic by x1.07-1.20 at rv >= 0.465.

DRYRUN=${PORTAL_PIPELINE_DRYRUN:-0}

if [ "$DRYRUN" -eq 0 ]; then
    module load python
    source "$HOME/cosmo_env/bin/activate"
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
fi

cd /global/u1/g/gab97/secluded_DM || exit 1

SCRIPTS=GCE/scripts
LOGS=$SCRIPTS/logs
JOBID=${SLURM_JOB_ID:-dryrun}
WORKERS=128                      # physical cores; see relic_portal.sbatch
PORTALS="BL LiLj B"
PRE_LOG=$LOGS/relic_portals_preflight_${JOBID}.log

banner () { printf '\n========== %s ==========\n' "$1"; }
die     () { echo; echo "PIPELINE HALTED: $1" >&2; exit 1; }

# run_slice <portal> <only> -- one output of one portal, on one node.
run_slice () {
    local portal=$1 only=$2
    local label=${portal}_${only}
    local log=$LOGS/relic_portals_${JOBID}_${label}.out
    local rc=0 dry
    # Must end in `return`, never `exit`: callers background these and an exit
    # would kill the caller's subshell.
    {
        echo "[$label] $(date '+%F %T') starting --portal $portal --only $only"
        if [ "$DRYRUN" -eq 1 ]; then
            echo "[$label] DRYRUN: would srun make_relic_grid.py --portal $portal" \
                 "--xi-ini 1 --only $only -n $WORKERS"
        else
            dry=$(python "$SCRIPTS/make_relic_grid.py" --portal "$portal" \
                  --xi-ini 1 --only "$only" --dry-run 2>&1)
            echo "$dry"
            if grep -qE "was built with Phase 3|starting fresh" <<<"$dry"; then
                # "starting fresh" is expected for a brand-new portal; only the
                # Phase-3 mismatch is fatal, and load_state exits non-zero on it.
                grep -q "was built with Phase 3" <<<"$dry" && {
                    echo "[$label] ABORT: Phase-3 setting disagrees with the stored grid."
                    rc=1; }
            fi
            if [ $rc -eq 0 ]; then
                srun --nodes=1 --ntasks=1 --cpus-per-task=256 --exclusive \
                     --cpu-bind=none \
                     python "$SCRIPTS/make_relic_grid.py" --portal "$portal" \
                     --xi-ini 1 --only "$only" -n "$WORKERS"
                rc=$?
                echo "[$label] $(date '+%F %T') finished rc=$rc"
            fi
        fi
    } >"$log" 2>&1
    return $rc
}

banner "Stage 0: provenance"
git rev-parse HEAD
git status --porcelain "$SCRIPTS" GCE/output/relic | sed 's/^/  dirty: /'
if [ "$DRYRUN" -eq 0 ]; then
    echo "  allocation: $SLURM_JOB_NUM_NODES nodes -- $SLURM_JOB_NODELIST"
    [ "${SLURM_JOB_NUM_NODES:-0}" -ge 9 ] \
      || die "this job holds ${SLURM_JOB_NUM_NODES} nodes but runs 9 slices at once.
       Resubmit with -N 9."
fi

# ---- Stage 1: preflight ------------------------------------------------
# G2 is the gate: at the worst corner (rv=0.950, the mX where the vP grid
# peaks) the relic crossing must still be bracketed below alpha_X = 1.  A
# portal whose Gamma_Y is 14x the vector portal's dilutes less, so Omega rises
# and the crossing moves UP -- if it leaves the axis, production would write
# NaN exactly where the grid matters most.
banner "Stage 1: portal preflight (BL, LiLj, B)"
if [ "$DRYRUN" -eq 1 ]; then
    echo "  DRYRUN: skipping preflight; assuming every portal is GO"
    : >"$PRE_LOG"
    for p in $PORTALS; do echo "  $p    production: GO" >>"$PRE_LOG"; done
else
    python "$SCRIPTS/portal_preflight.py" --portals "BL,LiLj,B" -n 48 \
        2>&1 | tee "$PRE_LOG"
fi

GO_LIST=""
for p in $PORTALS; do
    if grep -qE "^  $p +production: GO" "$PRE_LOG"; then
        GO_LIST="$GO_LIST $p"
    else
        echo "SKIPPING $p: preflight did not return GO (see $PRE_LOG)."
    fi
done
[ -n "$GO_LIST" ] || die "no portal passed the preflight.  Nothing was run.
       A NO-GO means the crossing left the alpha_X axis at the worst corner;
       raise ALPHAX_LIMS in make_relic_grid.py and resubmit."

# ---- Stage 2: every slice at once, one per node ------------------------
banner "Stage 2: launch slices for:$GO_LIST"
declare -A RC
PIDS=""
for p in $GO_LIST; do
    for o in secluded 1e-9 1e-10; do
        run_slice "$p" "$o" & PIDS="$PIDS $!:${p}_${o}"
        echo "  ${p}_${o}  2083 nodes"
    done
done
echo "  per-slice logs: $LOGS/relic_portals_${JOBID}_<portal>_<output>.out"

for entry in $PIDS; do
    pid=${entry%%:*}; name=${entry#*:}
    wait "$pid"; RC[$name]=$?
done

banner "Finished"
fail=0
for k in "${!RC[@]}"; do
    printf '  %-16s rc=%s\n' "$k" "${RC[$k]}"
    [ "${RC[$k]}" -ne 0 ] && fail=1
done
if [ "$fail" -eq 0 ]; then
    echo; echo "All slices completed.  New grids in GCE/output/relic."
else
    echo; echo "A slice returned non-zero -- check its log before trusting its .npz."
    echo "Slices checkpoint per node, so resubmitting resumes where it stopped."
fi
exit $fail
