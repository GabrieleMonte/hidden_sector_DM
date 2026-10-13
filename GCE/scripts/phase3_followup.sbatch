#!/bin/bash
#SBATCH -J relic-p3-followup
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_p3fu_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_p3fu_%j.err
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
# Phase-3 follow-up (docs/phase3_followup.md), unattended.  Submit ONE job:
#
#   sbatch GCE/scripts/phase3_followup.sbatch
#
# Rehearse the control flow without spending anything (no preflight, no
# invalidation, no srun; the drivers' read-only --dry-run gates still run):
#
#   PHASE3_FOLLOWUP_DRYRUN=1 bash GCE/scripts/phase3_followup.sbatch
#
# What it does:
#   Stage 0  provenance; phase3_followup_variants.py --check
#   Stage 1  preflight F1-F3 on the head node, ~10 min        } skipped once the
#   Stage 2  phase3_followup_invalidate.py (writes its marker) } marker exists
#   Stage 3  four slices on three nodes, each its own .npz:
#              node 1  xi=1   eps=1e-9   cols  0-66  2083 nodes
#              node 2  xi=2.5 eps=1e-10  cols  0-56  1623 nodes
#              node 3  xi=2.5 eps=1e-9   cols  0-49  1301 nodes,
#                      then xi=1 eps=1e-10 cols 57-66  460 nodes
#   Stage 4  phase3_followup_report.py: new vs old crossings (read-only)
#
# REQUEUE-SAFE.  The invalidation must run exactly once (a second run would
# clear nodes this job already re-solved), so it leaves
# GCE/output/relic/phase3_followup_invalidated.txt and a requeued job skips
# Stages 1-2 and resumes every slice from its checkpoint mask.
# 48 h: a job killed at its time limit ends TIMEOUT and is NOT requeued, and the
# slices exit as soon as their queues drain, so the long limit costs nothing.

DRYRUN=${PHASE3_FOLLOWUP_DRYRUN:-0}
REPO=${PHASE3_FOLLOWUP_REPO:-/global/u1/g/gab97/secluded_DM}

if [ "$DRYRUN" -eq 0 ]; then
    module load python
    source "$HOME/cosmo_env/bin/activate"
fi
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd "$REPO" || exit 1

SCRIPTS=GCE/scripts
LOGS=$SCRIPTS/logs
MARKER=GCE/output/relic/phase3_followup_invalidated.txt
JOBID=${SLURM_JOB_ID:-dryrun}
WORKERS=128                        # physical cores; see relic_phase3.sbatch for why not 256
mkdir -p "$LOGS"

banner () { printf '\n========== %s ==========\n' "$1"; }
die     () { echo; echo "FOLLOW-UP HALTED: $1" >&2; exit 1; }

# run_slice <label> <driver.py> <expected nodes> -- one slice on one node.
# Gate: the driver's own --dry-run must not say "starting fresh" (the stored
# axis/eps do not match, and the first _save would drop every other column) and
# must queue at most the slice's node count (more means the wrong columns).
# Ends in `return`, never `exit`: node 3 chains two slices in one subshell.
run_slice () {
    local label=$1 driver=$2 expect=$3
    local log=$LOGS/relic_p3fu_${JOBID}_${label}.out
    local rc=0 dry left
    {
        echo "[$label] $(date '+%F %T') starting, driver=$driver"
        dry=$(python "$SCRIPTS/$driver" --dry-run 2>&1)
        echo "$dry"
        left=$(sed -nE 's/.* ([0-9]+) of [0-9]+ nodes still to solve.*/\1/p' <<<"$dry")
        if grep -q "starting fresh" <<<"$dry"; then
            echo "[$label] ABORT: the driver would start fresh -- refusing to overwrite."
            rc=1
        elif [ -z "$left" ] || [ "$left" -gt "$expect" ]; then
            echo "[$label] ABORT: queue '${left}' is not within the slice's $expect nodes."
            rc=1
        elif [ "$DRYRUN" -eq 1 ]; then
            echo "[$label] DRYRUN: queue $left/$expect; would srun -N1 -n1 -c256 python $driver -n $WORKERS"
        else
            echo "[$label] queue $left of $expect nodes"
            srun --nodes=1 --ntasks=1 --cpus-per-task=256 --exclusive --cpu-bind=none \
                 python "$SCRIPTS/$driver" -n "$WORKERS"
            rc=$?
            echo "[$label] $(date '+%F %T') finished rc=$rc"
        fi
    } >"$log" 2>&1
    return $rc
}

banner "Stage 0: provenance"
git rev-parse HEAD
git status --porcelain "$SCRIPTS" GCE/output/relic | sed 's/^/  dirty: /'
python "$SCRIPTS/phase3_followup_variants.py" --check \
  || die "driver variants are stale -- rerun phase3_followup_variants.py and commit"
if [ "$DRYRUN" -eq 0 ]; then
    echo "  allocation: $SLURM_JOB_NUM_NODES nodes -- $SLURM_JOB_NODELIST"
    [ "${SLURM_JOB_NUM_NODES:-0}" -ge 3 ] \
      || die "this job holds ${SLURM_JOB_NUM_NODES} nodes; the follow-up runs 3 concurrently."
fi

if [ -f "$MARKER" ]; then
    banner "Stages 1-2: already done -- resuming"
    cat "$MARKER"
else
    banner "Stage 1: preflight"
    if [ "$DRYRUN" -eq 1 ]; then
        echo "  DRYRUN: skipping preflight; assuming GO"
    else
        python "$SCRIPTS/phase3_followup_preflight.py" -n 32 \
          2>&1 | tee "$LOGS/relic_p3fu_preflight_${JOBID}.log"
        [ "${PIPESTATUS[0]}" -eq 0 ] || die "preflight did not pass -- read the F1-F3 lines
       in $LOGS/relic_p3fu_preflight_${JOBID}.log.  Nothing was invalidated."
    fi

    banner "Stage 2: invalidate"
    if [ "$DRYRUN" -eq 1 ]; then
        echo "  DRYRUN: skipping phase3_followup_invalidate.py"
    else
        python "$SCRIPTS/phase3_followup_invalidate.py" || die "invalidation failed"
    fi
fi

banner "Stage 3: four slices on three nodes"
echo "  per-slice logs: $LOGS/relic_p3fu_${JOBID}_<label>.out"
declare -A RC
run_slice xi1_e09  make_relic_grid_vP_fu_e09.py            2083 & P1=$!
run_slice xi25_e10 make_relic_grid_vP_xi_ini_2.5_fu_e10.py 1623 & P2=$!
( run_slice xi25_e09 make_relic_grid_vP_xi_ini_2.5_fu_e09.py 1301; r1=$?
  run_slice xi1_e10  make_relic_grid_vP_fu_e10.py             460; r2=$?
  exit $(( r1 != 0 || r2 != 0 )) ) & P3=$!
wait $P1; RC[xi1_e09]=$?
wait $P2; RC[xi25_e10]=$?
wait $P3; RC[xi25_e09+xi1_e10]=$?

banner "Stage 4: report"
fail=0
for k in "${!RC[@]}"; do
    printf '  %-20s rc=%s\n' "$k" "${RC[$k]}"
    [ "${RC[$k]}" -ne 0 ] && fail=1
done
if [ "$DRYRUN" -eq 0 ]; then
    python "$SCRIPTS/phase3_followup_report.py" || echo "  (report failed -- the grids are unaffected)"
fi
if [ "$fail" -eq 0 ]; then
    echo; echo "Follow-up complete.  Commit GCE/output/relic (incl. the marker) on cluster-relic-vP."
else
    echo; echo "A slice returned non-zero -- read its log.  Resubmitting resumes from the checkpoints."
fi
exit $fail
