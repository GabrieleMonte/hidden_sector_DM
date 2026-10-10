#!/bin/bash
#SBATCH -J relic-p3-pipeline
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_p3_pipeline_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_p3_pipeline_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 6
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH -q regular
#SBATCH -t 48:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#
# Phase-3 rerun campaign, Stages 1-5, unattended (docs/phase3_rerun_campaign.md).
# Submit this ONE job and log out:
#
#   sbatch GCE/scripts/phase3_pipeline.sbatch
#   sbatch --mail-type=END,FAIL --mail-user=<you> GCE/scripts/phase3_pipeline.sbatch
#
# Rehearse the control flow without spending anything (skips the preflight, the
# invalidation and the in-place slice; prints the sbatch lines instead of running
# them):
#
#   PHASE3_PIPELINE_DRYRUN=1 bash GCE/scripts/phase3_pipeline.sbatch
#
# WHY ONE BIG JOB INSTEAD OF A PREFLIGHT JOB PLUS PRODUCTION JOBS.  The preflight
# gates everything, so a separate preflight job would put a second full queue wait
# (measured 1.7-3.1 h on regular_1) on the critical path.  -q interactive would
# have avoided it, but NERSC rejects sbatch to that QOS ("Cannot submit batch jobs
# to interactive_ss11"), and -q shared did not start promptly either.  So this job
# holds one node for the whole campaign front end: it runs the preflight itself,
# and on GO it submits the other five slices and then runs the pass-A xi=2.5 slice
# in place.  The campaign therefore pays exactly ONE queue wait.  The preflight
# leaves ~95 of 128 cores idle for ~2 h, which costs ~2 node-hours out of ~2,200
# remaining -- a cheap trade for cutting hours off the wall clock.
#
# What this job does, in order:
#   Stage 1  preflight P1/P3/P4 (pooled) with P2 running concurrently on 1 core
#   Stage 2  phase3_invalidate.py, only if P1 and P4 passed
#   Stage 3  submit pass A xi=1 and xi=0.5 (782 nodes each)
#   Stage 4  submit pass B x3 (460 nodes each), iff P3 showed Phase 3 engaging
#   Stage 5  submit the xi=1 legacy backfill, afterok on the pass-A xi=1 job
#            (same output file, so it must not run concurrently with it)
#   Stage 6  run pass A xi=2.5 on THIS node, so the campaign starts without
#            waiting in the queue again
# Nothing below here needs a human.

DRYRUN=${PHASE3_PIPELINE_DRYRUN:-0}

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
WORKERS=128                       # physical cores; see relic_phase3.sbatch for why not 256
MAIN_LOG=$LOGS/relic_p3_preflight_${JOBID}.log
P2_LOG=$LOGS/relic_p3_preflight_P2_${JOBID}.log

banner () { printf '\n========== %s ==========\n' "$1"; }
die     () { echo; echo "PIPELINE HALTED: $1" >&2; echo "Nothing further was submitted." >&2
             exit 1; }

# run_slice <label> <driver.py> -- one slice on one node of this allocation.
# Each call blocks until its slice finishes, so callers background it.  The dry-run
# gate is the same hard check relic_phase3.sbatch applies: "starting fresh" means
# Stage 2 did not take, and continuing would drop the finished low-rv columns.
run_slice () {
    local label=$1 driver=$2
    local log=$LOGS/relic_p3_${JOBID}_${label}.out
    local rc=0 dry
    # NB: this must end in `return`, never `exit`.  run_slice is called inside a
    # ( ... && ... ) subshell to chain stage 5 after pass-A xi=1, and an `exit`
    # here would kill that subshell and silently skip the chained slice.
    {
        echo "[$label] $(date '+%F %T') starting, driver=$driver"
        if [ "$DRYRUN" -eq 1 ]; then
            echo "[$label] DRYRUN: would srun -N1 -n1 -c256 python $driver -n $WORKERS"
        else
            dry=$(python "$SCRIPTS/$driver" --dry-run 2>&1)
            echo "$dry"
            if grep -q "starting fresh" <<<"$dry"; then
                echo "[$label] ABORT: driver would start fresh -- Stage 2 did not take."
                rc=1
            else
                srun --nodes=1 --ntasks=1 --cpus-per-task=256 --exclusive --cpu-bind=none \
                     python "$SCRIPTS/$driver" -n "$WORKERS"
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
    python "$SCRIPTS/phase3_make_variants.py" --check \
      || die "driver variants are stale -- rerun phase3_make_variants.py and commit"
    echo "  allocation: $SLURM_JOB_NUM_NODES nodes -- $SLURM_JOB_NODELIST"
    [ "${SLURM_JOB_NUM_NODES:-0}" -ge 6 ] \
      || die "this job holds ${SLURM_JOB_NUM_NODES} nodes but the campaign runs 6 slices
       concurrently.  Resubmit with -N 6."
fi

# ---- Stage 1: preflight ------------------------------------------------
# P2 is serial and may need up to 7200 s on its rtol-1e-8 leg, so it runs
# alongside the pooled probes rather than after them.  It decides KEEP/REVERT on
# phase3_kw's rtol, not GO/NO-GO, so it never blocks: see the policy below.
# This runs on the head node only; the other five sit idle for ~2 h, which costs
# ~10 node-hours out of ~2,200 and buys the campaign a single queue wait.
banner "Stage 1: preflight"
if [ "$DRYRUN" -eq 1 ]; then
    echo "  DRYRUN: skipping preflight; assuming P1/P4 PASS and P3 GO"
    GATE=0
    : >"$MAIN_LOG"
    : >"$P2_LOG"
else
    python "$SCRIPTS/phase3_preflight.py" --only P2 >"$P2_LOG" 2>&1 &
    P2_PID=$!

    python "$SCRIPTS/phase3_preflight.py" --only P1,P3,P4 -n 32 2>&1 | tee "$MAIN_LOG"
    GATE=${PIPESTATUS[0]}

    wait $P2_PID || true
    echo
    echo "--- P2 (ran concurrently) ---"
    cat "$P2_LOG"
fi

[ "$GATE" -eq 0 ] || die "preflight P1/P4 did not pass (exit $GATE).  The anchors or the
       alpha_X axis are wrong -- read the P1/P4 lines in $MAIN_LOG before rerunning."

# P2 policy.  An explicit FAIL means both legs finished and disagreed by >0.5%,
# contradicting validation D -- that is a real result and it stops the campaign.
# An incomplete P2 (the 1e-8 leg exceeding its own 7200 s cap) is the expected
# behaviour of the pathological rv=0.95 node and is evidence FOR loosening the
# Phase-3 rtol, not against it; validation D already covers that ground, so the
# campaign proceeds with a loud note rather than burning another queue cycle.
if grep -q "^P2: incomplete" "$P2_LOG"; then
    echo
    echo "WARNING: P2 did not complete (see $P2_LOG)."
    echo "         Proceeding on validation D; phase3_kw rtol stays 1e-6."
elif grep -qE "^P2:.*FAIL" "$P2_LOG"; then
    die "P2 FAILED: Omega moved by >0.5% between phase3_kw rtol 1e-6 and 1e-8, which
       contradicts validation D.  Set phase3_kw=dict(rtol_value=1e-8) in PHASE3_KW in all
       three pass-A drivers, rerun phase3_make_variants.py, commit, and resubmit."
fi

# P3 gates pass B only.  "DROP" means Phase 3 never engages at eps=1e-10, so those
# slices would only reproduce the legacy numbers at full cost.
RUN_PASS_B=1
grep -q "Stage 4 (pass B.*DROP" "$MAIN_LOG" && RUN_PASS_B=0

# ---- Stage 2: invalidate ----------------------------------------------
banner "Stage 2: invalidate"
if [ "$DRYRUN" -eq 1 ]; then
    echo "  DRYRUN: skipping phase3_invalidate.py"
else
    python "$SCRIPTS/phase3_invalidate.py" || die "phase3_invalidate.py failed"
fi

# ---- Stages 3-5: every slice at once, one per node ---------------------
# Each slice writes its own .npz, so all six are safe to run concurrently and the
# campaign's wall clock is the slowest slice, not the sum.  Stage 5 writes
# relic_grid_vP_eps1e-09.npz, the SAME file as pass-A xi=1, so it is chained after
# that slice on the node it frees rather than run alongside it.
banner "Stages 3-5: launch all slices (one node each)"
echo "  passA xi=2.5  782 nodes"
echo "  passA xi=1.0  782 nodes  -> then stage5 (662 nodes, same output file)"
echo "  passA xi=0.5  782 nodes"
if [ "$RUN_PASS_B" -eq 1 ]; then
    echo "  passB xi=2.5  460 nodes"
    echo "  passB xi=1.0  460 nodes"
    echo "  passB xi=0.5  460 nodes"
else
    echo "  passB SKIPPED: preflight P3 reported DROP (Phase 3 never engaged)."
fi
echo "  per-slice logs: $LOGS/relic_p3_${JOBID}_<label>.out"

declare -A RC

run_slice passA_xi2.5 make_relic_grid_vP_xi_ini_2.5.py & P_A25=$!
run_slice passA_xi0.5 make_relic_grid_vP_xi_ini_0.5.py & P_A05=$!
( run_slice passA_xi1.0 make_relic_grid_vP.py \
    && run_slice stage5 make_relic_grid_vP_stage5.py ) & P_A10=$!

P_B25=""; P_B10=""; P_B05=""
if [ "$RUN_PASS_B" -eq 1 ]; then
    run_slice passB_xi2.5 make_relic_grid_vP_xi_ini_2.5_passB.py & P_B25=$!
    run_slice passB_xi1.0 make_relic_grid_vP_passB.py            & P_B10=$!
    run_slice passB_xi0.5 make_relic_grid_vP_xi_ini_0.5_passB.py & P_B05=$!
fi

wait $P_A25; RC[passA_xi2.5]=$?
wait $P_A05; RC[passA_xi0.5]=$?
wait $P_A10; RC[passA_xi1.0+stage5]=$?
for pid_name in "P_B25 passB_xi2.5" "P_B10 passB_xi1.0" "P_B05 passB_xi0.5"; do
    set -- $pid_name
    pid=${!1}
    [ -n "$pid" ] && { wait "$pid"; RC[$2]=$?; }
done

banner "Campaign finished"
fail=0
for k in "${!RC[@]}"; do
    printf '  %-22s rc=%s\n' "$k" "${RC[$k]}"
    [ "${RC[$k]}" -ne 0 ] && fail=1
done
if [ "$fail" -eq 0 ]; then
    echo
    echo "All slices completed.  Grids updated in GCE/output/relic."
else
    echo
    echo "At least one slice returned non-zero -- check its log before trusting its .npz."
    echo "Slices checkpoint per node, so resubmitting this job resumes where it stopped."
fi
exit $fail
