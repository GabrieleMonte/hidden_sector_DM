#!/bin/bash
#SBATCH -J relic-hp
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hp_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hp_%j.err
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
# Higgs-portal relic grids at xi_ini = 1.  Submit this ONE job and log out:
#
#   sbatch GCE/scripts/relic_hp_pipeline.sbatch
#
# Rehearse the control flow without spending anything:
#
#   RELIC_HP_DRYRUN=1 bash GCE/scripts/relic_hp_pipeline.sbatch
#
# ONE QUEUE WAIT, as in portal_pipeline.sbatch: the preflight gates the slices,
# so running it separately would put a second full queue wait (measured 1.6-3.1
# h) on the critical path, and sbatch is refused for -q interactive here.
#
# THREE SLICES, ONE PER NODE: secluded, eps=1e-9, eps=1e-10.  All three write
# different .npz files, so the wall clock is the slowest slice, not the sum.
#
# THE HIGGS PORTAL IS A DIFFERENT HIDDEN SECTOR, not another Gamma_Y.  X is a
# Majorana fermion and Y a real scalar: gY = 1 and include_antiparticlesX =
# False, per the HiggsPortal docstring, set through _HIDDEN_KW in the driver.
# The swept coupling is `lam`, not alphaX -- the preflight prints it under the
# name alpha_relic because it shares relic_alphaX, but it is a Yukawa-like
# coupling and sigmav ~ lam^4, so the crossing moves as the 4th root.
#
# Gamma_phi comes from HDECAY (hidden_sector_DM/phi_decay).  Every solved node
# has mY >= 5 GeV, so only the HDECAY regime is used, never the Winkler or
# pQCD ones.  run_hdecay copies the binary into a private TemporaryDirectory,
# so the 128 workers cannot race; a cold call is ~376 ms and cached ones ~1 us,
# about 6 s of overhead across the whole grid.

DRYRUN=${RELIC_HP_DRYRUN:-0}

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
PRE_LOG=$LOGS/relic_hp_preflight_${JOBID}.log

banner () { printf '\n========== %s ==========\n' "$1"; }
die     () { echo; echo "PIPELINE HALTED: $1" >&2; exit 1; }

# run_slice <only> -- one output, on one node.
run_slice () {
    local only=$1
    local log=$LOGS/relic_hp_${JOBID}_${only}.out
    local rc=0 dry
    # Must end in `return`, never `exit`: callers background these and an exit
    # would kill the caller's subshell.
    {
        echo "[$only] $(date '+%F %T') starting --portal hP --only $only"
        if [ "$DRYRUN" -eq 1 ]; then
            echo "[$only] DRYRUN: would srun make_relic_grid.py --portal hP" \
                 "--xi-ini 1 --only $only -n $WORKERS"
        else
            dry=$(python "$SCRIPTS/make_relic_grid.py" --portal hP \
                  --xi-ini 1 --only "$only" --dry-run 2>&1)
            echo "$dry"
            if grep -q "was built with Phase 3" <<<"$dry"; then
                # "starting fresh" is expected for a brand-new portal; only the
                # Phase-3 mismatch is fatal, and load_state exits non-zero on it.
                echo "[$only] ABORT: Phase-3 setting disagrees with the stored grid."
                rc=1
            fi
            if [ $rc -eq 0 ]; then
                srun --nodes=1 --ntasks=1 --cpus-per-task=256 --exclusive \
                     --cpu-bind=none \
                     python "$SCRIPTS/make_relic_grid.py" --portal hP \
                     --xi-ini 1 --only "$only" -n "$WORKERS"
                rc=$?
                echo "[$only] $(date '+%F %T') finished rc=$rc"
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
    [ "${SLURM_JOB_NUM_NODES:-0}" -ge 3 ] \
      || die "this job holds ${SLURM_JOB_NUM_NODES} nodes but runs 3 slices at once.
       Resubmit with -N 3."
fi

# ---- Stage 1: preflight ------------------------------------------------
# G2 is the gate: at the worst corner the relic crossing must still be
# bracketed below the top of the coupling axis (1.0).  For hP that axis is
# `lam`; lam ~ 1 is already non-perturbative, so a NO-GO means the model is
# off the sampled axis rather than that the axis is merely too short.
banner "Stage 1: hP preflight"
if [ "$DRYRUN" -eq 1 ]; then
    echo "  DRYRUN: skipping preflight; assuming hP is GO"
    : >"$PRE_LOG"
    echo "  hP    production: GO" >>"$PRE_LOG"
else
    python "$SCRIPTS/portal_preflight.py" --portals hP -n 48 2>&1 | tee "$PRE_LOG"
fi

grep -qE "^  hP +production: GO" "$PRE_LOG" \
  || die "hP did not pass the preflight.  Nothing was run (see $PRE_LOG).
       A NO-GO means the relic crossing left the lam axis at the worst corner."

# ---- Stage 2: every slice at once, one per node ------------------------
banner "Stage 2: launch slices"
declare -A RC
PIDS=""
for o in secluded 1e-9 1e-10; do
    run_slice "$o" & PIDS="$PIDS $!:$o"
    echo "  hP_${o}  2083 nodes"
done
echo "  per-slice logs: $LOGS/relic_hp_${JOBID}_<output>.out"

for entry in $PIDS; do
    pid=${entry%%:*}; name=${entry#*:}
    wait "$pid"; RC[$name]=$?
done

banner "Finished"
fail=0
for k in "${!RC[@]}"; do
    printf '  %-10s rc=%s\n' "$k" "${RC[$k]}"
    [ "${RC[$k]}" -ne 0 ] && fail=1
done
if [ "$fail" -eq 0 ]; then
    echo; echo "All slices completed.  New grids in GCE/output/relic."
else
    echo; echo "A slice returned non-zero -- check its log before trusting its .npz."
    echo "Slices checkpoint per node, so resubmitting resumes where it stopped."
fi
exit $fail
