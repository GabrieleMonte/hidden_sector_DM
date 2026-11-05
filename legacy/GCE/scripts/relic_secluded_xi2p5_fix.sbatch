#!/bin/bash
#SBATCH -J relic-xi2p5-fix
#SBATCH -o legacy/logs/%x_%j.out
#SBATCH -e legacy/logs/%x_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 256
#SBATCH -q regular
#SBATCH -t 12:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#
# Re-solve the vP xi_ini = 2.5 SECLUDED grid with the entropy-correction fix
# (202aeb250d, solver.py: S = (1 + 4I/3)^(3/4) and xi0 reset after the reheating).
#
#   vP, xi_ini 2.5  ->  GCE/output/relic/relic_grid_vP_xi_ini_2.5_eps_secluded.npz
#
# The grid on main is the 2026-11-04 Perlmutter rerun (abda3baf22, trimmed to eps
# 1e-11 and 1e-12), solved with the old S = exp(I).  Re-solving the plotted nodes
# with the fix (a local-slope estimate of each new crossing) grows the eps = 1e-11
# 99% region of the xi_ini figure from 734 to 746 points (top edge, mX 33.5-42.2 GeV,
# rv 0.842-0.856); the other four Z' secluded grids and hP keep their plotted
# regions, so only this one is re-solved
# (legacy/GCE/probes_2026-11-04_secluded_rerun/check_panel_edges.py).
#
# Submit from the root of a clone on main that holds both the fix and the two-eps
# make_relic_grid.py (eps_secluded = (1e-11, 1e-12)); the 2026-11-04 clone does,
# after `git switch main && git pull --ff-only`:
#
#   mkdir -p legacy/logs
#   sbatch legacy/GCE/scripts/relic_secluded_xi2p5_fix.sbatch
#
# The job refuses to start from a tree whose solver or relic script differs from
# HEAD, so the grid is the output of committed code.  The committed grid is moved
# aside first, to legacy/GCE/output/relic_pre_fix/: a run resumes from the file in
# GCE/output/relic/, and the committed one holds a value at every node, so left in
# place it would be reused node for node and nothing would be solved.  Only a file
# still byte-identical to the committed one is moved, so a resubmit after a
# time-out or requeue keeps the partly finished file and resumes from it (every
# output checkpoints per (mX, rv) node).
#
# Cost: one configuration of the 2026-11-04 job, 2083 nodes on 128 workers, ~4.5 h
# at the script's SEC_PER_NODE estimate; a stiff node is capped at 8 couplings x
# 600 s plus one retry.  Anything unfinished at 12 h resumes on resubmission.
#
# Afterwards: commit the new grid on a branch, push, and pull it locally.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd "$SLURM_SUBMIT_DIR"
[ -f GCE/scripts/make_relic_grid.py ] || {
    echo "HALTED: submit from the repository root." >&2
    exit 1; }

DRIVER=GCE/scripts/make_relic_grid.py
RELIC=GCE/output/relic
OLD=legacy/GCE/output/relic_pre_fix
LOGS=legacy/logs
WORKERS=128          # physical cores; 256 logical workers roughly double the
                     # per-solve latency (see relic_portal.sbatch)
PORTAL=vP
XI_INI=2.5
TAG=vP_xi_ini_2.5

# ---- the committed code, with the fix and two secluded eps -----------

git diff --quiet HEAD -- hidden_sector_DM "$DRIVER" || {
    echo "HALTED: hidden_sector_DM/ or $DRIVER differs from HEAD; run committed code." >&2
    exit 1; }
git merge-base --is-ancestor 202aeb250d HEAD 2>/dev/null || {
    echo "HALTED: HEAD does not contain the entropy-correction fix 202aeb250d." >&2
    exit 1; }
grep -q 'eps_secluded=(1e-11, 1e-12),' "$DRIVER" || {
    echo "HALTED: $DRIVER does not carry eps_secluded = (1e-11, 1e-12)." >&2
    exit 1; }
echo "code: $(git rev-parse --short HEAD) ($(git log -1 --format=%cs))"

# ---- move the committed grid aside, once -----------------------------

f=$RELIC/relic_grid_${TAG}_eps_secluded.npz
if [ -f "$f" ] && git ls-files --error-unmatch "$f" >/dev/null 2>&1 \
        && git diff --quiet HEAD -- "$f"; then
    [ -e "$OLD/$(basename "$f")" ] && {
        echo "HALTED: $f is the committed grid, but $OLD already holds one." >&2
        exit 1; }
    mkdir -p "$OLD"
    mv "$f" "$OLD/"
    echo "  $TAG: moved the committed grid to $OLD/"
fi

# ---- dry-run gate: mesh, reuse and queue -----------------------------

DRY=$(python "$DRIVER" --portal "$PORTAL" --xi-ini "$XI_INI" --only secluded --dry-run 2>&1) || {
    echo "$DRY"
    echo "ABORT: dry run refused for $TAG -- see the message above." >&2
    exit 1; }
echo "$DRY"

# ---- solve -----------------------------------------------------------

log=$LOGS/relic_secluded_${TAG}_${SLURM_JOB_ID}.out
{
    echo "[$TAG] $(date '+%F %T') starting"
    srun --nodes=1 --ntasks=1 --cpus-per-task=256 --exclusive --cpu-bind=none \
         python "$DRIVER" --portal "$PORTAL" --xi-ini "$XI_INI" --only secluded -n "$WORKERS"
    rc=$?
    echo "[$TAG] $(date '+%F %T') finished rc=$rc"
} >>"$log" 2>&1

echo "  $TAG rc=$rc (log: $log)"
exit "$rc"
