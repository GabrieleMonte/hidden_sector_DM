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
#SBATCH -t 06:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#
# TEMPORARY (remove once the new grid is in).  Re-solve the top rv rows of the
# vP xi_ini = 2.5 SECLUDED grid with the entropy-correction fix (202aeb250d,
# solver.py: S = (1 + 4I/3)^(3/4) and xi0 reset after the reheating):
#
#   GCE/output/relic/relic_grid_vP_xi_ini_2.5_eps_secluded.npz, rv >= 0.82
#   (rows 0.831, 0.869, 0.909, 0.950; 184 of 2083 nodes, both eps)
#
# The grid on main (2026-11-04 Perlmutter rerun, old S = exp(I)) changes there
# only: Omega at the stored crossing moves by 0.4% (rv 0.831) and 67-127% (rv >=
# 0.869), and by <= 0.1% below (legacy/GCE/probes_2026-11-04_secluded_rerun/).
# Row 0.831 is taken along because it sits under the plotted region's top edge and
# costs no wall time (184 nodes are two passes of 128 workers, like 138).
#
# How: the nodes of those rows are marked not done (coupling_relic and sigmav NaN,
# status unrecorded) in the committed file, and the public driver then resumes it
# as usual, solving exactly those nodes.  The marking runs only while the file is
# still byte-identical to the committed one, so a resubmit after a time-out or
# requeue keeps what was solved and carries on.
#
# Submit from the root of a clone on main (HEAD must contain 202aeb250d):
#
#   sbatch legacy/GCE/scripts/relic_secluded_xi2p5_fix.sbatch
#
# Afterwards: commit the grid on a branch, push, and pull it locally.

module load python
source "$HOME/cosmo_env/bin/activate"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd "$SLURM_SUBMIT_DIR"
[ -f GCE/scripts/make_relic_grid.py ] || {
    echo "HALTED: submit from the repository root." >&2
    exit 1; }
mkdir -p legacy/logs

DRIVER=GCE/scripts/make_relic_grid.py
GRID=GCE/output/relic/relic_grid_vP_xi_ini_2.5_eps_secluded.npz
LOGS=legacy/logs
WORKERS=128          # physical cores; 256 logical workers roughly double the
                     # per-solve latency (see relic_portal.sbatch)
RV_MIN=0.82          # re-solve the rv columns at or above this

# ---- the committed code, with the fix --------------------------------

git diff --quiet HEAD -- hidden_sector_DM "$DRIVER" || {
    echo "HALTED: hidden_sector_DM/ or $DRIVER differs from HEAD; run committed code." >&2
    exit 1; }
git merge-base --is-ancestor 202aeb250d HEAD 2>/dev/null || {
    echo "HALTED: HEAD does not contain the entropy-correction fix 202aeb250d." >&2
    exit 1; }
echo "code: $(git rev-parse --short HEAD) ($(git log -1 --format=%cs))"

# ---- mark the top rv rows not done, once -----------------------------

if git diff --quiet HEAD -- "$GRID"; then
    python - "$GRID" "$RV_MIN" <<'EOF' || { echo "ABORT: marking failed." >&2; exit 1; }
import sys

import numpy as np

path, rv_min = sys.argv[1], float(sys.argv[2])
with np.load(path) as saved:
    grid = {k: saved[k] for k in saved.files}

valid = np.isfinite(grid['mX_grid'])                    # off-grid nodes stay done
redo = valid & (10**grid['log_rv'] >= rv_min)[None, :]
grid['done_mask'][redo] = False
grid['coupling_relic'][:, redo] = np.nan
grid['sigmav'][:, redo] = np.nan
grid['status'][:, redo] = -1                            # STATUS_UNRECORDED
grid['n_done'] = np.array(int((grid['done_mask'] & valid).sum()))
np.savez(path, **grid)
print(f'  marked {int(redo.sum())} nodes at rv >= {rv_min} not done; '
      f'{int(grid["n_done"])} of {int(valid.sum())} kept')
EOF
else
    echo "  $GRID already differs from HEAD: resuming it"
fi

# ---- dry-run gate: reuse and queue -----------------------------------

DRY=$(python "$DRIVER" --portal vP --xi-ini 2.5 --only secluded --dry-run 2>&1) || {
    echo "$DRY"
    echo "ABORT: dry run refused -- see the message above." >&2
    exit 1; }
echo "$DRY"

# ---- solve -----------------------------------------------------------

log=$LOGS/relic_secluded_vP_xi_ini_2.5_${SLURM_JOB_ID}.out
{
    echo "[vP_xi_ini_2.5] $(date '+%F %T') starting"
    srun --nodes=1 --ntasks=1 --cpus-per-task=256 --exclusive --cpu-bind=none \
         python "$DRIVER" --portal vP --xi-ini 2.5 --only secluded -n "$WORKERS"
    rc=$?
    echo "[vP_xi_ini_2.5] $(date '+%F %T') finished rc=$rc"
} >>"$log" 2>&1

echo "  rc=$rc (log: $log)"
exit "$rc"
