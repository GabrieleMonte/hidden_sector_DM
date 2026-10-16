"""The one SOLVE_FAILED node of the Higgs-portal eps=1e-10 grid: is it real?

relic_grid_hP_eps1e-10.npz records STATUS_SOLVE_FAILED at exactly one of its
2083 nodes, mX=95.872 mY=24.977 (rv=0.2605).  Seven of its eight coupling
points timed out at SOLVE_TIMEOUT=600, leaving one finite sample, and
_crossing needs two.  All eight of its neighbours are STATUS_NO_SOLUTION, so
the expected answer is that this node has no crossing either and the failure
is bookkeeping rather than lost physics.  This settles it by solving the node
again with a generous timeout and timing every point.

Each coupling point is solved in its own worker, so the wall time is the
slowest single point rather than their sum.

    cosmo_env
    sbatch GCE/scripts/probe_failed_node_hP.sbatch

Writes GCE/output/diag/probe_failed_node_hP.npz: omega and seconds per
coupling point, plus the crossing and status they imply.
"""

import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_relic_grid as base                                    # noqa: E402
import make_relic_grid_hP as hp                                   # noqa: E402

# ---- configuration ---------------------------------------------------

NODE_IJ = (44, 37)       # the SOLVE_FAILED node, as stored in the eps=1e-10 grid
PROBE_EPS = 1e-10
PROBE_TIMEOUT = 3600     # s per solve; the campaign ran 600 and before that 7200

DIAG_DIR = base.OUT_DIR.parent / "diag"
OUT = DIAG_DIR / "probe_failed_node_hP.npz"


def one_point(arg):
    """(k, coupling) -> (k, Omega_c h^2, seconds).  NaN Omega on a timeout."""
    k, mX, mY, coupling = arg
    t0 = time.perf_counter()
    omega = base.omega_h2_joint(mX, mY, coupling, PROBE_EPS)
    return k, omega, time.perf_counter() - t0


def main():
    ncores = int(sys.argv[1]) if len(sys.argv) > 1 else base.NCORES

    hp._configure(1.0)
    base._set_coupling(PROBE_EPS)
    base.SOLVE_TIMEOUT = PROBE_TIMEOUT
    base.COUPLING_SWEEP_STOP = False          # we want every point, not a bracket

    i, j = NODE_IJ
    mX, mY = float(base.mX_grid[i, j]), float(base.mY_grid[i, j])
    axis = np.asarray(base.COUPLING)

    print(f"node (i={i}, j={j}): mX={mX:.3f} mY={mY:.3f} rv={mY / mX:.4f}")
    print(f"eps={PROBE_EPS:g}, {axis.size} coupling points, "
          f"timeout {PROBE_TIMEOUT} s, {ncores} workers", flush=True)

    tasks = [(k, mX, mY, float(a)) for k, a in enumerate(axis)]
    omega = np.full(axis.size, np.nan)
    secs = np.full(axis.size, np.nan)

    with mp.Pool(min(ncores, axis.size)) as pool:
        for k, om, dt in pool.imap_unordered(one_point, tasks):
            omega[k], secs[k] = om, dt

    print(f"\n{'coupling':>12} {'Omega_c h^2':>14} {'seconds':>10}  note")
    for k, a in enumerate(axis):
        note = "TIMEOUT" if not np.isfinite(omega[k]) else ""
        if np.isfinite(omega[k]) and omega[k] >= base.OCH2_TARGET:
            note = "above target"
        print(f"{a:12.5g} {omega[k]:14.6g} {secs[k]:10.1f}  {note}")

    crossing, status = base._crossing(omega)
    finite = int(np.isfinite(omega).sum())
    print(f"\nfinite samples : {finite} of {axis.size}")
    print(f"peak Omega     : {np.nanmax(omega):.6g}   target {base.OCH2_TARGET}")
    print(f"crossing       : {crossing!r}")
    print(f"status         : {base.STATUS_NAMES[int(status)]}")
    print(f"slowest solve  : {np.nanmax(secs):.1f} s")

    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, omega=omega, seconds=secs, coupling=axis, mX=mX, mY=mY,
             eps=PROBE_EPS, timeout=PROBE_TIMEOUT, crossing=crossing,
             status=np.int8(status), node_ij=np.array(NODE_IJ))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
