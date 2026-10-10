"""Reproduce the NaN nodes logged by the 2026-10 Phase-3 campaign, one case per
run, and print the full traceback (docs/phase3_rerun_campaign.md, job 59610830).

The campaign's drivers logged each failed solve as a `NaN node` line but kept
only the exception type, which is not enough to locate the raising expression.
This re-runs a single (mX, mY, alphaX, eps, xi_ini) point with the campaign's
solver keywords and lets the exception propagate.

    python GCE/scripts/phase3_repro_nan.py --case stage5_valueerror
    python GCE/scripts/phase3_repro_nan.py --list
    python GCE/scripts/phase3_repro_nan.py --mX 23.85 --mY 9.707 \
        --alpha 9.382e-05 --eps 1e-9 --xi-ini 1.0 --legacy

`--legacy` runs with PHASE3_KW = {} (Phase 3 off), which is what the stage-5
backfill used; otherwise the campaign's Phase-3 switches are applied.
"""

import argparse
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from hidden_sector_DM.HiddenSectorDM import (          # noqa: E402
    BoltzmannSolver, Cosmology, VectorPortal, rhoc, s0_cosmo,
)

# Campaign solver settings (make_relic_grid_vP*.py as of tag
# phase3-campaign-complete).  Kept here explicitly so this script does not
# depend on the drivers having their campaign constants in place.
PHASE3_KW = dict(rho_switch_thresh=0.3, Sdot_switch_thresh=0.05,
                 enable_phase3=True, phase3_kw=dict(rtol_value=1e-6))
GOH_EXIT_THRESH = 1e99
XMAX_JOINT = 5e3

# One representative per unexplained failure class, taken from the per-slice
# logs.  (mX GeV, mY GeV, alphaX, epsX, xi_ini, legacy_path)
CASES = {
    "stage5_valueerror":   (23.850,  9.707, 9.382e-05, 1e-9, 1.0, True),
    "passA05_valueerror":  (17.022, 12.938, 9.382e-05, 1e-9, 0.5, False),
    "passA05_valueerror2": (100.00, 83.100, 4.401e-04, 1e-9, 0.5, False),
    "passA10_indexerror":  (16.320, 10.850, 4.401e-04, 1e-9, 1.0, False),
    "passA10_indexerror2": (46.821, 26.041, 2.065e-03, 1e-9, 1.0, False),
    "passB05_indexerror":  (21.922, 17.422, 2.132e-01, 1e-10, 0.5, False),
    "passB25_overflow":    (15.000, 13.034, 1.000e+00, 1e-10, 2.5, False),
}


def run(mX, mY, alphaX, eps, xi_ini, legacy):
    """Solve one node, returning (outcome, seconds).  Exceptions are caught and
    their traceback printed -- that is the point of the script."""
    extra = {} if xi_ini == 1.0 else {"xi_ini": xi_ini}
    kw = {} if legacy else dict(PHASE3_KW)

    print(f"mX={mX} GeV  mY={mY} GeV  rv={mY / mX:.4f}  alphaX={alphaX:g}  "
          f"eps={eps:g}  xi_ini={xi_ini}  {'LEGACY (no Phase 3)' if legacy else 'Phase 3 ON'}")
    print("-" * 78)

    solver = BoltzmannSolver(Cosmology(), VectorPortal(
        mX=mX, mY=mY, gX=2, gY=3, alphaX=alphaX,
        include_antiparticlesX=True, include_antiparticlesY=False, **extra))

    t0 = time.perf_counter()
    try:
        res = solver.solve_boltzmann_joint(
            epsX=eps, return_bg_ICs=True, min_eps_floor_ratio=0.0,
            GoH_exit_thresh=GOH_EXIT_THRESH, xmax=XMAX_JOINT,
            phase2_form="Y", verbose=False, **kw)
    except Exception as exc:
        secs = time.perf_counter() - t0
        traceback.print_exc()
        print("-" * 78)
        print(f"RAISED {type(exc).__name__} after {secs:.0f} s")
        return type(exc).__name__, secs

    secs = time.perf_counter() - t0
    pf = res["phase_final"]
    print(f"joint solve OK after {secs:.0f} s: phase_final={pf}")

    # omega_h2_joint returns here only for an Omega-ready exit; every other
    # phase_final makes a SECOND call, and that is where some campaign NaNs
    # were actually raised.
    if pf in ("1.5S", 3):
        omega = res["YX_relic"] * mX * s0_cosmo / rhoc
        print(f"SUCCEEDED (no background step needed): Omega_c h^2 = {omega:.4g}")
        return "ok", secs

    print(f"phase_final={pf} is not Omega-ready -> running solve_background, "
          f"as omega_h2_joint does")
    t1 = time.perf_counter()
    try:
        yield_X = solver.solve_background(
            epsX=eps, bg_ICs=res["bg_ICs"], correct_Y_decay=True,
            traj_has_decay=True, stop_on_rhoY=True, verbose=False)["YX_final"]
    except Exception as exc:
        secs2 = time.perf_counter() - t1
        traceback.print_exc()
        print("-" * 78)
        print(f"RAISED {type(exc).__name__} in solve_background after {secs2:.0f} s")
        print("NOTE: the joint solve itself succeeded; this is the SECOND call.")
        return f"bg:{type(exc).__name__}", secs + secs2

    omega = yield_X * mX * s0_cosmo / rhoc
    print(f"SUCCEEDED after background step: Omega_c h^2 = {omega:.4g}")
    return "ok", secs + (time.perf_counter() - t1)


def main():
    ap = argparse.ArgumentParser(description="reproduce a campaign NaN node")
    ap.add_argument("--case", choices=sorted(CASES))
    ap.add_argument("--list", action="store_true", help="list the known cases")
    ap.add_argument("--mX", type=float)
    ap.add_argument("--mY", type=float)
    ap.add_argument("--alpha", type=float)
    ap.add_argument("--eps", type=float)
    ap.add_argument("--xi-ini", type=float, default=1.0)
    ap.add_argument("--legacy", action="store_true",
                    help="PHASE3_KW = {} (the stage-5 backfill path)")
    args = ap.parse_args()

    if args.list or not (args.case or args.mX):
        for k, v in sorted(CASES.items()):
            print(f"  {k:<22} mX={v[0]:<7} mY={v[1]:<7} alpha={v[2]:<10.4g} "
                  f"eps={v[3]:<6g} xi={v[4]}{'  legacy' if v[5] else ''}")
        return

    if args.case:
        run(*CASES[args.case])
    else:
        run(args.mX, args.mY, args.alpha, args.eps, args.xi_ini, args.legacy)


if __name__ == "__main__":
    main()
