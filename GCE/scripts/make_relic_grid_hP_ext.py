"""lam that reproduces Omega_c h^2 = 0.12 for the Higgs-portal cascade over the
mX = 100-130 GeV extension, the region the main (15, 100] GeV grid cannot close.

    cosmo_env
    python GCE/scripts/make_relic_grid_hP_ext.py --dry-run   # grid and ETA
    python GCE/scripts/make_relic_grid_hP_ext.py -n 64

Four solver files land in GCE/output/relic/, beside the main grid and never into
it, so neither driver's checkpoint can clobber the other's:

    relic_grid_hP_ext_eps1e-07.npz    relic_grid_hP_ext_eps1e-09.npz
    relic_grid_hP_ext_eps1e-08.npz    relic_grid_hP_ext_eps1e-10.npz

Those are checkpoints, not the grid anyone reads.  `_merge_into_main` then pads each
one onto the main grid's mesh and writes the result to GCE/output/sv_relic_fits/ as
relic_grid_hP_eps{eps}.npz, 52 x 67 where the main grid is 46 x 67.  That is the
single file per eps the notebooks load, so nothing downstream concatenates anything.

The model, the lam axis, the solver settings, the checkpointing and the
crossing extraction are make_relic_grid_hP.py's and make_relic_grid.py's, used
unchanged.  This script supplies only the mesh.

WHY
---
The hP favoured region (south, free <sigma v>) still reaches the mX = 100 GeV
edge of the GCE grid along the mY ~ mX diagonal, so both the 2- and 3-sigma
contours have to be closed from above; docs/hP_extension_campaign.md is the
campaign.  Stage 2 builds the photon spectra there; this is its relic half.
The Higgs portal only: the vector portals' contours close inside the main grid.

GRID
----
The extension is its own file on its own mesh rather than a widened main grid,
so the main grid's 2083 solved nodes are neither re-read nor rewritten.  The
mesh continues the main one's spacing exactly, so the two concatenate with no
interpolation across the seam at mX = 100 GeV:

  columns   log_mX continued at the main step du = log10(100/15)/45 = 0.0183091
            dex for N_MX_EXT columns -> mX = 104.31 ... 128.78 GeV.  A seventh
            column would land at 134.3, past the 130 GeV the campaign asks for.
  rows      the main log_rv axis restricted to rv >= RV_MIN_EXT, i.e. 11 of its
            67 rows, rv = 0.6081 ... 0.9500 at the same dv = 0.0193751 dex.

Why the rows stop at rv = 0.6: Stage 2 generates the spectra for this region
over rv 0.6-1 only, because for mX >= 92 the 2- and 3-sigma GCE regions sit at
rv >= 0.70, narrowing to >= 0.82 at mX = 100.  A relic solved below that could
never be paired with a spectrum.

mY is deliberately NOT capped at the main grid's 100 GeV ceiling -- 15 of the
66 nodes sit above it, up to mY = 122.34 GeV -- which is Stage 2's choice too.

66 nodes x 14 lam samples per eps, so ~1 node-hour per output at most.

EPSILONS
--------
1e-7 and 1e-8 go through the joint solver, one file each.  1e-9 and 1e-10 go
through the secluded route, as on the main grid and for the reason set out in
make_relic_grid_hP.EPS_SECLUDED: at that coupling the sectors have already
decoupled and the joint route only costs more.  run_secluded writes them as one
file with eps on the leading axis, so _split_secluded then writes the per-eps
copies downstream reads -- the main grid's relic_grid_hP_eps1e-09.npz and
_eps1e-10.npz are exactly these slices (checked 2026-10-27: identical to slices
0 and 1 of relic_grid_hP_eps_secluded.npz).
"""

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_relic_grid as base                                   # noqa: E402
import make_relic_grid_hP as hP                                  # noqa: E402

# ---- configuration ---------------------------------------------------

N_MX_EXT = 6             # continuation columns past the main grid's mX = 100 GeV
RV_MIN_EXT = 0.6         # lowest mY/mX kept, matching the Stage-2 spectra; see GRID

EPS_JOINT = (1e-7, 1e-8)             # joint route, one output file each
EPS_SECLUDED = hP.EPS_SECLUDED       # (1e-9, 1e-10), secluded route, one file then split

SV_RELIC_DIR = base.OUT_DIR.parent / "sv_relic_fits"       # what the notebooks read
EMBED_TOL = 1e-9             # the tolerance make_relic_grid._embed uses on a log axis

# How a node outside the extension's rv range is marked, matching make_relic_grid:
# off the domain in mass, already "done" so no sweep picks it up, status unrecorded.
PAD_VALUE = {"mX_grid": np.nan, "mY_grid": np.nan, "coupling_relic": np.nan,
             "sigmav": np.nan, "done_mask": True, "status": base.STATUS_UNRECORDED}

TAG_STEM = "hP_ext"      # output-name stem; xi_ini != 1 is appended


def _set_mesh():
    """Replace make_relic_grid's module-level mesh with the extension's.

    Must run before load_state and before any worker is forked: _pending,
    load_state and _save all read these globals, and the pool inherits them.
    """
    du = base.log_mX[1] - base.log_mX[0]
    log_mX = base.log_mX[-1] + du*np.arange(1, N_MX_EXT + 1)
    log_rv = base.log_rv[base.log_rv >= np.log10(RV_MIN_EXT)]

    mX_grid, rv_grid = np.meshgrid(10.0**log_mX, 10.0**log_rv, indexing="ij")
    mY_grid = rv_grid*mX_grid
    off_grid = mY_grid < base.MY_LIM[0]      # empty at RV_MIN_EXT = 0.6; kept so the
    mX_grid[off_grid] = np.nan               # mesh still obeys the mY floor if that
    mY_grid[off_grid] = np.nan               # configuration is ever lowered

    base.log_mX, base.log_rv = log_mX, log_rv
    base.mX_grid, base.mY_grid, base.off_grid = mX_grid, mY_grid, off_grid
    base.N_MX, base.N_RV = log_mX.size, log_rv.size
    base.N_VALID = int(np.isfinite(mX_grid).sum())
    base.MX_LIM = (10.0**log_mX[0], 10.0**log_mX[-1])
    base.RV_MAX = float(10.0**log_rv[-1])


def _split_secluded():
    """Write the per-eps copies of the secluded output.

    run_secluded writes one file with eps on the leading axis; the notebooks
    read one file per eps.  Pure reshaping -- nothing is re-solved.
    """
    path = base.OUT_DIR / f"relic_grid_{base.TAG}_eps_secluded.npz"
    if not path.exists():
        raise RuntimeError(f"{path.name} not written; nothing to split")

    base._set_coupling("secluded")           # _save records COUPLING as the axis
    with np.load(path, allow_pickle=True) as d:
        for k, eps in enumerate(d["eps"]):
            out = base.OUT_DIR / f"relic_grid_{base.TAG}_eps{eps:.0e}.npz"
            base._save(out, float(eps), d["coupling_relic"][k], d["sigmav"][k],
                       d["done_mask"], d["status"][k])
            n_ok = int(np.isfinite(d["coupling_relic"][k]).sum())
            print(f"  {out.name}: {n_ok}/{base.N_VALID} nodes with a crossing")


def _merge_into_main():
    """Pad each extension file onto the main grid's mesh and write the merged grid.

    The two meshes share their log-mX step and the extension's rows ARE main rows,
    so this is padding, not interpolation: the extension's six columns are placed at
    the rv rows they were solved on and the rest of those columns are marked
    off-domain the way make_relic_grid marks one -- mX_grid and mY_grid NaN,
    done_mask True, status unrecorded.  Nothing is re-solved.

    Writes into SV_RELIC_DIR rather than over the solver's own file, so re-running
    either driver still resumes from an untouched checkpoint.
    """
    for eps in EPS_JOINT + EPS_SECLUDED:
        main_path = base.OUT_DIR / f"relic_grid_hP_eps{eps:.0e}.npz"
        ext_path = base.OUT_DIR / f"relic_grid_{base.TAG}_eps{eps:.0e}.npz"
        if not (main_path.exists() and ext_path.exists()):
            raise RuntimeError(f"need both {main_path.name} and {ext_path.name}")

        with np.load(main_path) as m, np.load(ext_path) as e:
            for key in ("coupling_axis", "eps"):
                if not np.array_equal(m[key], e[key]):
                    raise RuntimeError(f"{key} differs between the two {eps:.0e} files")
            log_mX, log_rv = m["log_mX"], m["log_rv"]
            if e["log_mX"].min() <= log_mX.max():
                raise RuntimeError("the extension overlaps the main mesh in mX")

            rows = [int(np.argmin(np.abs(log_rv - v))) for v in e["log_rv"]]
            if max(abs(log_rv[r] - v) for r, v in zip(rows, e["log_rv"])) > EMBED_TOL:
                raise RuntimeError("the extension's rv rows are not main grid rows")

            n_mx = log_mX.size + e["log_mX"].size
            out = {"log_mX": np.concatenate([log_mX, e["log_mX"]]), "log_rv": log_rv,
                   "coupling_axis": m["coupling_axis"], "eps": m["eps"],
                   "xi_ini": m["xi_ini"], "phase3": m["phase3"], "portal": m["portal"],
                   "n_done": np.int64(int(m["n_done"]) + int(e["n_done"]))}
            for key, pad in PAD_VALUE.items():
                a = np.full((n_mx, log_rv.size), pad, dtype=m[key].dtype)
                a[:log_mX.size] = m[key]
                a[log_mX.size:, rows] = e[key]
                out[key] = a

        SV_RELIC_DIR.mkdir(parents=True, exist_ok=True)
        np.savez(SV_RELIC_DIR / main_path.name, **out)
        n_domain = int(np.isfinite(out["mX_grid"]).sum())
        print(f"  {main_path.name}: {n_mx}x{log_rv.size}, "
              f"{int(out['n_done'])}/{n_domain} nodes done")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--xi-ini", type=float, default=1.0,
                    help="T_D/T_SM at Tinf (default 1.0)")
    ap.add_argument("--only", default=None, choices=["secluded", "1e-7", "1e-8"],
                    help="run just one output.  Each writes a different file, "
                         "so they can be separate concurrent slices.  "
                         "'secluded' is 1e-9 and 1e-10 together, as on the "
                         "main grid")
    ap.add_argument("-n", "--ncores", type=int, default=base.NCORES,
                    help=f"worker processes (default {base.NCORES})")
    # Run control, not configuration: like -n these say how to run, not what
    # to compute.
    ap.add_argument("--retry-failed", action="store_true",
                    help="also re-solve stored nodes whose status is "
                         "solve_failed (pair with --solve-timeout)")
    ap.add_argument("--solve-timeout", type=int, default=None, metavar="S",
                    help=f"seconds per solve (default {base.SOLVE_TIMEOUT})")
    ap.add_argument("--dry-run", action="store_true",
                    help="report the grid, what is reusable and the ETA, then exit")
    args = ap.parse_args()

    hP._configure(args.xi_ini)
    base.TAG = TAG_STEM if base.XI_INI == 1.0 else f"{TAG_STEM}_xi_ini_{base.XI_INI:g}"
    base.EPS_JOINT = EPS_JOINT
    base.EPS_SECLUDED = EPS_SECLUDED
    base.RUN_SECLUDED = True
    _set_mesh()

    if args.retry_failed:
        base.RETRY_STATUSES = (base.STATUS_SOLVE_FAILED,)
    if args.solve_timeout is not None:
        base.SOLVE_TIMEOUT = args.solve_timeout
    if args.only is not None:
        base.RUN_SECLUDED = args.only == "secluded"
        base.EPS_JOINT = () if base.RUN_SECLUDED else (float(args.only),)

    base.OUT_DIR.mkdir(parents=True, exist_ok=True)

    eps_run = ("secluded " + ", ".join(f"{e:g}" for e in base.EPS_SECLUDED)
               if base.RUN_SECLUDED else "")
    eps_run = " + ".join(x for x in [eps_run, ", ".join(f"joint {e:g}" for e in
                                                        base.EPS_JOINT)] if x)
    print(f"portal: HiggsPortal  xi_ini={base.XI_INI:g}  "
          f"gX={hP.GX_MAJORANA} gY={hP.GY_REAL_SCALAR} Majorana X  |  "
          f"Phase 3 {'ON' if base.PHASE3_KW else 'OFF'}  |  tag '{base.TAG}'")
    print(f"grid: {base.N_MX} x {base.N_RV} = {base.N_VALID} nodes, "
          f"mX in [{base.MX_LIM[0]:.2f}, {base.MX_LIM[1]:.2f}] GeV, "
          f"rv in [{10.0**base.log_rv[0]:.4f}, {base.RV_MAX:.4f}], "
          f"mY in [{np.nanmin(base.mY_grid):.2f}, {np.nanmax(base.mY_grid):.2f}] GeV")
    print(f"outputs: {eps_run}  |  lam axis {hP.LAM_LIMS[None]} "
          f"x {hP.N_LAM[None]} points  |  out -> {base.OUT_DIR}")

    if args.dry_run:
        base.dry_run()
        return

    print(f"workers: {args.ncores}")
    if base.RUN_SECLUDED:
        base.run_secluded(args.ncores)
        _split_secluded()
    for eps in base.EPS_JOINT:
        base.run_joint(eps, args.ncores)

    # Only a run that produced every eps can write the merged grids; a --only slice
    # leaves the others untouched and merges once the last slice has landed.
    if args.only is None:
        print("merging onto the main grid's mesh:")
        _merge_into_main()


if __name__ == "__main__":
    main()
