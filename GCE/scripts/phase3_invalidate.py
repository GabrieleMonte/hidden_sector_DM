"""Phase-3 rerun campaign, Stage 2: one-time state surgery on the six joint
grids (docs/phase3_rerun_campaign.md).  Run from anywhere AFTER preflight
probes P1 + P4 pass:

    python GCE/scripts/phase3_invalidate.py

Per file it clears done_mask / alpha_relic / sigmav in the rerun columns so
load_state re-queues them, and for the eps=1e-9 files it also swaps the
stored alpha_X axis to the new ALPHAX_LIMS[1e-9] = (2e-5, 1.0) one -- the
kept crossings are axis-independent results, so after the swap load_state's
exact-equality axis check passes and the finished low-rv columns survive.

Idempotent (a second run re-clears NaN columns and re-writes the same axis).
Old values stay recoverable in cluster-relic-vP's git history.
"""
from pathlib import Path

import numpy as np

OUT_DIR = Path(__file__).resolve().parents[1] / "output" / "relic"
NEW_AXIS = np.logspace(np.log10(2e-5), 0.0, 8)   # == driver ALPHAX_LIMS[1e-9]

#       file                                     clear rv >=   swap axis
JOBS = [("relic_grid_vP_eps1e-09.npz",            0.46,        True),
        ("relic_grid_vP_xi_ini_0.5_eps1e-09.npz", 0.46,        True),
        ("relic_grid_vP_xi_ini_2.5_eps1e-09.npz", 0.46,        True),
        ("relic_grid_vP_eps1e-10.npz",            0.63,        False),
        ("relic_grid_vP_xi_ini_0.5_eps1e-10.npz", 0.63,        False),
        ("relic_grid_vP_xi_ini_2.5_eps1e-10.npz", 0.63,        False)]


def invalidate(out_dir=OUT_DIR, jobs=JOBS):
    for name, rvmin, swap in jobs:
        path = Path(out_dir) / name
        if not path.exists():
            print(f"{name}: MISSING -- skipped")
            continue
        d = dict(np.load(path))
        cols = 10.0 ** d["log_rv"] >= rvmin
        d["done_mask"][:, cols] = False
        d["coupling_relic"][:, cols] = np.nan
        d["sigmav"][:, cols] = np.nan
        if swap:
            d["coupling_axis"] = NEW_AXIS
        off_grid = ~np.isfinite(d["mX_grid"])
        d["n_done"] = int((d["done_mask"] & ~off_grid).sum())
        np.savez(path, **d)
        print(f"{name}: cleared {int(cols.sum())}/{d['log_rv'].size} rv columns"
              f" (rv >= {rvmin}), {d['n_done']} finished nodes kept"
              + (", alpha axis -> (2e-5, 1.0)" if swap else ""))


if __name__ == "__main__":
    invalidate()
