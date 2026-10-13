"""Phase-3 follow-up, ROUND 2 (docs/phase3_followup.md): one-time state surgery
on the two xi=1 joint grids it re-solves.  Run AFTER phase3_followup_preflight.py
passes:

    python GCE/scripts/phase3_followup_invalidate.py

Per file it clears done_mask / coupling_relic / sigmav (and phase_path /
p3_sdot, where a file already has them) in the follow-up's rv columns, so the
variant drivers re-queue them.  The stored alpha_X axis and eps
are checked, never rewritten: both are already what the drivers expect, and a
mismatch would make load_state start fresh and drop every other column.

RUN ONCE.  Unlike the campaign's phase3_invalidate.py this is NOT idempotent in
effect: a second run would also clear the nodes the job has re-solved since.
So it writes MARKER after the last file and refuses to run while MARKER exists
-- a requeued phase3_followup.sbatch therefore resumes from the drivers'
checkpoints.  --force overrides that (only for a deliberate fresh start).
Round 1 left its own marker (phase3_followup_invalidated.txt); round 2 has a
separate one, so the two never block each other.  Old values stay recoverable
from this branch's git history.
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

OUT_DIR = Path(__file__).resolve().parents[1] / "output" / "relic"
MARKER_NAME = "phase3_followup2_invalidated.txt"
AXIS = np.logspace(np.log10(2e-5), 0.0, 8)   # == driver ALPHAX_LIMS[1e-9] and [1e-10]

#       file                                    eps    clear cols [lo, hi)   nodes
JOBS = [("relic_grid_vP_eps1e-09.npz", 1e-9,  (0, 67),  2083),
        ("relic_grid_vP_eps1e-10.npz", 1e-10, (57, 67), 460)]


def _check(d, name, eps):
    """Refuse a file the variant driver would not resume from."""
    problems = []
    if d["coupling_relic"].shape != (46, 67):
        problems.append(f"shape {d['coupling_relic'].shape}, expected (46, 67)")
    if not np.array_equal(d["coupling_axis"], AXIS):
        problems.append("stored alpha_X axis is not logspace(2e-5, 1.0, 8)")
    if not np.array_equal(d["eps"], np.atleast_1d(eps)):
        problems.append(f"stored eps {d['eps']} != {eps:g}")
    if problems:
        raise SystemExit(f"{name}: " + "; ".join(problems) + " -- nothing written")


def invalidate(out_dir=OUT_DIR, jobs=JOBS, force=False):
    out_dir = Path(out_dir)
    marker = out_dir / MARKER_NAME
    if marker.exists() and not force:
        print(f"{marker.name} exists -- already invalidated, nothing done:")
        print("  " + marker.read_text().strip().replace("\n", "\n  "))
        return

    loaded = []
    for name, eps, cols, nodes in jobs:          # check every file before touching any
        path = out_dir / name
        if not path.exists():
            raise SystemExit(f"{name}: MISSING -- nothing written")
        d = dict(np.load(path))
        _check(d, name, eps)
        loaded.append((path, d, cols, nodes))

    lines = []
    for path, d, (lo, hi), nodes in loaded:
        off_grid = ~np.isfinite(d["mX_grid"])
        sel = np.zeros(d["done_mask"].shape, bool)
        sel[:, lo:hi] = True
        sel &= ~off_grid
        if int(sel.sum()) != nodes:
            raise SystemExit(f"{path.name}: cols {lo}-{hi - 1} hold {int(sel.sum())} "
                             f"valid nodes, expected {nodes} -- stopping")
        d["done_mask"][sel] = False
        d["coupling_relic"][sel] = np.nan
        d["sigmav"][sel] = np.nan
        for key in ("phase_path", "p3_sdot"):
            if key in d:
                d[key][sel] = np.nan
        d["n_done"] = int((d["done_mask"] & ~off_grid).sum())
        tmp = path.with_name(path.stem + ".tmp.npz")
        np.savez(tmp, **d)
        tmp.replace(path)
        msg = (f"{path.name}: cleared cols {lo}-{hi - 1} ({nodes} nodes), "
               f"{d['n_done']} finished nodes kept")
        print(msg)
        lines.append(msg)

    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=out_dir,
                              capture_output=True, text=True).stdout.strip() or "unknown"
    except OSError:
        head = "unknown"
    marker.write_text(f"invalidated {time.strftime('%Y-%m-%d %H:%M:%S')} at commit {head}\n"
                      + "\n".join(lines) + "\n")
    print(f"wrote {marker.name}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Phase-3 follow-up round-2 state surgery (run once)")
    ap.add_argument("--force", action="store_true",
                    help="run even though the marker says it already ran "
                         "(discards any re-solved follow-up nodes)")
    invalidate(force=ap.parse_args().force)
