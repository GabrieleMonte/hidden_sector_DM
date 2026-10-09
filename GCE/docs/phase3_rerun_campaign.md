# Phase-3 rerun campaign (Oct 2026)

Recompute every joint-grid block whose stored `alpha_relic` is untrusted, using the
Phase-3 (entropy-injection era) solver — cluster-relic-vP commit `94a43d37b1` or later.
Secluded slices are NOT rerun (no in-era injection; unaffected).

## What and why

| Block (all three xi_ini: 1, 0.5, 2.5) | rv columns | Why |
|---|---|---|
| eps = 1e-9, rv >= 0.465 (cols 50–66, 17 cols) | 782 nodes/slice | Legacy conserved-S bias under injection: +9/+43/+47 % in Omega at rv = 0.556/0.869/0.909 (xi=1); at xi=2.5 the old high-rv values are outright seam garbage |
| eps = 1e-10, rv >= 0.636 (cols 57–66, 10 cols) | 460 nodes/slice | Currently masked (notebook ceiling 0.61); Phase 3 should recover it |

Total ~3.7k nodes x 8 alpha-solves; Phase-3 solves run 200–1500 s (tail ~6000 s).

**Preflight (cheap, do first):** Phase 3 has only been validated at eps=1e-9.
Before burning the eps=1e-10 block, run 2 nodes at eps=1e-10 (one rv~0.5 where the
legacy answer is trusted — must agree to ~1% — one rv~0.9 sanity).

## Driver edits (same in all three `make_relic_grid_vP*.py`)

1. In `omega_h2_joint`, add to the `solve_boltzmann_joint` call:
   `rho_switch_thresh=0.3, Sdot_switch_thresh=0.05, enable_phase3=True`
   and change the early return to `if joint["phase_final"] in ("1.5S", 3):`
   (Phase 3 returns an Omega-ready `YX_relic`; skip `solve_background`).
2. `SOLVE_TIMEOUT = 7200` (observed 6200 s Phase-3 node).
3. `RUN_SECLUDED = False`. Do NOT touch `ALPHAX_LIMS` or `N_RV` (load_state resets).
4. Per pass: `EPS_JOINT = (1e-9,)`  with `ONLY_RV_COLUMNS = tuple(range(50, 67))`,
   then      `EPS_JOINT = (1e-10,)` with `ONLY_RV_COLUMNS = tuple(range(57, 67))`.

## Invalidate the stored nodes (resume mask skips done nodes otherwise)

Run once per output file before launching (works on both 34- and 67-wide files;
old values stay in git history on this branch):

```python
import numpy as np
for f, rvmin in [
    ("GCE/output/relic/relic_grid_vP_eps1e-09.npz",            0.46),
    ("GCE/output/relic/relic_grid_vP_xi_ini_0.5_eps1e-09.npz", 0.46),
    ("GCE/output/relic/relic_grid_vP_xi_ini_2.5_eps1e-09.npz", 0.46),
    ("GCE/output/relic/relic_grid_vP_eps1e-10.npz",            0.63),
    ("GCE/output/relic/relic_grid_vP_xi_ini_0.5_eps1e-10.npz", 0.63),
    ("GCE/output/relic/relic_grid_vP_xi_ini_2.5_eps1e-10.npz", 0.63),
]:
    d = dict(np.load(f))
    cols = 10.0 ** d["log_rv"] >= rvmin
    d["done_mask"][:, cols] = False
    d["alpha_relic"][:, cols] = np.nan
    d["sigmav"][:, cols] = np.nan
    np.savez(f, **d)
    print(f, int(cols.sum()), "columns cleared")
```

Then `--dry-run` each script and check the queue: 782 nodes per eps=1e-9 slice,
460 per eps=1e-10 slice.

## Sanity anchors (rerun must reproduce to ~1 %)

mX = 34.86, xi_ini = 2.5, eps = 1e-9: alpha_relic = 8.38e-3 at rv = 0.950,
6.40e-3 at rv = 0.909 (Phase-3 validated, local slope d ln Omega / d ln alpha = -2.0).
