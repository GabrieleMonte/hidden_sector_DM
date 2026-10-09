# Phase-3 rerun campaign — v2 (post cluster review)

v1 (commit 0225dc298f) was reviewed on-cluster 2026-10-09; this version resolves
its five findings and adopts the reviewer's staged plan. Solver: this branch at
`94a43d37b1` or later. Secluded slices are NOT rerun (no in-era injection).

## Status of the hold in 94a43d37b1 (finding 1)

Owner decision 2026-10-09: proceed. The hold asked for a scope map and an
independent review before re-deriving any grid; instead the xi=1/0.5 rerun
doubles as that test. B35 predicts alpha_relic shifts of x1.04 / x1.19 / x1.21
at rv = 0.556 / 0.869 / 0.909 (xi=1, local slope -2). Landing there confirms
the legacy conserved-S bias; landing at x1.00 falsifies it. Same compute.

## Blocks

A. eps = 1e-9, rv >= 0.465 (cols 50-66), xi_ini = 1 / 0.5 / 2.5: 782 nodes/slice.
   Why: B35 bias (+9/+43/+47 % in Omega at xi=1); at xi=2.5 the rv > 0.61
   values are seam garbage outright (x3.55 off at the validated anchor).
B. eps = 1e-10, rv >= 0.636 (cols 57-66), all three slices: 460 nodes/slice.
   Corrected rationale (finding 2): nothing is NaN there and the stored alpha
   axis already reaches 1.0 with ~7x headroom. "Masked" meant the VALIDITY
   ceiling RV_CEILING[1e-10] = 0.61 in notebooks/GCE_fit_vP_xi_test.ipynb —
   the cross-slice ordering test fails above it and the muX-gate audit flagged
   rv >= 0.73. The values exist but are untrusted. Block B is a recompute of
   untrusted values and runs ONLY if preflight P3 shows Phase 3 engages there.

## The alpha_X axis (finding 3 — resolved)

Corrected crossings move UP (x1.04-1.21 at xi=1, x3.55 at the xi=2.5 corner)
and ALPHAX_LIMS[None] = (2e-5, 2e-2) leaves only x1.3-1.9 headroom in the
block; relic_alphaX never extrapolates, so production would write NaN at
exactly the nodes we are paying for. Fix: add `1e-9: (2e-5, 1.0)` to
ALPHAX_LIMS (the axis eps=1e-10 and secluded already use; crossing-
interpolation error is negligible — measured |d2 lnOmega/d lnalpha2| ~ 0.001).
Raising the limit normally resets the file (load_state compares saved alphaX
exactly), so the invalidation snippet below ALSO rewrites the stored axis:
kept low-rv crossings are axis-independent results, and after the swap
load_state resumes with them intact.

## Driver edits (same in all three `make_relic_grid_vP*.py`)

1. `omega_h2_joint`, on the `solve_boltzmann_joint` call, add:
   `rho_switch_thresh=0.3, Sdot_switch_thresh=0.05, enable_phase3=True,
   phase3_kw=dict(rtol_value=1e-6)`
   (finding 5: _phase3_injection_era defaults to rtol 1e-8, two decades
   tighter than SOLVER_KW; validation D showed Omega stable over 1e-6..1e-8.
   Probe P2 re-checks this on one hot node and sets the real per-node cost.)
2. Early return: `if joint["phase_final"] in ("1.5S", 3):` — Phase 3 returns
   an Omega-ready YX_relic; skip `solve_background`.
3. Both `except Exception:` blocks in omega_h2_joint — log before the NaN
   (the S_corr_max=2.0 ValueError is otherwise swallowed silently, and the
   hot near-threshold nodes are exactly where it can trip):
   `except Exception as e:`
   `    print(f"NaN node mX={mX:.3f} mY={mY:.3f} a={alphaX:.2e} eps={eps:g}: "`
   `          f"{type(e).__name__}: {e}", flush=True)`
4. `ALPHAX_LIMS[1e-9] = (2e-5, 1.0)` (see above). Do NOT touch N_RV.
5. `SOLVE_TIMEOUT = 7200` (observed 6200 s at rtol 1e-8; revisit after P2).
   `RUN_SECLUDED = False`.
6. Pass A: `EPS_JOINT = (1e-9,)`,  `ONLY_RV_COLUMNS = tuple(range(50, 67))`.
   Pass B: `EPS_JOINT = (1e-10,)`, `ONLY_RV_COLUMNS = tuple(range(57, 67))`.

## Invalidate stored nodes + swap the 1e-9 axis (run once, after Stage 1)

```python
import numpy as np
NEW_AXIS = np.logspace(np.log10(2e-5), 0.0, 8)      # == new ALPHAX_LIMS[1e-9]
for name, rvmin, swap in [
    ("relic_grid_vP_eps1e-09.npz",            0.46, True),
    ("relic_grid_vP_xi_ini_0.5_eps1e-09.npz", 0.46, True),
    ("relic_grid_vP_xi_ini_2.5_eps1e-09.npz", 0.46, True),
    ("relic_grid_vP_eps1e-10.npz",            0.63, False),
    ("relic_grid_vP_xi_ini_0.5_eps1e-10.npz", 0.63, False),
    ("relic_grid_vP_xi_ini_2.5_eps1e-10.npz", 0.63, False),
]:
    f = "GCE/output/relic/" + name
    d = dict(np.load(f))
    cols = 10.0 ** d["log_rv"] >= rvmin
    d["done_mask"][:, cols] = False
    d["alpha_relic"][:, cols] = np.nan
    d["sigmav"][:, cols] = np.nan
    if swap:
        d["alphaX"] = NEW_AXIS   # kept crossings are axis-independent results
    np.savez(f, **d)
    print(name, int(cols.sum()), "cols cleared", "+ axis -> (2e-5, 1.0)" if swap else "")
```

Old values stay in this branch's git history.

## Staged plan (reviewer's, adopted)

- **Stage 1 — preflight, debug/shared QOS, ~1 h.** Four probes, each gating money:
  P1 anchors, end-to-end through the EDITED driver: mX=34.86, xi=2.5, eps=1e-9,
     rv=0.950 & 0.909 must give alpha_relic = 8.38e-3 / 6.40e-3 to ~1 %.
  P2 rtol: one hot node at phase3_kw rtol 1e-6 vs 1e-8 — Omega within 0.5 %,
     record the speedup (may justify lowering SOLVE_TIMEOUT).
  P3 eps=1e-10 trigger: rv ~ 0.9 node at xi=1 and xi=2.5 — does phase_final
     come back 3 for any alpha? If not, drop Stage 4 entirely (Phase 3 would
     reproduce legacy there and the 0.61 ceiling stays a separate question).
  P4 ceiling: alpha_relic at the worst corner (xi=0.5, high rv, high mX) on the
     new axis — bracketed below 1.0; record the ratio to the stored value.
- **Stage 2** — run the snippet above (only after P1 + P4 pass).
- **Stage 3** — production pass A, one xi at a time, xi=2.5 first (pathology
  proven, anchors in hand). Per slice: --dry-run, check queue = 782, launch.
- **Stage 4** — pass B only if P3 engaged. --dry-run queue = 460 per slice.
- **Stage 5** — backfill the xi=1 holes on the LEGACY path (no Phase-3 kwargs,
  matching the surrounding trusted columns): 662 nodes at rv < 0.465 in
  eps1e-09 (odd 67-mesh columns the 34->67 bisection never filled); optionally
  the empty midpoint columns below col 57 that the embed leaves in eps1e-10.

## File-shape caveats (finding 4)

relic_grid_vP_eps1e-10.npz is 34 columns wide: load_state embeds it into the
67 mesh, so 230 of pass B's 460 xi=1 nodes are new midpoints, not recomputes,
and after the rerun the file is 67 wide with empty low-rv odd columns until
Stage 5. Nothing the notebook used is lost (all 34 old columns land on even
indices) and griddata interpolates across the holes exactly as it did across
the coarser mesh. relic_grid_vP_eps1e-09.npz is itself only 2098/3082 done —
hence Stage 5.

## Resources

~30k solves x 700-1500 s = ~45-100 node-hours at 256 workers/node against
~2,200 left on m3166: budget is fine, a redo is not — hence the stages. No
relic sbatch script exists yet (only chanspec.sh); it needs --requeue, since a
worst-case grid node is ~32 h (8 alphas x 1 retry x 7200 s) and the
CHECKPOINT_EVERY=8 resume makes requeueing safe.
