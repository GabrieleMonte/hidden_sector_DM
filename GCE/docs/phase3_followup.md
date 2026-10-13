# Phase-3 follow-up

## Round 2 (2026-10-13): the xi=1 slices, every alpha on Phase 3

### What round 1 got wrong

Round 1 (job 59725646, commit `9098886f73`, described below) forced the Sdot
trigger to 1e-3. xi=2.5 came back clean: 0 isolated nodes, and the mX=34.86
ratios on prediction. xi=1 did not:
- isolated nodes went 18 -> 44 (eps=1e-9) and 2 -> 14 (eps=1e-10, cols 57-66);
- most of those are neighbours of about 20 genuinely bad nodes, every one moved
  DOWN, by 0.18-0.83x in alpha.

Cause (full alpha axes of 7 bad nodes re-solved locally, scratch
`fu_glitch_axis.py`): at one or more alphas per node, often the one nearest
the crossing, the 1e-3 trigger handed off early without running Phase 3. The
legacy path then gave Omega 2-7x too LOW (e.g. 0.018 vs 0.121), so that alpha
looked below target and the crossing fell one step. This is the 1e-7 failure
mode, only rarer. Lowering the trigger is not automatically safer.

Re-solving those alphas at Sdot 1e-4 engaged Phase 3 in all 13 misses. The
repaired crossings match the pre-follow-up (campaign) values to 0.1-1.6 %.

### What round 2 changes

Three driver switches, all off by default (the source drivers behave exactly
as before), all on in the generated `*_fu2_*` variants:
- `PHASE3_REQUIRE`: an alpha that ends on the legacy path is re-solved at Sdot
  1e-4, then 1e-5. If none engages, the legacy value is kept, logged as
  `LEGACY alpha ...` and flagged.
- `ALPHAX_SWEEP_STOP`: alphas are solved upward and the sweep stops one point
  past the first Omega below 0.12. The crossing is read off the four points
  around it (`relic_alphaX_local`), not a spline through all eight. The skipped
  alphas are the slowest solves (773 + 849 s of a 2162 s node) and the least
  reliable: alpha=1 has returned Omega = 1.7e-144, which pulled the global
  spline by ~1 %. On complete axes the local crossing matches the global one
  to <= 0.2 % (14 nodes).
- Every alpha's path is saved per node: `phase_path` (3 Phase 3, 4 SM-locked,
  1 legacy, 0 failed, -1 skipped) and `p3_sdot` (the trigger that engaged).

The retry runs at every alpha, not only near the crossing: at mX=24.88
rv=0.509 three alphas in a row missed (9.4e-5 .. 2.1e-3). With the stop rule
alone, a legacy miss would look like "crossing found" and lock in the wrong
answer.

### What is re-solved

| slice | rv columns | nodes |
|---|---|---|
| xi=1, eps=1e-9 | 0-66 (all) | 2083 |
| xi=1, eps=1e-10 | 57-66 | 460 |

The isolated-node test only catches shifts > 20 %, so a smaller miss could sit
anywhere in these slices; hence the whole slices, not just the flagged nodes.
xi=2.5 and xi=0.5 are not touched.

### Run it

```bash
sbatch GCE/scripts/phase3_followup.sbatch      # one job, 2 nodes, requeue-safe
```

Preflight G1 runs `node_joint` of the generated drivers on three round-1 glitch
nodes: the crossing must come back to the repaired value (0.5 %), with no
legacy or failed alpha, at least one alpha via the retry, and the tail skipped.
Then the invalidation (marker `phase3_followup2_invalidated.txt`, separate from
round 1's), the two slices, and the report. Rehearse with
`PHASE3_FOLLOWUP_DRYRUN=1 bash GCE/scripts/phase3_followup.sbatch`. When it
ends, commit `GCE/output/relic` (including the round-2 marker) on this branch.

Expect from the report:
- 0 queued, 0 NaN, 0 (or a listed handful of) legacy/failed alphas;
- isolated nodes back to about the pre-follow-up level: eps=1e-9 had 18 before
  round 1, all from the campaign's own mixed paths, so near 0-1 now; eps=1e-10
  cols 57-66 near 2 (two rv >= 0.87 nodes at mX=95.87 were isolated before any
  follow-up);
- alpha_new/alpha_campaign ~ 1 over the campaign's columns, except the
  campaign's own glitch nodes.

Known limit, not addressed: at xi=1 eps=1e-10 mX=34.86 rv=0.95 the direct
solve at the stored alpha gives Omega/0.12 = 1.12, but the 8-point axis puts the
crossing where interpolation says 1.00: Omega alpha^2 has a ~20 % bump between
the sampled alphas 9.7e-3 and 4.5e-2. rv=0.95 is above the notebook's eps=1e-10
ceiling (0.9) and is masked there.

The round-1 variants (`*_fu_e09.py`, `*_fu_e10.py`) stay in the tree as the
record of what round 1 ran; the generator now writes only the round-2 ones.

---

## Round 1 (2026-10-12)

The 2026-10 campaign (`phase3_rerun_campaign.md`) re-solved only rv >= 0.465
(eps=1e-9) and rv >= 0.636 (eps=1e-10). Below that cut the grids still hold
Phase-3-OFF (legacy) crossings, and at xi_ini = 2.5 those are wrong by up to
40 % in <sigma v>. This job re-solves them, plus the xi=1 nodes whose alpha
axis straddled the trigger.

### Run it

DONE (job 59725646). The `phase3_followup*` scripts now run round 2; round 1's
versions are at commit `9098886f73`.

```bash
sbatch GCE/scripts/phase3_followup.sbatch      # one job, 3 nodes, requeue-safe
```

The job runs the preflight, then the invalidation, then four slices, then a report.
Rehearse first with `PHASE3_FOLLOWUP_DRYRUN=1 bash GCE/scripts/phase3_followup.sbatch`.
When it ends, commit `GCE/output/relic` (including
`phase3_followup_invalidated.txt`) on this branch.

### What is re-solved (Phase 3 forced: Sdot_switch_thresh 1e-3)

| slice | rv columns | nodes | why |
|---|---|---|---|
| xi=2.5, eps=1e-9 | 0-49 (rv < 0.465) | 1301 | legacy below the cut: -17 % .. +13 % |
| xi=2.5, eps=1e-10 | 0-56 (rv < 0.636) | 1623 | legacy below the cut: +23 % .. +40 % |
| xi=1, eps=1e-9 | 0-66 (all) | 2083 | below the cut: -8 % .. +4 %; above it, 18 glitch nodes (e.g. rv=0.760, 0.7x their neighbours) |
| xi=1, eps=1e-10 | 57-66 (rv >= 0.636) | 460 | mixed path at rv=0.95: +12 % |

The percentages are the change in <sigma v>_relic: Omega_forced/0.12 at the
stored crossing, measured at mX=34.86 (rv 0.159-0.950). They are not
re-solved: xi=0.5 (at most 4 %, inside the trigger's design tolerance),
xi=1 eps=1e-10 below col 57 (under 1 %), and every secluded grid.

### Why the trigger is fine and forcing it is still right

Phase 3 fires when dlnS/dN >= 0.05. The legacy error grows with the entropy
injected by Y decays, which is ~xi^3, so at xi=2.5 the trigger fires at every
rv. Re-solving the stored crossings with the UNCHANGED campaign trigger engages
Phase 3 and agrees with the forced trigger to < 0.2 %. The large errors exist
only because the campaign's cut left those columns on the legacy path.

The trigger only does damage where it fires for some alphas of a node and not
others. The crossing then splines across two different solvers: those are the
xi=1 glitch nodes, plus the eps=1e-10 rv=0.95 column. Forcing (Sdot 1e-3) puts
every alpha on one path, and it reproduces the campaign wherever the campaign
engaged. That is why it is used for all four slices.

Full-axis check (8 alphas, both triggers, 2026-10-12):
- At xi=2.5 (eps=1e-9 rv=0.326; eps=1e-10 rv=0.05 and 0.426) the campaign
  trigger engages at all 8 alphas and agrees with forced to <= 0.2 %.
- At xi=1, eps=1e-9, rv=0.326 it leaves 2 of 8 alphas on the legacy path,
  including alpha=2.06e-3 next to the crossing: <sigma v> off by 5.5 %. So
  below the cut at xi=1 the campaign setting would itself make mixed-path
  nodes.
- At the glitch node mX=24.9, rv=0.760 (xi=1, eps=1e-9, inside the GCE 1-sigma
  region) the one legacy alpha puts the crossing at 0.77x the forced one:
  <sigma v> 41 % too low.
- Forced engaged at all 56 alphas tested, with no failures, at 30-95 s mean
  per solve (max 265 s).

### Gates

- `phase3_followup_preflight.py` (all through the generated drivers):
  - F1: the validated anchor alpha_relic = 8.386e-3 (xi=2.5, eps=1e-9,
    mX=34.86, rv=0.950) to 2 %.
  - F2: seven stored-crossing Omega/0.12 values reproduce the diagnosis
    to 1 %.
  - F3: the glitch node at alpha=2.065e-3 comes back on the Phase-3 branch
    (Omega 0.308), not the legacy one (0.181).
- `phase3_followup_invalidate.py` is run-once. It checks shape, alpha axis and
  eps first, then writes a marker, so a requeued job resumes instead of
  re-clearing finished nodes.
- Each slice aborts if its driver's `--dry-run` says "starting fresh" or queues
  more nodes than its block.

### Cost

5467 nodes x 8 alphas, about 20-200 s per solve. Expect ~4-6 h of wall clock
on 3 nodes, about 15-25 node-hours.

### After the run

`phase3_followup_report.py` (Stage 4 of the job) compares every file with the
pre-invalidation commit. Expect:
- 0 queued nodes;
- at mX=34.86, alpha_new/alpha_old ~ sqrt(Omega_forced/0.12), e.g. x1.18 at
  xi=2.5 eps=1e-10 rv=0.159;
- the xi=1 eps=1e-9 isolated-node count dropping from 18 toward 0.
