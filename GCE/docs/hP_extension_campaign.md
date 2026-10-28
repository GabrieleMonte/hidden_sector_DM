# Higgs-portal extension campaign — off-shell WW*/ZZ* and mX = 100–130 GeV

> **Status (2026-10-20).** Code for the off-shell WW*/ZZ* channels is committed on this
> branch (`GCE/pythia_runner.py`, `GCE/spectrum.py`, `GCE/scripts/make_channel_spectra.py`).
> Stage 1 is ready to run. Stage 2 needs its driver written first (spec below). Stage 3
> runs on the local machine, **not** here — see "Why the fits stay local".

The hP favoured region (south, free <sigma v>) reaches the mX = 100 GeV edge of the grid
along the mY ~ mX diagonal: Delta chi2 profiled over mY is 5.76 at mX = 100 against a
3-sigma (2-parameter) level of 11.83. Both the 2- and 3-sigma regions must be closed, so
the grid is extended to mX = 130 GeV for the Higgs portal only. Above mY ~ 100 GeV the
scalar decays to W W* / Z Z* (HDECAY: BR(WW*) 1.1% at 100, 5% at 110, 21% at 125, 50% at
140), which the cache could not represent — those channels were skipped below 2 m_V with
no renormalisation. They are now generated.

## Quickstart (on the cluster, after pulling this branch)

```bash
source $HOME/cosmo_env/bin/activate
cd /global/u1/g/gab97/secluded_DM

# Stage 1 — the 28 WW files the main grid now needs (hP points with mY >= 90.379)
# GATE: chanspec.sh writes to $PSCRATCH/channel_spectra unconditionally.  If that staged
# copy was purged, the driver would regenerate all ~157k files.  This MUST print
# "157874/157902 files present, 28 missing" (all WW, mX 91-100) before you submit:
CHANNEL_SPECTRA_DIR=$PSCRATCH/channel_spectra python GCE/scripts/make_channel_spectra.py --check
sbatch -N 1 -q regular -t 01:00:00 GCE/scripts/chanspec.sh      # one node, ~20 min of work
# then copy the 28 back into the repo (chanspec_finish.sh does this, or rsync by hand),
# re-run --check WITHOUT CHANNEL_SPECTRA_DIR (repo copy) until it says 0 missing,
# commit GCE/channel_spectra/WW_*N400000_s12345.npz on this branch and push.
#
# If the $PSCRATCH copy is incomplete, run against the repo copy instead on an
# interactive node (28 files is nothing for DVS):
#   salloc -C cpu -q interactive -N 1 -t 01:00:00 -A m3166
#   python GCE/scripts/make_channel_spectra.py --check        # must say 28 missing
#   python GCE/scripts/make_channel_spectra.py -n 28

# Stage 2 — the hP extension grid (driver not written yet; spec below)
```

Everything below is the why and the spec.

## What is committed (no hand edits needed for Stage 1)

- `pythia_runner.py`: `WIDTHS = {23: 2.4952, 24: 2.085}` and `M_OFFSHELL_MIN = 10.0`
  (Pythia's own 23:mMin = 24:mMin). A channel whose first daughter has a width is a
  V V* channel: `_inject_decay` draws both masses with `_offshell_masses` (exact
  accept-reject from BW(m1^2) BW(m2^2) sqrt(lam)(lam + 12 x1 x2), the CP-even scalar's
  matrix element) and injects the pair exactly like every other channel — same
  isotropic decay, boost, shower scale, hadronisation; Pythia decays each V.
  `min_mY('WW') = 90.379`, `min_mY('ZZ') = 101.188` (one on-shell V + 10 GeV); below
  that HDECAY's WW*/ZZ* share is < 0.2% / 0.1% and stays skipped.
- `spectrum.py`: only the comment on the `min_mY` skip in `cascade_spectrum`.
- `make_channel_spectra.py`: unchanged logic (it already includes WW/ZZ in `CHANNELS`
  and skips below `min_mY`); `--check` and a corrected comment.
- Every other channel's random stream is untouched: bb, tautau, gg came out
  bit-identical to the previous module at the same seed.

## Validation done before committing (local, 2026-10-15 .. 20)

- Sampler vs the exact distribution (per-row quadrature in t = atan((m^2-m0^2)/(m0 w))):
  chi2/bin 0.45-1.52, max |pull| <= 3.2, W at mY = 91-170, Z at 102-130.
- Above threshold (mY = 170): N_gamma per decay within 0.1% (at rest) / 0.6% (boosted)
  of the old fixed-mass on-shell injection; the 2-2.5% <sigma v> shift at BR = 1 is the
  W width that injection ignored.
- Below threshold vs Pythia's own resonance route (inject h0 at mY, 25:onIfMatch = 24 -24):
  W* mass medians agree to 0.1 GeV (22.9/23.0 at 110, 33.1/33.0 at 125, 45.9/45.8 at
  140); N_gamma within 1%; hP-mix Delta chi2 between the two routes -0.02 / -0.11 / +0.49
  at mY = 110 / 125 / 140, i.e. MC noise.
- Effect of including WW* in the hP fit (vs the old silent drop), at mY = 110 / 125 / 140:
  Delta chi2 +0.16 / +0.50 / +1.03, fitted <sigma v> lower by 4% / 18% / 46%
  (the old drop inflated it by 1.04 / 1.22 / 1.84).
- Cost: 1.2-1.4 ms/event for WW/ZZ, but ~2.9 ms/event just above the floor (W at 91-95,
  Z at ~102), where the plain accept-reject keeps few pairs — a 400k file there is
  ~19 min instead of ~9. Accepted deliberately: the public code is what generates the
  public files (a faster, more complex envelope sampler existed and was dropped).

## Stage 2 spec — the hP extension grid (driver to be written)

Decided with the user:
- **Range**: mX 100 -> 130 GeV, Higgs portal only; the shared 16384-point grid and its
  cache are not touched (changing `MX_LIM` would move every node).
- **Columns**: continue the main grid's log spacing, du = log10(100/15)/123 = 0.0066984
  dex: mX = 10^(2 + k du), k = 1..17 -> 101.55 ... 129.98 GeV.
- **Rows**: mY/mX from 0.6 to 1 (mY up to mX, NOT capped at 100), at the main grid's
  density in log mY (its mX = 100 column: 193 points over log10(100/5), dv = 0.0067762
  dex) -> 34 points per column, 578 points. Round both to 2 decimals like the main grid
  (DECIMALS), since the cache filenames carry them.
  Why 0.6: for mX >= 92 the 2- and 3-sigma regions sit at mY/mX >= 0.70, narrowing to
  >= 0.82 at mX = 100.
- **Channels**: the Higgs portal's photon channels only — bb, cc, tautau, gg, ss, mumu,
  gamgam, WW, ZZ (WW at 272 of the points, ZZ at 153). Same N_EVENTS = 400_000,
  SEED = 12345, N_BINS = 180, `default_ebins(mX)` as the main grid.
- **Cost**: ~1 CPU-h per point (400k events at mY ~ 100: bb 6, cc 5, gg 8, tautau 10 min;
  WW/ZZ 9-19 min) -> ~600 CPU-h, a few node-hours on Perlmutter via chanspec.sh-style
  sharding (--worker/--nworkers).
- **Open design choice — ask the user before writing it**: (a) add the extension to
  `make_channel_spectra.py` as configuration constants, so one public driver builds main
  grid + extension and stays resumable, or (b) a separate public driver. Either way it
  must be public code in the repo's style: configuration constants, not CLI flags
  (only -n, --worker/--nworkers, --check say how to run). No patched-up one-off scripts:
  the GCE/ code will be released and must reproduce every file.
- **Check before production**: `--check` count must equal the expected file count; run a
  handful of points first and compare N_gamma per decay with the main grid's mX = 100
  column (continuity across the seam).

Expectation from a b bbar-only proxy along the diagonal (direct bb at m_DM = mX/2, which
ignores WW*): 3-sigma closes near mX ~ 112; with WW* included the full hP mix is already
at Delta chi2 ~ 16 at (mX, mY) = (115.5, 110). 130 GeV is the margin.

## Stage 3 — downstream (hP DONE 2026-10-28; vP/bP/blP still to refit)

1. Copy the 28 Stage-1 WW files and the Stage-2 files from this branch into the local
   `GCE/channel_spectra/` (new files; nothing local is overwritten).  Done: `--check`
   reports 162373/162373 files present, 0 missing, over 16384 main + 578 extension
   points.
2. Re-run the hP free-<sigma v> fit on the main grid (`GCE/scripts/make_GCE_sv_grids.py`):
   the 28 points with mY >= 90.38 now include WW*.  Under 'paper', like everything
   else from now on.  DONE 2026-10-28.  The full refit costs ~2.1 h, but almost none
   of it was needed: `cascade_template` depends on (mX, mY) and the halo alone, never
   on the covariance, so the 16356 templates already stored in the grid are still
   exact (rebuilt and compared bit-for-bit at 327 of them: worst difference 0.000e+00)
   and only had to be run through `gls_fit` again.  The 28 WW* templates were the only
   ones rebuilt.  Minutes, not hours.
3. Fit the extension grid.  DONE 2026-10-28, same pass, 578 points, ~4 min.
   `make_GCE_sv_grids.py` now carries `PORTALS` / `EXT_PORTALS` and writes
   `GCE_sv_grid_hP_ext_south.npz` beside the main file, so this is reproducible from
   the production driver alone.

   RESULT. Best fit moves to chi2 = 15.095 at (mX, mY) = (73.46, 65.85), from 15.603
   at (80.58, 76.89) under the old covariance.  Against it the extension holds no
   point within 1 or 2 sigma and 17 within 3 sigma, and the 3-sigma contour closes at
   **mX = 106.4 GeV** -- 24 GeV inside the grid edge, and below the ~112 the bb-only
   proxy predicted.  2 sigma was already closed at the main grid's own mX = 100 edge
   (min dchi2 there is 8.01, against 6.18).  The seam is smooth: min dchi2 goes 8.01
   at mX = 100 to 8.60 at mX = 101.55, rising monotonically to 26.51 at 129.98.

4. Refit vP, bP and blP under 'paper' too.  DONE 2026-10-28, same template reuse,
   14 s for all three.  They had to be done: `GCE_fits.ipynb` builds its free-<sigma v>
   regions straight from the stored `chi2`/`sigmav`/`sigmav_err` (cells 4, 6, 8, 12-15)
   with nothing recomputed, and `sv_relic_region` in cell 5, although it rebuilds
   tCt/tCphi/phiCphi with the live Cinv, anchors them on `mod_info['chi2'][k_best]`
   from the file.  A stale file therefore biased both kinds of region.

   chi2_min, stored -> 'paper', with the best fit moving:

       vP    18.045 (32.94, 20.01) -> 17.281 (33.97, 22.66)
       bP    18.853 (40.25, 21.60) -> 18.101 (44.16, 30.40)
       blP   18.948 (27.80, 16.88) -> 17.888 (26.95, 16.36)
       hP    15.603 (80.58, 76.89) -> 15.095 (73.46, 65.85)

   The relic regions had been drawn with dchi2 too small by those shifts, so slightly
   too generous.  Inside 3 sigma <sigma v> itself moves by at most 7%; the ratios of
   several hundred that the rebuild reports all sit at dchi2 ~ 111, where the old fit
   had sv_hat ~ 1e-30.  No point is clamped at sv = 0 under 'paper' any more.

   `GCE_fits.ipynb` must be re-run from cell 3 to pick any of this up -- its saved
   outputs predate all four files.  Note also that the extension is a SEPARATE npz, so
   cell 3's "the same grid in every file" still holds; plotting the hP extension needs
   a small addition there, which is the notebook's business, not the driver's.
4. Thermal relic: `GCE/scripts/make_GCE_sv_thermal_relic_grids.py` already
   has the 2026-10-20 coverage fix (RV_MAX 0.995, one coarse row below mY = 5,
   <sigma v> = 0 at mY = mX) and needs a rerun (~2.5 h); its coarse grid stops at
   mX = 100, so the hP relic over 100-130 needs the same extension.

### The covariance prescription (RESOLVED 2026-10-27)
`sigma_mod_scaling='paper'` is the prescription: this branch's `GCE/gce_data.py` /
`GCE/likelihood.py` as of `d4e4419c99` (Sigma_mod x 0.45^2 in the south). The user's
decision; do not change either file back. The paragraph this replaces said the fits had
to stay on the local machine because the two trees disagreed -- that is settled now, and
a fit may run wherever is convenient.

What follows from it. `GCE/output/sv_fits/` was produced by the LOCAL tree under the
rescaled-correlation prescription, so it is superseded: no chi2 read from it is on the
'paper' scale, and the main grid and the extension must be fitted under the same one
before their contours go on one figure. Measured size of the gap at the hP node
(mX, mY) = (80.58, 76.89): chi2 = 15.47 under 'paper' against 15.60 from the
local-produced grid.

Stage 3 items 2-3 are DONE for the Higgs portal as of 2026-10-28; the three vector
portals are not, and are flagged under item 3. The spectra never depended on the
covariance, so Stages 1-2 were always safe here.

## Reproducibility note
The cache is reproducible statistically, not bit-for-bit: `make_channel_spectra.py` builds
one `PythiaRunner(seed=12345)` per worker process and never reseeds between files, and
`imap_unordered` makes the order vary. Per-file seeding from (mX, mY, channel, seed) would
make each file exactly regenerable; it was discussed and deliberately NOT done now.

## Related decisions
- Contour levels: 2.30 / 6.18 / 11.83 (68.3 / 95.4 / 99.73% for 2 parameters). The fit
  notebooks' `9.12` level looks like a typo for 9.21 (99%) — flagged to the user, not
  changed.
- `GCE/notebooks/GCE_fits.ipynb` is the user's combined-portal notebook; do not edit it
  unasked.
