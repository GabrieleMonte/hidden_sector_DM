# Why the Higgs-portal GCE region splits in two at eps = 1e-7

*2026-10-27. Measured from `output/sv_relic_fits/relic_grid_hP_eps1e-0{7,8}.npz` (the
Phase-3 rebuilds, `649bc7e` and `1a3f644`) through `notebooks/GCE_fits.ipynb` at
`c6d72e1`. Companion to `relic_sigmav_regions_note.md`, which covers vP, bP and blP, and
to `sigmav_dip_near_degeneracy.md`, whose mechanism this is.*

- **Settings:** south ROI; xi_ini = 1; 99% c.l. means Delta chi^2 <= 9.21 (2 dof);
  rho_sun profiled over 0.44 +- 0.13 GeV/cm^3.
- **Notation:** r = m_Y/m_X, as in the companion notes.
- **The hP epsilon ladder** is its own: cell 6 sets `EPS_JOINT_hP = (1e-7, 1e-8)` and
  `EPS_SECLUDED_hP = (1e-9, 1e-10)`, against (1e-9, 1e-10) and (1e-11, 5e-12, 1e-12)
  for the vector portals. The Higgs-portal mediator decays through its mixing with the
  Higgs, so the same physics sits two decades higher in eps.

## 1. The observation

At eps = 1e-7 the relic-matched 99% region for the Higgs portal is **two disjoint
pieces**. At eps = 1e-8 it is one. Counted on the native 16384-point fit grid with no
smoothing of any kind:

| r band | eps = 1e-7 | best Delta chi^2 | eps = 1e-8 | best Delta chi^2 |
|---|---|---|---|---|
| 0.20-0.25 | 0 / 1407 | 9.58 | 0 / 1407 | 9.89 |
| 0.25-0.30 | 247 | 4.13 | 226 | 4.67 |
| 0.30-0.55 | 1202 | **3.14** | 1170 | **4.67** |
| 0.55-0.65 | 310 | 3.36 | 196 | 5.58 |
| 0.65-0.70 | 124 | 4.12 | 9 | 8.81 |
| 0.70-0.75 | 63 | 5.75 | 0 | 11.88 |
| 0.75-0.80 | **0** / 496 | 9.90 | 0 | 15.38 |
| 0.80-0.85 | **0** / 496 | 10.32 | 0 | 32.68 |
| 0.85-0.90 | 10 / 496 | 8.53 | 0 | 153.71 |
| 0.90-0.95 | 38 / 372 | 7.23 | 0 | 1306.26 |

1994 of 15759 finite points are inside at 1e-7, in two groups separated by an empty band
at r = 0.75-0.86:

- **Main lobe:** m_X 30.0-79.3, r 0.25-0.75, Delta chi^2_min = 3.14.
- **Island:** m_X 69.1-85.7, r >= 0.86, Delta chi^2_min = 7.23.

On the smoothed field cell 11 actually contours (`relic_mX_mY_smooth_contours`,
SIGMA = 3), `scipy.ndimage.label` finds exactly two components and no speckle: 13187
cells at m_X 30.2-78.5, r 0.25-0.75, min 3.39; and 706 cells at m_X 67.4-83.5, min 7.55.

At eps = 1e-8 there is a single component, m_X 27.4-70.7, r 0.25-0.66,
Delta chi^2_min = 4.90. Its right-hand column above is the point: the fit degrades
monotonically with r, through 11.9, 15.4, 32.7, 154 to 1306. Nothing comes back.

**The split is shallow.** The gap floor is Delta chi^2 = 9.90-10.32 against a 9.21
threshold, so the two lobes separate by about one unit. At 95% they merge. "Disjoint" is
true at 99% and should always be quoted with the level.

## 2. Why: the dip, in the currency the GCE pays in

Relic <sigma v> divided by the GCE-preferred value at the same masses -- the diagnostic
of `relic_sigmav_regions_note.md` section 4 -- at eps = 1e-7:

| m_X | r=0.61 | 0.73 | **0.79** | 0.87 | 0.91 |
|---|---|---|---|---|---|
| 60.1 | 3.1 | 5.9 | **8.0** | 7.3 | 6.6 |
| 70.1 | 2.7 | 4.9 | **7.2** | 6.6 | 5.9 |
| 80.6 | 2.6 | 4.3 | **6.9** | 6.2 | 5.5 |
| 92.6 | 2.5 | 4.0 | **6.6** | 6.0 | 5.2 |

It rises, peaks at r ~ 0.79, and falls -- at every mass. Against the ~4x that the
rho_sun prior can absorb (section 1 of the regions note):

- **r <~ 0.73** the ratio is 2.5-5, affordable: the main lobe.
- **r ~ 0.79** it peaks at 6.6-8.0, more than rho_sun can take: the gap.
- **r ~ 0.91** the dip has pulled it back to 5.2-6.6, and at high m_X that is close
  enough for the region to reopen: the island.

The non-monotonicity is the one `sigmav_dip_near_degeneracy.md` establishes for vP at
eps = 1e-9: as r -> 1 the annihilation migrates from the s-wave (a_v ~ beta^3) into the
p-wave (b_v ~ beta), so the freeze-out rate stays flat while the part that survives to
v -> 0 collapses. It shows through only when the coupling does not have to climb steeply
to hold the relic fixed, which is the case when the mediator already decays faster than
Hubble at freeze-out and brakes the cannibal heating.

**eps = 1e-7 is the Higgs portal's 1e-9.** It is the eps nearest this portal's
thermalisation floor, so it is the one with the decay brake. At eps = 1e-8 the brake is
gone, the ratio climbs monotonically -- 2.8 -> 63.9 across r = 0.30 -> 0.91 at
m_X = 60 -- and the collapse is swamped. No dip, no island. Same mechanism, same
conclusion as the dip note's "the dip belongs to the eps closest to the thermalization
floor".

## 3. Why hP's island survives where vP's dip did not

The regions note records vP at 1e-9 having the same peak-and-dip shape with the dip
bottoming at 4.7x: just outside 99%, so only its left side survives. hP's dip bottoms no
lower in ratio -- 5.2-5.5 at the top -- yet its island is inside. The difference is
**where** the dip lands.

The unconstrained hP GCE best fit is m_X = 80.58, m_Y = 76.89 (r = 0.954),
chi^2 = 15.60 / 11 dof -- the best of the four portals. The island sits on it. The
spectral penalty is smallest exactly where the dip has made the relic value affordable,
and the two combine to get under 9.21. The island is not a corner artefact; it is where
the GCE already prefers the Higgs portal, with the relic line passing close enough at
this one eps for the two to meet.

This is also why the hP panel is drawn on m_X = 15-100 rather than the vector portals'
narrower window (`portal_axes_info['hP']` = mxlims [15, 100], mylims [5, 100]): both
components fall inside it, the island included.

## 4. Caveats

- **r > 0.95 in the smoothed field is extrapolation.** The second component is reported
  up to r = 1.02; the coarse relic grid stops at rv = 0.95 and there is nothing beyond
  it. The native-grid numbers in section 1 do not have this problem.
- **Unsolved nodes sit nearby, but not under the island.** The six eps = 1e-7 nodes that
  never converge are all at rv = 0.95, m_X <= 65.6, below the island's m_X 69.1-85.7.
  The r = 0.90-0.95 band carrying the island has 38 of 372 points inside 99%, so it does
  not rest on a few interpolated nodes.
- **The Phase-3 rebuild did not create this.** Both components, their extents and their
  minima are the same when the quarantined legacy-path grids are substituted
  (`phase3_path_superseded_2026-10-26/`). The topology predates the fix.
- **The ratios use the s-wave <sigma v>**, which near degeneracy understates the
  freeze-out annihilation by up to ~4.5x. That is the right number to compare with the
  GCE and the wrong one to read freeze-out dynamics from; see the dip note's corollary.
- **No portal-specific xi, Gamma_Y/H or dilution was measured for hP here.** Section 2
  applies the vP mechanism to hP on the strength of the measured ratio curve having the
  same shape. The underlying quantities have not been isolated in the solver for this
  portal, exactly as `relic_sigmav_regions_note.md` section 5 flags for bP and blP. This
  is the one genuinely soft step in the argument.

## 5. How this was checked

`GCE_fits.ipynb` was neither edited nor imported. Every code cell was executed in order
in one namespace and the notebook's own functions called, so each number above comes
from the notebook's machinery rather than a reimplementation. Cells 0-7, which produce
all of it, run clean. Two environment substitutions were made and neither touches a
number:

- `plt.savefig` and `Figure.savefig` disabled, so the figure cells cannot overwrite
  `paper_figs/`;
- `plt.style.use` disabled, because `notebooks/mine.mplstyle` is not on the cluster and
  this box has no LaTeX. That leaves `legend.fontsize` at its default string `'medium'`,
  so the figure cells 8-12 raise `TypeError` on `mpl.rcParams['legend.fontsize']*1.1`.
  That is an artefact of running headless, not a defect in those cells. Cell 13
  separately needs `np.trapz`, removed in NumPy 2.

Regions were read two ways: `dchi2_profiled` straight off the 16384-point fit grid for
section 1's counts, and `relic_mX_mY_smooth_contours(fit, SIGMA=3)` -- the call cell 11
makes -- for the component labelling, with `scipy.ndimage.label` (4-connectivity)
discarding groups below 20 cells. The legacy comparison swapped only
`relic_grid_hP_eps1e-0{7,8}.npz` for the quarantined copies, leaving the other 21 files
in `sv_relic_fits/` untouched.

**Known gap, inherited:** as the dip note records, the relic grids store only
`coupling_relic` and `sigmav`, not the Omega h^2(coupling) array they were read from, so
the crossing that places each point of these regions cannot be audited after the fact.
