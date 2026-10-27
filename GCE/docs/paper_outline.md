# Paper outline — the GCE from hidden-sector dark matter without thermalization

*Follow-up to arXiv:2605.20186 (Hooper, Krnjaic & Montefalcone, "WIMP-like Dark Matter
Without Thermalization At Freeze-Out"), which closes by naming this analysis as future
work. Re-analysis and extension of arXiv:1912.08821 (Hooper & Witte). PRD format.
Drafted 2026-10-02.*

## What is new relative to the two predecessors

1. **A likelihood replaces a peak-position criterion.** 1912.08821 accepted any model
   whose `E^2 dN/dE` peaked in 1.4-4.0 GeV. We fit the 14-bin SED of arXiv:2112.09706
   (Cholis et al.) against its full Eq. 18 covariance, with `<sigma v>` profiled out in
   closed form.
2. **Favoured regions, not benchmarks.** A 16k-point `(mX, mY)` triangle grid per portal,
   with the `rho_sun` prior marginalised into the chi^2.
3. **The relic condition is imposed through the non-thermalized Boltzmann evolution**, so
   `eps` becomes a label on the favoured region rather than an assumption. This is the
   part 1912.08821 could not do, and the observational statement of the Letter's thesis.
4. **Spectral corrections**: open-flavour thresholds, decay-in-flight treatment, forced
   time shower with per-daughter scale setting.

## Scope decisions (settled)

- **Antiprotons**: gamma-only paper. One discussion paragraph noting the AMS-02 excess,
  that `pbar` spectra are already cached per channel, and that a joint fit is deferred.
- **Portals**: all four (hypercharge, B-L, baryonic, Higgs) carry the full chain —
  GCE region, relic overlay, eps slices, BBN and DD — as in the previous paper. Grids
  currently exist for the hypercharge portal only; the rest is pending compute.
- **Dwarf spheroidals / CMB**: a discussion paragraph in Sec. V, not a dedicated section.

---

# Structure

## I. Introduction

Status of the GCE in 2026; one compressed paragraph on morphology and point-source
systematics with citations, not 1912.08821's two pages. Secluded DM and why cascade
annihilation shifts the favoured masses. What 1912.08821 did and its two limitations
(peak-window acceptance; thermalization assumed). What the Letter established. The four
deliverables above. Roadmap.

## II. Hidden Sector Dark Matter and Portal Models

Short unnumbered preamble before subsection A defining `X`, `Y`, `alpha_X`, `eps` and the
generic `X Xbar -> Y Y -> SM` topology — otherwise the reader meets `eps` in the portal
subsections without knowing what it does.

- **A. Hypercharge Portal**
- **B. B-L Portal**
- **C. Baryonic Portal**
- **D. Higgs Portal**
- **E. Hidden Sector Freeze-Out & Entropy Dilution**

A-D stay tight: Lagrangian, the `Z'`/scalar couplings, the BR structure vs `mY`, and a
single table collecting `g_fV`/`g_fA`, representative BRs and `kappa`. Everything else —
the 2->2 and 3->2 cross sections, `Delta_1/2/3`, `sigma_SI` — is a pointer to the Letter's
supplement and 1912.08821. Target ~2 pages for A-D combined. This section is review; keep
it that way.

E is the recap: the Boltzmann system, the two stages, `xi` through the cannibal phase,
`S_f/S_i`, and `eps_therm` (define it here, once, with the analytic heuristic; it is
needed again in IV B 2 and in V).

**One thing in E must not be relegated:** the GCE sees the `v -> 0` s-wave
`<sigma v>(XX -> YY)`, whereas the relic abundance fixes the thermally averaged value at
freeze-out. Near `mY -> mX` those differ by up to 4.5x, and the whole of Sec. IV B depends
on the distinction. One short paragraph here, recalled once in the results.
See `GCE/docs/sigmav_dip_near_degeneracy.md` Part 1.

## III. The Galactic Center Gamma-Ray Excess

### A. The Excess Emission and Its Systematics

arXiv:2112.09706: the 40x40 deg ROI with the `|b| <= 2 deg` mask, the 14-bin Table III
SED, the 80 interstellar-emission-model variations and the five ranked, and the Eq. 18
covariance `C_ij = sigma_i^2 delta_ij + Sigma_ij,mod^trunc` with `n_pc = 3` (98% of the
variance).

State the point physically, not procedurally: the model systematic is large, strongly
correlated (mean off-diagonal 0.82), and its leading principal component carries 86% of
the variance with the same sign in all 14 bins — exactly the mode a single free
normalisation cannot distinguish from a signal. That is *why* the full covariance is
inverted rather than its diagonal, and it is the main methodological advance over a
peak-position criterion.

State and justify the hemisphere choice (the fit notebooks run `REGION = "south"`), and
show what moves with `both` / `north`.

### B. Generation of the Energy Spectra

Pythia 8.317 hadron-level standalone: one injected two-body `Y -> f fbar` decay, mediator
boosted by `p = sqrt(mX^2 - mY^2)`, forced time shower (silent no-op without
`ev[i].scale(pTmax)` on every daughter), photons histogrammed per decay, 400k events per
channel. BR-weighted summation x2 mediators; the direct arm as the `mY = 2 mX`
(mediator-at-rest) limit of the same code path.

Then the flux assembly: the gNFW profile of their Eq. 6 (`gamma = 1.2`, `r_c = 20` kpc,
`r_sun = 8.5` kpc), the ROI-averaged `Jbar`, and `kappa = 8` (Majorana reference) vs 16
(Dirac portals).

**The open-flavour threshold correction belongs here**, not in Sec. II: `BR(bb)` turns on
at `2 m_B = 10.56` GeV, not `2 m_b = 8.36` GeV, rescaled by the hadronic phase space
(`beta_H(3-beta_H^2)/beta_q(3-beta_q^2)` for spin-1, `(beta_H/beta_q)^3` for the Higgs
portal). It is new work on how the spectra are built, and keeping it out of Sec. II
preserves that section as pure review.

Mediator-spin independence: one sentence citing 1912.08821 Fig. 3, not redone.

### C. Spectral Analysis and the Favoured Cross Section

The chi^2, `chi2 = (phi - s t)^T C^-1 (phi - s t)`; the closed-form GLS profiling of
`<sigma v>`: `s_hat = t C^-1 phi / t C^-1 t`, `sigma_s = (t C^-1 t)^-1/2`, chi^2 exactly
parabolic in `s`, so the entire three-parameter surface is linear algebra — no sampler and
no minimiser anywhere. The `s >= 0` clamp and what it means at the region boundary.
`Delta chi2 = 2.30 / 6.18 / 9.12` for two parameters.

The `rho_sun` prior folded into the chi^2 as a log-normal on `Jbar`
(`0.44 +/- 0.13 GeV/cm^3`, deliberately not the reference 0.40, so it both shifts and
widens); fixed-halo and marginalised contours shown together.

The Monte-Carlo noise floor: the measured per-point Pythia chi^2 scatter (0.35 / 0.16),
why the field is smoothed rather than interpolated, and what size of `Delta chi2`
difference is meaningful.

Closure test: `chi chi -> b bbar` against Fig. 18 of 2112.09706.

*Title alternatives if a shorter one is wanted: "Comparison with the Measured Excess",
"Spectral Fit to the Excess".*

## IV. Results

### A. Favoured Regions from the Excess Alone

The `(mX, mY)` and `(mX, <sigma v>)` regions per portal, `rho_sun` fixed and marginalised.

The physics result is the mass ordering across portals and *why*: the GeV bump is `pi0`
from light-quark and charm fragmentation at `E_q ~ mX/2`, so the tau fraction and the
invisible fraction set where each portal lands (B-L lowest ~35 GeV, hypercharge ~40,
baryonic ~45-55), with the Higgs portal as the control case that climbs back to ~80 GeV
because Yukawa weighting makes it ~85% `b bbar`. See `GCE/docs/gce_mass_range_note.md`.

This subsection is the shorter one and should read as setup for B.

### B. Imposing the Observed Relic Abundance

The headline. Give it room.

1. **Thermal equilibrium with the SM** (`xi = 1`): `alpha_X` at `Omega_c h^2 = 0.12`, the
   induced s-wave `<sigma v>`, and the subset of each region that survives.
2. **Non-thermalized freeze-out and entropy dilution**: the `eps` slices
   (`1e-9, 1e-10, 1e-11, 5e-12`). The result is that the excess is reproducible across
   ~3 decades of `eps`, with the favoured region *moving* as `eps` drops — so the signal
   selects a region but does not pin the portal coupling.

   Within 2, the portal dependence of `Gamma_Y(eps)` needs saying explicitly: the same
   `eps` means different dilution in each portal, or the four panels will look
   inconsistent. The Higgs portal is the sharp case.

   The near-degenerate corner (`mY -> mX`): the s-wave/p-wave split, the decay brake on
   cannibal heating, why the dip appears at `eps = 1e-9` and not at `1e-10`, and the
   solver-reliability caveat at high `r`.

## V. Complementary Constraints and Prospects

BBN hadronic injection (KKMT bound on `m_Y Y_Y` vs lifetime), per portal — the bound
tracks the hadronic BR, so baryonic and Higgs are constrained differently from B-L.

Direct detection: `sigma_SI(eps)` against LZ / XENONnT and the neutrino fog. The headline
is that a GCE-explaining, relic-compatible point can sit orders of magnitude below the
fog. (Watch: the DD helpers return `nan` outside the tabulated `mX` range.)

Then the closing `(mX, eps)` figure carrying the GCE band, relic curve, BBN, DD and
`eps_therm` together — the direct analogue of the Letter's Fig. 3 with the signal on it.

Final paragraph: Fermi dwarf spheroidals (limits sit near the thermal cross section, so
the favoured `<sigma v>` is right at the boundary — quote the number), CMB energy
injection, and antiprotons.

## VI. Summary and Conclusions

## Appendices

- **A.** Boltzmann system, solver phase ladder (joint vs secluded, the `Gamma_Y/H`
  stitching at `GOH_THRESH = 0.05`), and convergence of the relic grids.
- **B.** Monte-Carlo spectra: cache definition, statistical convergence, threshold guards
  (`MIN_ACCEPTED`, `WARN_ACCEPTED`), the decay-in-flight rebinning.
- **C.** Covariance construction and the hemisphere assumption — the multiplicative-in-flux
  rescaling of `Sigma_mod` is a leading-order guess and should be labelled one — with the
  north / south / both comparison of the favoured region.

---

# Figures (main text)

| # | Content | Status |
|---|---|---|
| 1 | Cosmological evolution benchmark (`rho`, `S/S_i`, `xi`), retuned to a GCE-favoured point | needs one re-run |
| 2 | Representative cascade spectra vs the Cholis band, 2x2 by portal | vP current; bP/blP/hP at 200k |
| 3 | `chi chi -> b bbar` closure test vs 2112.09706 Fig. 18 | done |
| 4 | GCE region, `(mX, mY)` + `(mX, <sigma v>)`, 2x2 by portal | vP current; others need the 400k re-run |
| 5 | GCE region + thermalized relic subset, 2x2 | vP only |
| 6 | GCE region at fixed `eps` slices, 2x2 | vP only |
| 7 | `(mX, eps)` with GCE + relic + BBN + DD + `eps_therm`, 2x2 | **not built for any portal** |

Tables: (I) portal couplings, representative BRs, `kappa`; (II) per portal best fit
`(mX, mY, <sigma v>)`, `chi2_min`/dof, 1-sigma ranges.

# Work queue

1. **Fix `decay_width_to_all` below `2 m_B`** in `hidden_sector_DM/model.py`. `Gamma_total`
   still counts partonic `bb` there: 1.2% for vP, 4.8% for bP, **53% for hP**. Blocks hP
   lifetimes, hence hP BBN and hP `eps` slices. Propose as a diff first.
2. Regenerate the bP / blP / hP fit notebooks at `N_EVENTS = 400_000` so all four portals
   sit on one cache.
3. Generalise `make_relic_grid_vP.py` to the other three portals (portal-agnostic apart
   from the hard-coded `VectorPortal` constructor and `ALPHAX_LIMS`).
4. Same for `make_bbn_grid_vP.py`.
5. Build the Fig. 7 assembler.

# Open items

- **Figure load in Sec. IV**: four portals x (region, +thermal relic, +`eps` slices) is
  three four-panel figures plus Figs. 2 and 3. Heavy, but it matches what 1912.08821
  carried, so keep all four portals in the main text rather than pushing three to an
  appendix the way the Letter did.
- Whether the collaborator-comparison material (the `ee` channel and the `pi±`/`K±` decay
  convention, +8% photons) warrants a fourth appendix or a footnote.
