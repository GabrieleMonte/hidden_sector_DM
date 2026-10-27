# What the relic-matched GCE regions show, and why (vP, bP, blP)

*2026-10-22.* The figures are "Regions with the observed relic abundance" in
`notebooks/GCE_fits.ipynb` (`paper_figs/GCE_fits_{vP,bP,blP}_secluded_v0`).

- **Settings:** south ROI; initial temperature ratio xi_ini = 1; 99% c.l. means
  Delta chi^2 <= 9.21 (2 dof); rho_sun profiled over 0.44 +- 0.13 GeV/cm^3.
- **Notation:** r = m_Z'/m_DM throughout.
- **Sources:** numbers come from the current `output/sv_fits` and `output/sv_relic_fits`
  files unless a source is named. Measurements quoted from `sigmav_dip_near_degeneracy.md`
  are for vP at m_DM = 60 GeV only; section 5 says where they are extrapolated.

## 1. What a panel shows

- **Top row:** the GCE Delta chi^2 when <sigma v> is fixed, at each (m_DM, m_Z'), to the
  value that gives Omega h^2 = 0.12.
- **Bottom row:** that <sigma v> inside the 99% region.

On its own, the GCE fit prefers <sigma v>_GCE = 3.4 (vP), 2.6 (bP) and 6.4 (blP)
x 1e-26 cm^3/s at its best fit (rho_sun = 0.44). The preferred value also moves somewhat
with the masses.

The gamma-ray signal scales as <sigma v> x rho_sun^2. A relic value k times too large can
therefore be compensated by lowering rho_sun by sqrt(k), at a cost from the rho_sun prior.
In practice a relic value up to roughly 4x the GCE value still lands inside 99%, with
rho_sun pulled to ~0.2-0.25. At 4.7x it already does not (vP's diagonal at eps = 1e-9;
numbers in that section). This is why the ends of each colour scale always sit on a region edge.

## 2. What sets the relic <sigma v>

The bottom row shows *today's* annihilation rate. It is the product of two things:
(a) how strong the dark coupling alpha_X must be to leave Omega h^2 = 0.12, which is set by
the cosmology, and (b) how much of that coupling's annihilation still operates at today's
halo velocities.

**(a) The cosmology.** Four effects, whose weight changes with eps:

1. **Thermal contact** (eps >= eps_therm).
   - *Mechanism:* the mediator decays into SM particles and is re-created from them fast
     enough that the dark sector stays at the SM temperature. Every DM annihilation
     therefore ends up as SM particles: ordinary freeze-out.
   - *Consequence:* the annihilation rate at freeze-out must be about 4e-26 cm^3/s (Dirac
     DM). Nothing here depends on how the mediator couples to the SM, so the value is
     identical for vP, bP and blP.

2. **A hotter dark sector** (eps < eps_therm).
   - *Mechanism:* without contact, once the mediators become non-relativistic, reactions
     in which three particles become two (Z'Z'Z' -> Z'Z', and the same with a DM
     particle) keep converting mediator rest mass into kinetic energy. The dark sector
     then cools much more slowly than the SM, and xi = T_dark/T_SM grows.
   - *Why that costs <sigma v>:* DM freezes out at a fixed *dark* temperature
     (about m_DM/20). The hotter the dark sector is relative to the SM at that moment, the
     more DM is left per unit of SM entropy. To first approximation the required
     <sigma v> grows in proportion to xi.
   - *Measured:* xi ~ 7-9 at freeze-out for eps = 1e-9 and 9-13 for eps = 1e-10 (r 0.73-0.95).
   - *Dependence on r:* a heavier mediator becomes non-relativistic earlier, so the heating
     runs longer before DM freezes out. xi grows with r.

3. **Can the dark sector shed particles during freeze-out?** Annihilation X X -> Z' Z' only
   lowers the DM abundance if the mediators do not turn back into DM. Two cases:
   - *The mediators decay faster than the universe expands* (Gamma/H > 1). They then leave
     the dark sector and carry energy away, which also slows the heating of point 2.
     Measured Gamma/H at DM freeze-out: 1.8-4.6 at eps = 1e-9 and 0.03-0.10 at 1e-10
     (r 0.73-0.95). It is smaller still at 1e-11.
   - *They do not decay in time.* They then stay away from the DM only if they lack the
     energy to make it: the reverse reaction Z' Z' -> X X is suppressed by roughly
     e^{-(1-r) m_DM/T_dark}.
     - For a light mediator (small r) that factor is tiny, and annihilation removes DM
       normally.
     - As r -> 1 the suppression disappears. DM and mediators keep converting into each
       other, so DM can only fall as fast as the whole dark sector shrinks, which only the
       three-to-two reactions do. The relic then responds weakly to the coupling:
       measured Omega ~ alpha_X^-1.65 at r = 0.4 and alpha_X^-1.14 at r = 0.9, against
       alpha_X^-2 for ordinary freeze-out. A much larger alpha_X is needed.

4. **Dilution afterwards.**
   - *Mechanism:* if the mediators live long enough to dominate the energy of the universe
     before decaying, their decays heat the SM and dilute the DM by a factor D, so a
     smaller <sigma v> suffices. A shorter lifetime means less domination and less
     dilution.
   - *Measured:* D = 1.00 at 1e-9, and 3.4 -> 2.3 at 1e-10 over r 0.73 -> 0.95 (a heavier
     mediator decays earlier). In the earlier secluded runs at 1e-11, D was 1e2-1e3.

**(b) What still operates today.**
- *Two partial waves.* DM annihilates in an s-wave, which works at any velocity, and a
  p-wave, which needs the DM to be moving. The p-wave matters at freeze-out (v ~ 0.3) but is
  negligible in today's halo (v ~ 1e-3).
- *Threshold behaviour.* Near r = 1 the mediators come out almost at rest. The s-wave then
  falls as beta^3 and the p-wave only as beta, with beta = sqrt(1-r^2) their speed. So as
  r -> 1, a growing share of the freeze-out annihilation is p-wave, which does nothing
  today.
- *Size of the effect.* The p-wave share is ~8% at r = 0.3 and ~60% at r = 0.9 in the
  thermal case. At eps = 1e-9 it grows from 35% to 83% between r = 0.73 and 0.95.

## 3. Why the three portals differ

The three portals share the same dark sector; at a given eps they differ only in how fast
the mediator decays. eps is normalised differently in each (kinetic mixing for vP, eps*g
for bP and blP), so the same eps means a different lifetime. At m_Z' = 17-21 GeV:

- Gamma_bP ~ 3.2 Gamma_vP (the baryon mediator couples only to quarks);
- Gamma_blP ~ 15 Gamma_vP (the B-L mediator couples to every fermion, neutrinos included).

Since Gamma ~ eps^2, blP at a given eps behaves roughly like vP at 3.9 eps, and bP like vP
at 1.8 eps. The portals also differ in their GCE fit: the preferred <sigma v>, and which
masses give a good spectrum.

## 4. Panel by panel

<sigma v>(r) at m_DM = 33 GeV, in 1e-26 cm^3/s (vP / bP / blP). It depends very little on
m_DM between 28 and 40 GeV:

| r | thermal | eps = 1e-9 | eps = 1e-10 | eps = 1e-11 |
|---|---|---|---|---|
| 0.30 | 4.1 (all) | 7.5 / 7.8 / 7.0 | 2.8 / 4.2 / 6.1 | 0.26 / 0.50 / 1.2 |
| 0.61 | 3.9 | 26 / 19 / 11 | 16 / 23 / 29 | 1.0 / 2.2 / 5.9 |
| 0.79 | 3.6 | 46 / 24 / 8.0 | 102 / 116 / 94 | 7.1 / 17 / 48 |
| 0.91 | 2.6 | 26 / 13 / 5.1 | 373 / 245 / 86 | 78 / 160 / 144 |

Best point of each panel: Delta chi^2_min with rho_sun profiled (rho_sun fixed at 0.44 in
brackets), the profiled rho_sun there, and the r range of the 99% region:

| | vP | bP | blP |
|---|---|---|---|
| thermal | 0.06 (0.23), 0.41, r 0.22-0.97 | 0.26 (0.69), 0.38, r 0.25-0.97 | 0.50 (1.70), 0.51, r 0.22-0.95 |
| 1e-9 | 6.65 (31.8), 0.22, r 0.26-0.54 | 6.97 (138), 0.23, r 0.33-0.94 | 0.57 (0.96), 0.39, r 0.22-0.94 |
| 1e-10 | 1.78 (8.17), 0.32, r 0.21-0.65 | 4.01 (36.8), 0.26, r 0.23-0.55 | 2.05 (12.5), 0.31, r 0.22-0.67 |
| 1e-11 | 0.32 (0.47), 0.42, r 0.46-0.86 | 0.05 (0.17), 0.45, r 0.22-0.79 | 0.02 (0.08), 0.46, r 0.24-0.78 |

### Thermal (eps >= eps_therm)

**What we see.** <sigma v> is 4.1e-26 at small r, falls to 2.6e-26 at r = 0.9 and to 0.8e-26
at r = 0.95. It is identical for the three portals, and the regions nearly fill the
GCE-only region.

**Why.**
- *Small r:* this is ordinary freeze-out (point 1). The annihilation is almost entirely
  s-wave, so today's rate equals the freeze-out rate, ~4e-26.
- *Towards the diagonal:* the freeze-out requirement barely changes, but more of it is
  p-wave (point b, 8% -> 60%). The coupling grows only modestly to compensate (alpha_X x1.6
  from r = 0.3 to 0.91), so today's rate falls.
- *Regions:* the thermal value is already close to what the GCE wants, so nearly all of the
  GCE region survives. The portal differences come only from the spectral fit. blP's GCE
  fit wants more than thermal (6.4e-26), so its best point raises rho_sun to 0.51.

### eps = 1e-9

**What we see.** For vP and bP, <sigma v> is above thermal everywhere. It rises to a peak
near r = 0.8 (vP 46e-26, bP 24e-26) and falls towards the diagonal (vP 21e-26, bP 10.5e-26
at r = 0.95). blP stays close to thermal (7-11.5e-26).

**Rising side, r < 0.8 — heating grows with r.**
- The dark sector is out of thermal contact, so it heats (point 2), and the heating grows
  with r because the mediator turns non-relativistic earlier.
- At small r, the mediator is light: annihilation removes DM normally and is almost pure
  s-wave. So today's rate equals the freeze-out rate, which is above thermal only by the
  modest heating: 1.8x at r = 0.3.
- By r = 0.79 the heating has had much longer to act, and the requirement is 13x thermal.
- Mediator decays happen promptly, so they never dominate and nothing dilutes the DM
  (D = 1.00).

**Falling side, r > 0.8 — the decays brake the heating and the annihilation moves into the
p-wave.**
- A heavier mediator decays faster (Gamma ~ m_Z'), and here it already decays faster than
  the expansion. The decays drain the dark sector: xi grows only x1.35 from r = 0.73 to
  0.95, against x1.5 at 1e-10.
- The required freeze-out rate therefore stops growing (flat at 4.6-6.1e-25 at m_DM 60),
  and alpha_X barely moves (x1.54).
- With the freeze-out requirement flat, today's rate follows the shrinking s-wave share
  (point b, 65% -> 17% of the annihilation). It falls x2.8 between r = 0.80 and 0.95 at
  m_DM 60, and from 46 to 21e-26 at m_DM 33.

**Portal differences.**
- A faster-decaying mediator drains more heat, so it needs less <sigma v>. The effect is
  larger where there is more heating to drain: bP/vP = 1.0 at r = 0.3, 0.75 at 0.61, 0.52
  at 0.79, 0.49 at 0.95.
- blP (15x vP's width) drains so fast that the dark sector stays close to the SM
  temperature. Its <sigma v> is nearly flat and near thermal, close to its preferred
  6.4e-26. It gets a large region and the best 1e-9 fit (0.57).

**The regions: why bP has two lobes and vP one.** Relic value divided by the GCE-preferred
value at the same masses (best case over the masses the GCE spectrum allows):

| r | 0.30-0.50 | 0.65-0.80 | 0.90-0.95 |
|---|---|---|---|
| vP | 2.8-3.7, inside | 7.3-9.6, outside | 4.7, outside |
| bP | 3.6-4.5, inside (left lobe) | 6.9-7.7, outside (gap) | 3.6, inside (right lobe) |

- **bP's left lobe** (r 0.30-0.55, m_DM 31-43): <sigma v> 7-16e-26, low because the heating
  has had little time to act.
- **bP's gap** (r 0.55-0.85): the heating is strongest and <sigma v> reaches 15-33e-26,
  about 7x what bP's spectrum wants. That is more than rho_sun can absorb.
- **bP's right lobe** (r 0.85-0.95, m_DM 46-61): the decays cap the heating and most of the
  freeze-out annihilation is p-wave, so today's value drops back to 12-22e-26 (3.6x). It
  holds the overall minimum (6.97).
- Both lobes need rho_sun ~ 0.23 (2.2 sigma below the prior mean).
- **vP** has the same peak-and-dip shape. Its diagonal values are about twice bP's (slower
  decays, weaker brake), and its GCE value is higher, so the dip only gets down to 4.7x:
  just outside 99%. Only its left side survives (r 0.26-0.54, rho_sun 0.22, Delta chi^2_min
  6.65).

### eps = 1e-10

**What we see.** <sigma v> is below the 1e-9 value at small r (vP 2.8 vs 7.5e-26 at
r = 0.3) and far above it at large r (373 vs 26e-26 at r = 0.91). It rises all the way to
the diagonal, with no peak; the 1e-9 and 1e-10 curves cross near r = 0.7.

**Why.** The decays are 100x slower than at 1e-9 (Gamma ~ eps^2). During freeze-out they
do nothing (Gamma/H = 0.03-0.10), so the heating runs unchecked: xi is 9-13 here against
7-9 at 1e-9. After freeze-out, the long-lived mediators briefly dominate and their decays
dilute the DM by D ~ 2-3.

- **Small r — dilution wins.** The mediator is light, so annihilation removes DM normally,
  and the heating is still modest. The main difference from 1e-9 is the dilution by ~3,
  which lowers the requirement by about that factor (2.8 vs 7.5e-26).
- **Large r — everything works against it.** Three things compound:
  - the heating grows unchecked;
  - DM and mediators keep converting into each other, so only the three-to-two reactions
    deplete the sector (point 3);
  - the dilution weakens (3.4 -> 2.3), because a heavier mediator decays sooner.

  alpha_X has to climb x12.8 from r = 0.73 to 0.95 at m_DM 60. That is far faster than the
  s-wave share falls, so today's rate keeps rising and there is no dip.

**Portal differences.** Faster decays mean less dilution, which raises the requirement at
small r (blP 6.1 > bP 4.2 > vP 2.8). They also give some draining during freeze-out, which
lowers it at large r (vP 373 > bP 245 > blP 86).

**Regions.** Only r <~ 0.55-0.67 is within reach, with rho_sun pulled to 0.26-0.32. bP fits
worst (4.0) because its preferred <sigma v> is the lowest.

### eps = 1e-11 (secluded solver)

**What we see.** <sigma v> is far below thermal at small r (vP 0.26e-26 at r = 0.3) and
rises steeply with r. At m_DM 33 it crosses the thermal value near r = 0.73 (vP), 0.66 (bP)
and 0.55 (blP), and approaches 1e-24 near the diagonal. The regions are diagonal bands.

**Why.** The mediators live long enough (for vP with m_Z' = 20 GeV they decay at
T_SM ~ 7 MeV; estimate) that the dark sector is completely isolated during freeze-out.
Afterwards the mediators dominate the universe for a long time, and their decays dilute
the DM by D ~ 1e2-1e3.

- **Small r:** freeze-out proceeds as at 1e-10, with DM removed normally, and is then
  followed by a large dilution. The required <sigma v> is ~16x below thermal.
- **Large r:** the mediators are as abundant as the DM (point 3) and nothing helps during
  freeze-out, so alpha_X must climb steeply. The dilution is set by the lifetime and hardly
  changes, so it can no longer keep up. <sigma v> crosses thermal and keeps rising.

**Why the regions are bands.** <sigma v> depends almost only on r: at r = 0.61 it is
1.04 / 1.01 / 0.98e-26 at m_DM = 28 / 33 / 40. Matching the GCE value therefore selects a
window in r, i.e. a wedge between two lines m_Z' = r m_DM. The colour runs from ~4e-27 at
the bottom edge to ~2-4e-25 at the top edge, because the rho_sun freedom (about 4x either
way) spans a wide r window on such a steep curve.

**Portal differences.** A faster-decaying mediator dominates for less time, so the DM is
diluted less and more <sigma v> is needed at a given r. At r = 0.3, blP : bP : vP =
4.5 : 1.9 : 1, which is what dilution ~ 1/sqrt(Gamma) predicts from the widths. The bands
of bP and blP therefore sit at smaller r (0.22-0.79 and 0.24-0.78) than vP's (0.46-0.86).
These are the best relic-matched fits of all (Delta chi^2_min 0.02-0.32), and rho_sun stays
at 0.42-0.46: the relic matches the GCE without help from the halo.

**Two consistency checks.** If mediator decays had already started during freeze-out at
1e-11, they would act as an extra particle sink and pull the 1e-11 values below both
predictions below. Neither shows that.
- *Across portals:* the ratios above (4.5 and 1.9 at r = 0.3) match the width prediction.
- *Against 1e-12* (10x more dilution): the requirement drops 13x at r = 0.3 and 62x at
  r = 0.91. That is what the measured weak coupling dependence of point 3 predicts
  (exponents 1.65 and 1.14).

## 5. Assumptions and open items

- **Extrapolated measurements.** xi, Gamma/H, D, the p-wave shares and the coupling
  exponents were measured for vP at m_DM = 60 and r >= 0.73 (D at 1e-11 comes from earlier
  secluded runs). Here they are applied at m_DM = 25-55, at small r, and to bP and blP via
  their widths.
  - The statement that the heating is weak at small r is physical reasoning (the mediator
    becomes non-relativistic later), not a measured xi.
  - bP's and blP's xi and D have not been measured at all.
- **"Required <sigma v> ~ xi"** is the standard scaling for a decoupled sector. It ignores
  the shift in the dark freeze-out temperature, so treat it as a trend, not a formula.
- **Point 3's near-degeneracy picture** is our explanation of the measured weak coupling
  dependence. It has not been isolated in the solver.
- **The top edge of the vP eps = 1e-11 band (r up to 0.86).** It relies on the secluded
  solver's assumption that the mediator is stable through freeze-out. That assumption fails
  if Gamma/H exceeds 0.05, i.e. once T_SM <~ 40 MeV for vP.
  - The checks above cover r = 0.79-0.91 and show no sign of a problem.
  - One secluded solve at (m_DM 33.4, r 0.87) would settle it directly.
- **Outside every plotted region** (they do not affect the figures):
  - eps = 1e-11 at r = 0.95 is erratic across m_DM and needs alpha_X up to ~0.6: the known
    folded-Omega / non-perturbative corner. The 1e-11 regions end at r <= 0.86.
  - One bP eps = 1e-12 node gives <sigma v> ~ 0 (1e-12 is not plotted).
  - eps = 1e-10 at r >~ 0.73 was flagged in earlier audits. The 1e-10 regions end at
    r <= 0.67.
- **Width ratios** were measured at m_Z' = 17-21 GeV only. bP couples only to quarks, so its
  ratio changes across quark thresholds.
