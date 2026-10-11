# Expected DM mass for the GCE fit in the secluded-cascade portals

**Question.** For `X Xbar -> Z' Z' -> 4f` near threshold (`r = mX/mY ~ 1`), where
should the GCE-favoured `mX` sit, and why does it come out ~35-55 GeV for the
vector / baryon / B-L portals instead of the naive ~80-100 GeV?

## Baseline expectation

Near threshold the two `Z'` are produced roughly at rest, so each SM fermion
carries `E_f ~ mX/2`. If the `Z'` decayed to `b bbar`, the cascade spectrum would
just be the direct `b bbar` spectrum at half the energy, and the fit would land
at `2x` the direct `b bbar` mass (~49 GeV) -> ~100 GeV.

## What actually sets the mass

The mediator does **not** decay mainly to `b`. The GCE 1-3 GeV bump is
`pi0 -> gamma gamma` from light-quark + charm fragmentation, a component common to
all portals with nearly identical shape, produced at `E_q ~ mX/2`.

Direct `DM DM -> light qqbar` fits the GCE at `m_DM ~ 24` GeV (charm ~38,
`b bbar` ~49). These are the per-channel best-fit masses from fitting the measured
GCE `E^2 dN/dE` to single-channel prompt spectra: Calore, Cholis, McCabe &
Weniger 2015, "A Tale of Tails", arXiv:1411.4647, Table I. The ordering
`light < c < b` is because a heavier quark gives a softer photon spectrum per unit
parton energy (harder fragmentation function, higher multiplicity, more energy
lost to neutrinos in semileptonic decays), so a larger `m_DM` is needed to push
the photon peak back to ~2 GeV.

Hence the cascade wants `mX/2 ~ 25-30` -> `mX ~ 40-60`, roughly
portal-independent. The naive `x2` (4-body final state) and a `/2` (light/charm
vs `b`) cancel.

## Taus split the portals

`tau+ tau-` is the hardest channel (direct GCE fit ~10 GeV). The larger the tau
fraction of the *visible* final state, the softer the spectrum per unit energy,
and the lower the preferred `mX`.

| portal      | Z' coupling            | ~`bb` / `tautau` / invisible-ish   | favoured `mX` |
|-------------|------------------------|------------------------------------|---------------|
| B-L (blP)   | `~ (B-L)_f`            | 4% / 13% / ~60% (`nu`, `e/mu`)     | lowest (~35)  |
| vector (vP) | `~ Q_f` (charge)       | 5% / 15% / ~30% (`e/mu`)           | ~40           |
| baryon (bP) | `~ B_f`, quarks only   | 20% / 0 / 0                        | ~45-55        |
| Higgs (hP)  | `~ y_f` (Yukawa)       | ~85% / ~8% / 0                     | **~80 (confirmed by run)** |

Branching fractions are approximate, evaluated near `mY ~ 40` GeV from the
coupling structure in `hidden_sector_DM/model.py` (`VectorPortal.gfv_gfa` reduces
to `~ e Q_f` for `mZ' << mZ`; bP is flavour-democratic over quarks; blP adds
charged leptons and neutrinos at `(B-L) = -1`).

The **Higgs portal is the control case**: Yukawa weighting makes it ~85%
`b bbar`, so the "cascade to `b bbar`" assumption holds and it does climb to
~80 GeV, as the baseline argument predicts. That the other three do not is the
final-state composition effect, not a fitting artefact.

## Normalisation aside

`e/mu` (vP, blP) and `nu` (blP) carry off annihilation energy without producing
GCE photons, so those portals need a larger `<sigma v>` than bP/hP to match the
GCE amplitude. The portals differ more in favoured `<sigma v>` than in favoured
`mX`.

## Caveat

The peak-position / band-flux test resolves `mX` only to a factor ~1.5-2 with
`mY` and `<sigma v>` floating, so the favoured bands are broad and overlapping.
The mass ordering above is real, but the band edges are soft.
