# Why the relic-compatible ⟨σv⟩ dips near m_Y → m_X, and why only at ε = 10⁻⁹

*Vector (hypercharge) portal. Written 2026-10-02 from `GCE/output/relic/relic_grid_vP_eps*.npz`
plus targeted re-solves; see "How this was checked" at the end.*

## The observation

At m_X = 60 GeV, ε = 10⁻⁹, the s-wave ⟨σv⟩ along the relic line rises to
3.8×10⁻²⁵ cm³/s at m_Y ≈ 48 GeV and then **falls** to 1.4×10⁻²⁵ at m_Y ≈ 57 GeV.
At ε ≤ 10⁻¹⁰ it rises monotonically all the way to the diagonal.

Both behaviours are physical. This is not a solver or interpolation artifact:
doubling the α_X sampling moves the answer by 0.4%, and extending `xmax` 8× moves
it by 0.1% (see below).

> Not to be confused with the **separate** turnover at m_X ≳ 85 GeV, which is the
> Z pole: there m_Y approaches m_Z = 91.19 GeV, Γ_Y explodes, k = ε/ε_therm runs
> from 0.3 to 77, the sector re-locks to the SM and the cannibal enhancement
> switches off. That one is an m_Y effect, not an m_Y/m_X effect, and it appears
> in the ε = 10⁻¹⁰ grid too (starting at m_X = 88 instead of 84.5, as it must).

## Part 1 — The dip is the s-wave dying at threshold, not the annihilation weakening

For `X Xbar -> Y Y` the two partial waves vanish at different rates as the
mediators become non-relativistic. With β = sqrt(1 − m_Y²/m_X²) the mediator
velocity (`model.py:423`, `model.py:436`):

    a_v ∝ β³        b_v ∝ β¹        =>        b_v/a_v ∝ 1/β²

The s-wave carries two extra powers of β. So as m_Y → m_X the annihilation does
not switch off — it **migrates into the p-wave**. At the measured freeze-out
(x_D ≃ 14) the p-wave overtakes the s-wave at β ≃ 0.57, i.e. m_Y/m_X ≃ 0.82.

| m_Y/m_X | 0.727 | 0.795 | 0.869 | 0.909 | 0.950 |
|---|---|---|---|---|---|
| β | 0.687 | 0.607 | 0.495 | 0.418 | 0.312 |
| p-wave share at freeze-out | 35% | 45% | **60%** | 70% | 83% |
| ⟨σv⟩ at freeze-out (x_D=20) | 4.7e-25 | 6.0e-25 | 4.6e-25 | 4.8e-25 | 6.1e-25 |
| ⟨σv⟩ s-wave (what the GCE sees) | 3.4e-25 | **3.8e-25** | 2.3e-25 | 1.8e-25 | **1.4e-25** |

The freeze-out cross section is flat across this range. What collapses by 2.8× is
only the piece that survives to v → 0.

**Why that matters today.** The relic calculation fixes a_v + b_v⟨v²⟩_FO, with
⟨v²⟩_FO ~ 0.4. The GCE sees only a_v: in the halo ⟨v²⟩ ~ 10⁻⁶ makes the p-wave
fraction (b_v/a_v)⟨v²⟩ ≈ 1.2×10⁻⁵ even at m_Y/m_X = 0.95. There is no Sommerfeld
rescue either — the mediator is heavy, α_X m_X/m_Y ≈ α_X ≪ 1.

The pipeline already handles this correctly: `GCE/spectrum.py:181` and
`GCE/scripts/make_relic_grid_vP.py:289` both call `sigmav_XX_to_YY_swave()`, never
the thermally averaged `sigmav_XX_to_YY(T)`. Nothing velocity-dependent reaches
the flux, and the spectrum *shape* is partial-wave blind anyway (same YY final
state) — this is purely about the normalisation.

So the degenerate corner buys its relic abundance in a currency that is worthless
today. The dip is a real suppression of the observable signal and belongs in the
overlay.

## Part 2 — Why only at ε = 10⁻⁹

The β³ collapse is pure kinematics and is identical at every ε. What differs is
**how fast α_X has to climb to keep the relic fixed.**

To hold the s-wave ⟨σv⟩ flat from m_Y/m_X = 0.73 to 0.95, α_X must rise ×2.43.
It rises **×12.8 at ε = 10⁻¹⁰** and only **×1.54 at ε = 10⁻⁹**. One beats the
collapse; the other does not.

Measured at fixed m_X = 60.3 GeV and fixed α_X = 7.455×10⁻³, varying only m_Y:

| | ξ at FO | Γ_Y/H at FO | dilution S | Ωh² growth |
|---|---|---|---|---|
| **ε = 1e-9**, rv 0.73 → 0.95 | 7.0 → 9.4 | **1.8 → 4.6** | 1.00 → 1.00 | **×2.5** |
| **ε = 1e-10**, rv 0.73 → 0.95 | 8.8 → 13.2 | 0.03 → 0.10 | **3.35 → 2.34** | **×25** |

Two mechanisms, both traceable to Γ_Y ∝ m_Y:

**1. The decay brake on the cannibal heating.** Raising m_Y makes the hidden
sector cannibalise over more e-folds — Y goes non-relativistic at x_D = m_X/m_Y,
which moves toward 1 — so ξ at freeze-out rises at both ε. More heating means more
overproduction at fixed coupling, hence more α_X needed. But raising m_Y *also*
speeds up the mediator decay. At ε = 10⁻⁹ the mediator is already decaying faster
than Hubble at freeze-out, and the drain grows from 1.8H to 4.6H exactly as the
heating would: the sector is throttled, ξ rises only ×1.35 and sits 25–40% below
the ε = 10⁻¹⁰ values throughout. At ε = 10⁻¹⁰ the drain is 3–10% of Hubble — no
brake at all, and ξ rises ×1.50 from a higher starting point.

**2. Losing dilution amplifies it, but only at low ε.** At ε = 10⁻¹⁰ the mediator
is long-lived enough to dominate and dilute (S = 3.35). A heavier mediator decays
*earlier*, so the dilution weakens as m_Y → m_X (3.35 → 2.34) — removing a
suppression and adding a further ×1.4 to the overproduction. At ε = 10⁻⁹ there is
no dilution to lose (S = 1.002).

**Net:** at ε = 10⁻¹⁰ the overproduction explodes ×25 across the row and α_X must
climb steeply; at ε = 10⁻⁹ the brake holds it to ×2.5 and α_X barely moves. Only
in the second case does the β³ collapse of the s-wave show through as a dip.

## Consequence for the GCE region

The dip belongs to the ε closest to the thermalization floor — the same regime
that produces the enhanced ⟨σv⟩ in the first place. Near the diagonal the two
effects compound: the hidden sector's heating is braked by the mediator's own
decay, *and* what annihilation remains is increasingly invisible at v → 0. Any
favoured region drawn along m_Y → m_X at ε = 10⁻⁹ inherits both penalties.

A corollary worth keeping in mind: `sigmav` in the relic grids is the s-wave
value, which near degeneracy understates the *freeze-out* annihilation by up to
4.5×. It is the right number to overlay on the GCE contours and the wrong number
to read freeze-out dynamics from.

## How this was checked

Convergence of the grid, m_X = 60.296 GeV, ε = 10⁻⁹, at the two rv columns that
straddle the dip (0.7947 and 0.8689):

| test | rv = 0.7947 | rv = 0.8689 | step |
|---|---|---|---|
| grid file (8 α-knots) | 0.008860 | 0.008471 | ×0.956 |
| 15 α-knots | 0.008863 | 0.008508 | ×0.960 |
| Ωh² at fixed α, xmax = 5e3 | 0.16655 | 0.15736 | ×1.0584 |
| Ωh² at fixed α, xmax = 4e4 | 0.16580 | 0.15691 | ×1.0567 |

Other things ruled out as the cause: Δ₁/Δ₂/Δ₃ and the 2→2 phase space are smooth
and monotone through rv = 0.87; k = ε/ε_therm is smooth there (0.261 → 0.283); the
`sv_handoff_tol = 1.5` coldness gate is inert in this grid because
`make_relic_grid_vP.py` sets `GOH_EXIT_THRESH = 1e99`, which disables `ev_GoH`
(`solver.py:1369`) entirely.

Solver calls matched `omega_h2_joint` exactly: `solve_boltzmann_joint(epsX, ...,
min_eps_floor_ratio=0.0, GoH_exit_thresh=1e99, xmax=5e3, phase2_form="Y")` followed
by `solve_background(..., correct_Y_decay=True, traj_has_decay=True,
stop_on_rhoY=True)`. Note that `SOLVER_KW` in that script is applied only to the
secluded path, so joint solves run at the solver default `rtol_value = 1e-8`.

**Known gap:** the grids store only `alpha_relic` and `sigmav`, not the Ωh²(α_X)
array they were read off, so the spline inversion in `relic_alphaX` cannot be
audited after the fact. Saving that array (or replacing the 8-knot spline read-off
with a Brent root on log α_X) would make this class of question answerable without
re-solving.
