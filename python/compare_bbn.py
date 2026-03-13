"""
Compare the three BBN constraint methods for the VectorPortal
with mX = 50 GeV, mY = mX / 1.5, alphaX = 1e-3.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.style.use('/home/gab/Desktop/PyCharm_env/mine.mplstyle')

from source import Cosmology, BoltzmannSolver
from source.model import VectorPortal

hbar = 6.582119514e-25  # GeV*s

# ── Model parameters ──
mX = 50.0
r  = 1.5
mY = mX / r
alphaX = 1e-3

print(f"mX = {mX:.1f} GeV,  mY = {mY:.2f} GeV,  r = {r},  alphaX = {alphaX}")
print("=" * 60)

# ── Build model and solve cannibal phase ──
cosmo = Cosmology(gstar_choice="standard")
model = VectorPortal(
    mX=mX, mY=mY, gX=2, gY=3,
    alphaX=alphaX,
    include_antiparticlesX=True,
    include_antiparticlesY=False,
)
solver = BoltzmannSolver(cosmo, model)

print("Running cannibal-phase solver...")
sol0 = solver.solve_boltzmann_chempot_QSSA(
    return_bg_ICs=True, verbose=False,
    Gamma_switch_QSSA=5e9, Gamma_switch_full=1e4,
    n_points=500, rtol_value=1e-8, atol_value=1e-18,
    convergence_threshold=5e-3,
)
bg_ICs = sol0['bg_ICs']
YY_frozen = sol0['YY_final']
YX_relic  = sol0['YX_relic']
print(f"YX_relic  = {YX_relic:.4e}")
print(f"YY_frozen = {YY_frozen:.4e}")
print(f"mY * YY   = {mY * YY_frozen:.4e} GeV")
print()

# ── Method 1: find_epsilon_BBN (temperature-based) ──
print("Method 1: find_epsilon_BBN (Y must decay before T_BBN = 10 MeV)...")
try:
    eps_T = solver.find_epsilon_BBN(
        bg_ICs, T_BBN=10e-3,
        log10eps_range=(-18, -6), root_rtol=1e-3,
    )
    tau_T = hbar / model.decay_width_to_SM(eps_T)
    print(f"  eps = {eps_T:.4e}   (log10 = {np.log10(eps_T):.2f})")
    print(f"  tau = {tau_T:.4e} s")
except Exception as e:
    print(f"  FAILED: {e}")
    eps_T = None

# ── Method 2: find_epsilon_BBN_Neff ──
print("\nMethod 2: find_epsilon_BBN_Neff (DNeff < 0.151 at T = 1 MeV)...")
try:
    eps_Neff = solver.find_epsilon_BBN_Neff(
        bg_ICs, T_eval=1e-3, DNeff_max=0.151,
        log10eps_range=(-18, -6), root_rtol=1e-3,
    )
    tau_Neff = hbar / model.decay_width_to_SM(eps_Neff)
    print(f"  eps = {eps_Neff:.4e}   (log10 = {np.log10(eps_Neff):.2f})")
    print(f"  tau = {tau_Neff:.4e} s")
except Exception as e:
    print(f"  FAILED: {e}")
    eps_Neff = None

# ── Method 3: find_epsilon_BBN_hadronic ──
print("\nMethod 3: find_epsilon_BBN_hadronic (Kawasaki et al. injection)...")
try:
    eps_had = solver.find_epsilon_BBN_hadronic(
        YY_frozen, log10eps_range=(-18, -6), root_rtol=1e-3,
    )
    tau_had = hbar / model.decay_width_to_SM(eps_had)
    print(f"  eps = {eps_had:.4e}   (log10 = {np.log10(eps_had):.2f})")
    print(f"  tau = {tau_had:.4e} s")
except Exception as e:
    print(f"  FAILED: {e}")
    eps_had = None

# ── Branching ratios for context ──
br = model.branching_ratios_to_SM()
br_had_frac = sum(v for ch, v in br.items()
                  if ch in {'uu','dd','ss','cc','bb','tt','gg','WW','ZZ','hh'})
print(f"\nBranching ratios:  BR_had = {br_had_frac:.3f}")
for ch, v in sorted(br.items(), key=lambda x: -x[1]):
    if v > 0.005:
        print(f"  {ch:>10s}: {v:.4f}")

# ── Plot: constraint strength vs eps ──
print("\nGenerating comparison plot...")

eps_arr = np.logspace(-16, -8, 300)
tau_arr = hbar / np.array([model.decay_width_to_SM(e) for e in eps_arr])

# For hadronic method: compute mY_bound(tau) for each eps
from source.bbn import mY_bound_model
mY_YY = mY * YY_frozen
bounds_had = np.array([mY_bound_model(mY, t, br) for t in tau_arr])
# ratio: mY*YY / bound.  >1 means excluded
ratio_had = mY_YY / bounds_had

fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex=True,
                         gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.08})

# --- Top panel: lifetime vs eps ---
ax = axes[0]
ax.loglog(eps_arr, tau_arr, 'k-', lw=1.5, label=r'$\tau(\varepsilon)$')
ax.axhspan(0.2, 50, color='red', alpha=0.10, label='BBN window (0.2–50 s)')

# Mark the three eps values
markers = []
if eps_T is not None:
    ax.axvline(eps_T, color='C0', ls='--', lw=1.2, label=rf'$T$-based: $\varepsilon={eps_T:.1e}$')
    markers.append((eps_T, tau_T, 'C0'))
if eps_Neff is not None:
    ax.axvline(eps_Neff, color='C1', ls='-.', lw=1.2, label=rf'$\Delta N_{{\rm eff}}$: $\varepsilon={eps_Neff:.1e}$')
    markers.append((eps_Neff, tau_Neff, 'C1'))
if eps_had is not None:
    ax.axvline(eps_had, color='C3', ls=':', lw=1.5, label=rf'Hadronic inj.: $\varepsilon={eps_had:.1e}$')
    markers.append((eps_had, tau_had, 'C3'))

for e, t, c in markers:
    ax.plot(e, t, 'o', color=c, ms=7, zorder=5)

ax.set_ylabel(r'$\tau_Y$ [s]')
ax.set_ylim(1e-8, 1e12)
ax.legend(fontsize=8, loc='upper right')
ax.set_title(
    rf'VectorPortal:  $m_X={mX:.0f}$ GeV,  $m_Y={mY:.1f}$ GeV,  '
    rf'$\alpha_X={alphaX}$,  $Y_Y={YY_frozen:.2e}$',
    fontsize=10,
)

# --- Bottom panel: BBN exclusion ratio ---
ax2 = axes[1]
ax2.loglog(eps_arr, ratio_had, 'C3-', lw=1.5,
           label=r'$m_Y Y_Y \,/\, (m_Y Y_Y)_{\rm bound}^{\rm had}$')
ax2.axhline(1, color='k', ls='-', lw=0.8)
ax2.fill_between(eps_arr, 1, ratio_had, where=(ratio_had > 1),
                 color='C3', alpha=0.15, label='Excluded (hadronic)')

if eps_T is not None:
    ax2.axvline(eps_T, color='C0', ls='--', lw=1.0)
if eps_Neff is not None:
    ax2.axvline(eps_Neff, color='C1', ls='-.', lw=1.0)
if eps_had is not None:
    ax2.axvline(eps_had, color='C3', ls=':', lw=1.5)

ax2.set_xlabel(r'$\varepsilon$')
ax2.set_ylabel(r'$m_Y Y_Y \,/\, {\rm bound}$')
ax2.set_ylim(1e-4, 1e4)
ax2.legend(fontsize=8, loc='upper right')

fig.savefig('compare_bbn_methods.pdf', bbox_inches='tight')
print("Saved: compare_bbn_methods.pdf")
