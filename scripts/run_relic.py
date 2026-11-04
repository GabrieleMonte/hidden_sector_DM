"""Single-point relic-abundance example: one secluded freeze-out and its late Y decay.

    cosmo_env
    python scripts/run_relic.py

1. Solve the hidden-sector freeze-out for a vector-portal benchmark (Y stable).
2. Find the portal coupling eps that gives Omega_c h^2 = 0.12 through the entropy
   injected by the late Y -> SM decay.
3. Evolve the background at that eps and plot the energy budget rho_j/rho_tot and the
   entropy S/S_i against the SM temperature.

Masses in GeV.  Output: the figure, shown on screen.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import LogLocator, NullFormatter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from hidden_sector_DM.HiddenSectorDM import BoltzmannSolver, Cosmology, VectorPortal  # noqa: E402

plt.style.use(ROOT / 'notebooks' / 'mine.mplstyle')

# ---- configuration -------------------------------------------------------

MXV = 50.0          # DM mass [GeV]
MYV = 20.0          # mediator mass [GeV]
ALPHAXV = 1e-4      # dark coupling alpha_X

# ---- 1. freeze-out, Y stable ---------------------------------------------

cosmo = Cosmology()
model = VectorPortal(mX=MXV, mY=MYV, gX=2, gY=3, alphaX=ALPHAXV,
                     include_antiparticlesX=True, include_antiparticlesY=False)
solver = BoltzmannSolver(cosmo, model)
sol0 = solver.solve_boltzmann_chempot_2phase(return_bg_ICs=True, verbose=True,
                                             Gamma_switch_threshold=1e6, n_points=500)

# ---- 2. eps for Omega_c h^2 = 0.12 ---------------------------------------

epsX = solver.find_epsilon_DMRelic(sol0['bg_ICs'])
print(f'eps for Omega_c h^2 = 0.12: {epsX:.4e}')

# ---- 3. background evolution at that eps ---------------------------------

sol_bg = solver.solve_background(epsX=epsX, bg_ICs=sol0['bg_ICs'], stop_on_rhoY=False,
                                 verbose=False)

fig1 = plt.figure(figsize=(4.5, 4.5))
gs1 = GridSpec(21, 23, figure=fig1)
ax1 = fig1.add_subplot(gs1[:13, :])
ax2 = fig1.add_subplot(gs1[13:, :])

# --- Energy budget ---
ax1.plot(sol_bg['T'], sol_bg['rhoR']/sol_bg['rhoTot'], lw=1, c='tab:red', label=r'SM')
ax1.plot(sol_bg['T'], sol_bg['rhoY']/sol_bg['rhoTot'], lw=1, c='tab:blue', label=r'$Y$')
ax1.plot(sol_bg['T'], sol_bg['rhoX']/sol_bg['rhoTot'], lw=1, c='tab:purple', label=r'$X$')
ax1.axvline(1e-2, ls='--', c='tab:gray', zorder=-1, lw=1)

# --- Entropy ---
ax2.plot(sol_bg['T'], sol_bg['SoverSi'], lw=1, c='k')

ax1.set_title(r'$m_X=' + f'{MXV:g}' + r'\,{\rm GeV},\,m_{Z^\prime}=' + f'{MYV:g}'
              + r'\,{\rm GeV},\,\alpha_X=' + f'{ALPHAXV:.0e}' + r',\,\epsilon='
              + f'{epsX:.1e}' + r',\,\xi_{\rm inf}=1$', fontsize=8)
for ax in (ax1, ax2):
    ax.set_xscale('log')
    ax.set_xlim(1e3, 1e-6)
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(1, 10)*0.1,
                                          numticks=2000))
    ax.xaxis.set_minor_formatter(NullFormatter())
ax1.set_yscale('log')
ax1.set_ylim(0.9e-9, 1.5)
ax1.set_yticks([1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
ax1.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(1, 10)*0.1, numticks=2000))
ax1.yaxis.set_minor_formatter(NullFormatter())
ax1.set_xticklabels([])
ax1.set_ylabel(r'$\rho_j/\rho_{\rm tot}$')
ax1.legend(loc='lower left', handletextpad=0.2)
ax2.set_xlabel(r'$T\,[{\rm GeV}]$')
ax2.set_ylabel(r'$S/S_i$')

plt.show()
plt.close()
