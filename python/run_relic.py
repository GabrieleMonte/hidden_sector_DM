import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
import scipy as sp
from scipy import special as scisp
from scipy.integrate import solve_ivp
from scipy.interpolate import InterpolatedUnivariateSpline,CubicSpline
import scipy.integrate
import scipy.optimize
import multiprocessing
from tqdm import trange
from scipy.optimize import curve_fit
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interp1d
from matplotlib.ticker import LogLocator, NullFormatter
import sys
from pathlib import Path
# add main_directory/ to the import path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from hidden_sector_DM.HiddenSectorDM import Cosmology, VectorPortal, BoltzmannSolver

plt.style.use('/home/gab/Desktop/PyCharm_env/mine.mplstyle')
plt.rcParams["axes.axisbelow"] = False



# 1) Initialise the SM background (point gstar_dir to wherever std.tab lives)
cosmo = Cosmology(gstar_choice="standard", gstarpath="./")

# 2) Build the hidden-sector model
model = VectorPortal(
    mX=5e1,       # 50 GeV DM
    mY=2e1,       # 20 GeV mediator
    gX=2, gY=3,
    alphaX=1e-4,
    include_antiparticlesX=True,
    include_antiparticlesY=False,
)

# 3) Create the solver and run
solver = BoltzmannSolver(cosmo, model)
sol0 = solver.solve_boltzmann_chempot_2phase(return_bg_ICs=True, verbose=True,Gamma_switch_threshold=1e6, n_points=500)

epsX = solver.find_epsilon_DMRelic(sol0['bg_ICs'])

# 4) Run background evolution for a particular epsilon
sol_bg = solver.solve_background(epsX=epsX, bg_ICs=sol0['bg_ICs'],
                          stop_on_rhoY=False, verbose=False)

cs=plt.cm.tab20([0,1,2,3,4,5,6,7,8])
fig1=plt.figure(figsize=(4.5, 4.5))
#fig1=plt.figure(figsize=(7.5, 7.5/1.61803398875*0.5))
gs1 = GridSpec(21,23,figure=fig1)
ax1 = fig1.add_subplot(gs1[:13,:])
ax2 = fig1.add_subplot(gs1[13:,:])
ax1.set_title(r'$m_{\chi}=50\,{\rm GeV},\,m_{Z^\prime}=20\,{\rm GeV},\, \alpha_X=0.2,\,\epsilon=10^{-12},\,\xi_{\rm inf}=1$',fontsize=9)

ax1.plot(sol_bg['T'],sol_bg['rhoR']/sol_bg['rhoTot'],lw=1,c=cs[4])
ax1.plot(sol_bg['T'],sol_bg['rhoY']/sol_bg['rhoTot'],lw=1,c=cs[0])
ax1.plot(sol_bg['T'],sol_bg['rhoX']/sol_bg['rhoTot'],lw=1,c=cs[2])

ax1.set_ylim(.9e-9,1.5)
ax1.set_yscale("log")
ax1.set_xscale("log")
ax2.set_xscale("log")
ax1.axvline(1e-2,ls='--',c='tab:gray',zorder=-1,lw=1)
ax2.set_xlabel(r"$T\,[\rm GeV]$")
ax1.set_ylabel(r"$\rho_j/\rho_{\rm tot}$")
ax1.invert_xaxis()
ax2.invert_xaxis()
ax1.set_xlim(1e3,1e-6)
ax2.set_xlim(1e3,1e-6)
ax1.set_xticklabels([])

ax2.plot(sol_bg['T'],sol_bg['SoverSi'],lw=1,c='k')
ax2.set_ylabel(r'$S/S_i$')
ax1.legend()
ax1.set_yticks([1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1])
ax1.yaxis.set_minor_locator(LogLocator(base=10.0, subs=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,.9], numticks=2000))
ax1.yaxis.set_minor_formatter(NullFormatter())
#ax2.set_yticks([1,1e1,1e2,1e3,1e4])
ax2.xaxis.set_minor_locator(LogLocator(base=10.0, subs=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,.9], numticks=2000))
ax2.xaxis.set_minor_formatter(NullFormatter())
ax1.xaxis.set_minor_locator(LogLocator(base=10.0, subs=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,.9], numticks=2000))
ax1.xaxis.set_minor_formatter(NullFormatter())

plt.show()
plt.close()
