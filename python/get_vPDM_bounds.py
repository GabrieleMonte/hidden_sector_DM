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
ROOT = Path.cwd().parent
sys.path.insert(0, str(ROOT))
from source.SecludedDM import Cosmology, VectorPortal, BoltzmannSolver, find_mXstar

plt.style.use('/home/gab/Desktop/PyCharm_env/mine.mplstyle')
plt.rcParams["axes.axisbelow"] = False

alphaXv=1e-3
rv=1.5
cosmo = Cosmology(gstar_choice="standard", gstarpath="./")
model_factory = lambda mX, mY: VectorPortal(mX=mX, mY=mY, gX=2, gY=3, alphaX=alphaXv,
                                             include_antiparticlesX=True, include_antiparticlesY=False)
res = find_mXstar(cosmo, model_factory, r=rv,verbose=True,verbose_solver=True)

print(res['mXstar'], res['epsX_min'])
mXs=np.linspace(res['mXstar']*(1+1e-2),10*res['mXstar'],200)
mYs=mXs/rv
epsXs=np.zeros(len(mXs))
epsXs_bbn=np.zeros(len(mXs))
for i in trange(len(mXs)):
    model = VectorPortal(
        mX=mXs[i],       
        mY=mYs[i],       
        gX=2, gY=3,
        alphaX=alphaXv,
        include_antiparticlesX=True,
        include_antiparticlesY=False,
    )
    solver = BoltzmannSolver(cosmo, model)
    sol0 = solver.solve_boltzmann_chempot_QSSA(return_bg_ICs=True, verbose=False, cannibal_switch_full=1,
                                                 convergence_threshold=2e-2)
    epsXs[i]= solver.find_epsilon_DMRelic(sol0['bg_ICs'],log10eps_range=(-18,-8),root_rtol= 1e-4)
    epsXs_bbn[i]=solver.find_epsilon_BBN_hadronic(sol0['YY_final'],log10eps_range=(-20,-9),root_rtol=1e-3)
    del model, solver
data = np.column_stack([mXs, mYs, epsXs, epsXs_bbn])
header = (f"alphaX={alphaXv:.1e}  r={rv}  mXstar={res['mXstar']:.6e}  "
          f"epsX_min={res['epsX_min']:.6e}  SRatio_plateau={res['SRatio_plateau']:.4f}\n"
          f"mX [GeV]    mY [GeV]    epsX_relic    epsX_bbn_Neff")
np.savetxt(f"../output/vP/relic_aX{alphaXv:.0e}_r{rv}.dat", data, header=header, fmt='%.8e')