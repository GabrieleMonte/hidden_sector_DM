import numpy as np
from tqdm import trange
import sys
from pathlib import Path

ROOT = Path.cwd().parent
sys.path.insert(0, str(ROOT))
from hidden_sector_DM.HiddenSectorDM import Cosmology, BLPortal, BoltzmannSolver, find_mXstar

cosmo = Cosmology(gstar_choice="standard", gstarpath="./")

cases = [
    (1e-3, 1.5),
    (1e-3, 5.0),
    (5e-3, 1.5),
    (5e-3, 5.0),
]

for alphaX, r in cases:
    print(f"\n=== BLPortal: alphaX={alphaX:.0e}, r={r} ===")
    model_factory = lambda mX, mY, aX=alphaX: BLPortal(
        mX=mX, mY=mY, gX=2, gY=3, alphaX=aX,
        include_antiparticlesX=True, include_antiparticlesY=False)

    res = find_mXstar(cosmo, model_factory, r=r, verbose=True, verbose_solver=True)
    print(f"mXstar={res['mXstar']:.6e}, epsX_min={res['epsX_min']:.6e}")

    mXs = np.linspace(res['mXstar'] * (1 + 1e-2), 10 * res['mXstar'], 200)
    mYs = mXs / r
    epsXs = np.zeros(len(mXs))
    epsXs_bbn = np.zeros(len(mXs))

    for i in trange(len(mXs)):
        model = BLPortal(
            mX=mXs[i], mY=mYs[i], gX=2, gY=3, alphaX=alphaX,
            include_antiparticlesX=True, include_antiparticlesY=False)
        solver = BoltzmannSolver(cosmo, model)
        sol0 = solver.solve_boltzmann_chempot_3phase(
            return_bg_ICs=True, verbose=False,
            cannibal_switch_full=1, convergence_threshold=2e-2)
        epsXs[i] = solver.find_epsilon_DMRelic(
            sol0['bg_ICs'], log10eps_range=(-18, -8), root_rtol=1e-4)
        epsXs_bbn[i] = solver.find_epsilon_BBN_hadronic(
            sol0['YY_final'], log10eps_range=(-20, -9), root_rtol=1e-3)
        del model, solver

    data = np.column_stack([mXs, mYs, epsXs, epsXs_bbn])
    header = (f"alphaX={alphaX:.1e}  r={r}  mXstar={res['mXstar']:.6e}  "
              f"epsX_min={res['epsX_min']:.6e}  SRatio_plateau={res['SRatio_plateau']:.4f}\n"
              f"mX [GeV]    mY [GeV]    epsX_relic    epsX_bbn")
    np.savetxt(f"../output/bLP/relic_aX{alphaX:.0e}_r{r}.dat",
               data, header=header, fmt='%.8e')
