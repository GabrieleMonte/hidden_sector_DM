import numpy as np
from tqdm import trange
import sys
from pathlib import Path

ROOT = Path.cwd().parent
sys.path.insert(0, str(ROOT))
from hidden_sector_DM.HiddenSectorDM import Cosmology, HiggsPortal, BoltzmannSolver, find_mXstar

cosmo = Cosmology(gstar_choice="standard", gstarpath="./")

cases = [
    (0.05, 5.0),
    (0.05, 1.5),
    (0.1,  5.0),
    (0.1,  1.5),
]

for lam, r in cases:
    print(f"\n=== HiggsPortal: lam={lam}, r={r} ===")
    model_factory = lambda mX, mY, l=lam: HiggsPortal(
        mX=mX, mY=mY, gX=2, gY=1, lam=l,
        include_antiparticlesX=True, include_antiparticlesY=False)

    res = find_mXstar(cosmo, model_factory, r=r, verbose=True, verbose_solver=True)
    print(f"mXstar={res['mXstar']:.6e}, epsX_min={res['epsX_min']:.6e}")

    mXs = np.linspace(res['mXstar'] * (1 + 1e-1), 10 * res['mXstar'], 100)
    mYs = mXs / r
    epsXs = np.zeros(len(mXs))
    epsXs_bbn = np.zeros(len(mXs))

    for i in trange(len(mXs)):
        model = HiggsPortal(
            mX=mXs[i], mY=mYs[i], gX=2, gY=1, lam=lam,
            include_antiparticlesX=True, include_antiparticlesY=False)
        solver = BoltzmannSolver(cosmo, model)
        sol0 = solver.solve_boltzmann_chempot_3phase(
            return_bg_ICs=True, verbose=False,
            cannibal_switch_full=1, convergence_threshold=2e-2)
        epsXs[i] = solver.find_epsilon_DMRelic(
            sol0['bg_ICs'], log10eps_range=(-18, -2), root_rtol=1e-4)
        epsXs_bbn[i] = solver.find_epsilon_BBN_hadronic(
            sol0['YY_final'], log10eps_range=(-20, -2), root_rtol=1e-3)
        del model, solver

    data = np.column_stack([mXs, mYs, epsXs, epsXs_bbn])
    header = (f"lam={lam:.2e}  r={r}  mXstar={res['mXstar']:.6e}  "
              f"epsX_min={res['epsX_min']:.6e}  SRatio_plateau={res['SRatio_plateau']:.4f}\n"
              f"mX [GeV]    mY [GeV]    epsX_relic    epsX_bbn")
    np.savetxt(f"../output/hP/relic_lam{lam:.0e}_r{r}.dat",
               data, header=header, fmt='%.8e')
