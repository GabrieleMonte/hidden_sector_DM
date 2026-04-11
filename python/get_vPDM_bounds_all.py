import numpy as np
from tqdm import trange
import sys
from pathlib import Path

ROOT = Path.cwd().parent
sys.path.insert(0, str(ROOT))
from hidden_sector_DM.HiddenSectorDM import Cosmology, VectorPortal, BoltzmannSolver, find_mXstar

cosmo = Cosmology(gstar_choice="standard", gstarpath="./")

cases = [
    (1e-4, 1.5, 4),
    (1e-4, 2.0, 4),
    (1e-4, 5.0, 4),
    (1e-3, 1.5, 4),
    (1e-3, 2.0, 4),
    (1e-3, 5.0, 4),
    (1e-2, 1.5, 5),
    (1e-2, 2.0, 5),
    (1e-2, 5.0, 5),
    (5e-4, 1.5, 5),
    (5e-4, 2.0, 5),
    (5e-4, 5.0, 5),
    (5e-3, 1.5, 5),
    (5e-3, 2.0, 5),
    (5e-3, 5.0, 5),
    (5e-2, 1.5, 5),
    (5e-2, 2.0, 5),
    (5e-2, 5.0, 5),
]

for alphaX, r, log10mXmax in cases:
    print(f"\n=== VectorPortal: alphaX={alphaX:.0e}, r={r} ===")
    model_factory = lambda mX, mY, aX=alphaX: VectorPortal(
        mX=mX, mY=mY, gX=2, gY=3, alphaX=aX,
        include_antiparticlesX=True, include_antiparticlesY=False)

    res = find_mXstar(cosmo, model_factory, r=r, verbose=True, verbose_solver=False)
    print(f"mXstar={res['mXstar']:.6e}, epsX_min={res['epsX_min']:.6e}")

    MZ = 91.1876  # GeV
    mX_res = r * MZ  # mX where mZ' hits the Z pole
    mX_lo = res['mXstar'] * (1 + 5e-2)
    mX_hi = 10**log10mXmax

    # Near mXstar: log-spaced for good resolution at the lower limit
    mXs_near = np.logspace(np.log10(mX_lo),
                           np.log10(min(10 * res['mXstar'], mX_hi)), 150)
    # Around Z-pole resonance (mZ' = MZ): linear cluster
    if mX_lo < mX_res * 1.3 and mX_res * 0.7 < mX_hi:
        res_lo = max(mX_res * 0.8, mX_lo)
        res_hi = min(mX_res * 1.2, mX_hi)
        mXs_res = np.linspace(res_lo, res_hi, 50)
    else:
        mXs_res = np.array([])
    # High-mass tail: log-spaced
    mXs_high = np.logspace(np.log10(max(10.01 * res['mXstar'], mX_lo)),
                           np.log10(mX_hi), 50)
    mXs = np.unique(np.concatenate([mXs_near, mXs_res, mXs_high]))
    mYs = mXs / r
    epsXs = np.zeros(len(mXs))
    epsXs_bbn = np.zeros(len(mXs))

    for i in trange(len(mXs)):
        model = VectorPortal(
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
    np.savetxt(f"../output/vP/relic_aX{alphaX:.0e}_r{r}.dat",
               data, header=header, fmt='%.8e')
