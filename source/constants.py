"""
Physical constants and SM fermion table.
All quantities in GeV unless stated otherwise.
"""

import numpy as np

# --- Fermion masses ---
Me   = 0.51099895e-3
Mmu  = 105.6583755e-3
Mtau = 1.77686
Mu   = 2.16e-3
Mc   = 1.27
Mt   = 172.69
Md   = 4.67e-3
Ms   = 93.4e-3
Mb   = 4.18

# --- Meson / boson masses ---
Mpip = 139.57039e-3
MKp  = 493.677e-3
MZ   = 91.1876
MW   = 80.379

# --- Planck mass (reduced) ---
Mpl = 2.435e18

# --- Weinberg angle ---
sW2 = 0.23121
sW  = np.sqrt(sW2)
cW  = np.sqrt(1.0 - sW2)
tW  = np.sqrt(sW2 / (1.0 - sW2))
s2W = np.sqrt(sW2 - sW2 * sW2)

# --- Couplings ---
alphaS  = 0.1179                              # alpha_s(mZ)
gW      = np.sqrt(4.0 * np.pi / 127.93) / sW  # weak coupling at mZ
gY      = gW * sW / cW                         # hypercharge coupling at mZ
WZ      = 2.4952                               # Z width  (GeV)
alphaEM = 1.0 / 137.0
eC      = np.sqrt(alphaEM * 4 * np.pi)

# --- Cosmological observables ---
s0_cosmo = 2891.0          # entropy density today  (cm^-3)
rhoc     = 1.05375e-5      # critical density       (GeV h^2 cm^-3)
Och2     = 0.12            # Omega_c h^2
mY_relic = Och2 * rhoc / s0_cosmo   # reference relic yield (GeV)

# --- SM fermion table ---
# key : (mass, Q, YL, YR, Nc)
SM_MYQ = {
    # charged leptons
    "e":   (Me,   -1.0, -0.5, -1.0, 1),
    "mu":  (Mmu,  -1.0, -0.5, -1.0, 1),
    "tau": (Mtau, -1.0, -0.5, -1.0, 1),
    # up-type quarks
    "u": (Mu, 2.0/3, 1.0/6,  2.0/3, 3),
    "c": (Mc, 2.0/3, 1.0/6,  2.0/3, 3),
    "t": (Mt, 2.0/3, 1.0/6,  2.0/3, 3),
    # down-type quarks
    "d": (Md, -1.0/3, 1.0/6, -1.0/3, 3),
    "s": (Ms, -1.0/3, 1.0/6, -1.0/3, 3),
    "b": (Mb, -1.0/3, 1.0/6, -1.0/3, 3),
}

# --- B-L fermion table ---
# key : (mass, (B-L) charge, Nc)
BL_FERMIONS = {
    # quarks: (B-L) = +1/3
    "u": (Mu, 1.0/3, 3), "d": (Md, 1.0/3, 3),
    "s": (Ms, 1.0/3, 3), "c": (Mc, 1.0/3, 3),
    "b": (Mb, 1.0/3, 3), "t": (Mt, 1.0/3, 3),
    # charged leptons: (B-L) = -1
    "e": (Me, -1.0, 1), "mu": (Mmu, -1.0, 1), "tau": (Mtau, -1.0, 1),
    # neutrinos: (B-L) = -1
    "nu_e": (0.0, -1.0, 1), "nu_mu": (0.0, -1.0, 1), "nu_tau": (0.0, -1.0, 1),
}
