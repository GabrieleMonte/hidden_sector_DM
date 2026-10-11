"""
hdecay_interface.py — Python wrapper for the HDECAY Fortran code.

Provides SM Higgs decay widths and branching ratios by calling the HDECAY
executable (Djouadi, Kalinowski, Muhlleitner, Spira) as a subprocess.

For the scalar portal dark matter application:
    Gamma_phi(m_phi) = sin^2(eps) * Gamma_SM(m_H = m_phi) + Gamma(phi->hh)

where Gamma_SM is computed by HDECAY and Gamma(phi->hh) is model-dependent.

References:
    [1] A. Djouadi, J. Kalinowski, M. Spira, CPC 108 (1998) 56
    [2] A. Djouadi et al., CPC 238 (2019) 214 (HDECAY: Twenty++ years after)
"""

import os
import re
import shutil
import subprocess
import tempfile
import numpy as np
from pathlib import Path

# ── Location of compiled HDECAY ──
_HDECAY_DIR = Path(__file__).resolve().parent / "2HDECAY" / "HDECAY"
_HDECAY_EXE = _HDECAY_DIR / "run"

# ── Default SM parameters (YR4 / PDG 2024) ──
_DEFAULT_PARAMS = dict(
    alpha_s  = 1.18000e-01,
    m_u      = 9.50000e-02,     # MSbar at 2 GeV
    m_c      = 0.98600e+00,     # MSbar at 3 GeV
    m_b      = 4.18000e+00,     # MSbar at m_b
    m_t      = 1.73200e+02,     # pole mass
    m_tau    = 1.77682e+00,
    m_mu     = 1.056583715e-01,
    inv_alpha = 1.37036e+02,
    GF       = 1.1663787e-05,
    M_Z      = 9.11876e+01,
    M_W      = 8.0385e+01,
)


def _parse_br_files(workdir):
    """Parse br.sm1 and br.sm2 output files from HDECAY SM mode.

    br.sm1 columns: MHSM, BB, TAU TAU, MU MU, SS, CC, TT
    br.sm2 columns: MHSM, GG, GAM GAM, Z GAM, WW, ZZ, WIDTH

    Returns list of dicts, one per mass point.
    """
    results = []

    sm1_path = os.path.join(workdir, "br.sm1")
    sm2_path = os.path.join(workdir, "br.sm2")

    # Parse br.sm1
    sm1_data = []
    with open(sm1_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('M') or line.startswith('_'):
                continue
            vals = line.split()
            if len(vals) >= 7:
                sm1_data.append([float(v) for v in vals])

    # Parse br.sm2
    sm2_data = []
    with open(sm2_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('M') or line.startswith('_'):
                continue
            vals = line.split()
            if len(vals) >= 7:
                sm2_data.append([float(v) for v in vals])

    for i in range(len(sm1_data)):
        r1 = sm1_data[i]
        r2 = sm2_data[i]
        m = r1[0]
        total_width = r2[6]  # GeV

        results.append({
            'mass':    m,
            'bb':      r1[1],
            'tautau':  r1[2],
            'mumu':    r1[3],
            'ss':      r1[4],
            'cc':      r1[5],
            'tt':      r1[6],
            'gg':      r2[1],
            'gamgam':  r2[2],
            'Zgam':    r2[3],
            'WW':      r2[4],
            'ZZ':      r2[5],
            'total_width': total_width,   # GeV
        })

    return results


def _prepare_input(workdir, mass_beg, mass_end, nma, params, hdecay_dir):
    """Prepare hdecay.in by copying the original and modifying key flags."""

    # Read original input file (has all 2HDM data blocks intact)
    orig = os.path.join(hdecay_dir, "hdecay.in")
    with open(orig, 'r') as f:
        text = f.read()

    # SM mode flags
    replacements = {
        r'^SLHAIN\s*=.*':    'SLHAIN   = 0',
        r'^SLHAOUT\s*=.*':   'SLHAOUT  = 0',
        r'^COUPVAR\s*=.*':   'COUPVAR  = 0',
        r'^HIGGS\s*=.*':     'HIGGS    = 0',
        r'^2HDM\s*=.*':      '2HDM     = 0',
        r'^OMIT ELW\s*=.*':  'OMIT ELW = 0',
        r'^OMIT ELW2=.*':    'OMIT ELW2= 0',
        r'^SM4\s*=.*':       'SM4      = 0',
        r'^FERMPHOB\s*=.*':  'FERMPHOB = 0',
        r'^MABEG\s*=.*':     f'MABEG    = {mass_beg:.5e}',
        r'^MAEND\s*=.*':     f'MAEND    = {mass_end:.5e}',
        r'^NMA\s*=.*':       f'NMA      = {nma:d}',
    }

    # SM parameter overrides
    param_map = {
        'alpha_s':  (r'^ALS\(MZ\)\s*=.*',   'ALS(MZ)  = {:.5e}'),
        'm_c':      (r'^MCBAR\(3\)\s*=.*',   'MCBAR(3) = {:.5e}'),
        'm_b':      (r'^MBBAR\(MB\)=.*',     'MBBAR(MB)= {:.5e}'),
        'm_t':      (r'^MT\s*=.*',           'MT       = {:.5e}'),
        'm_tau':    (r'^MTAU\s*=.*',         'MTAU     = {:.5e}'),
        'GF':       (r'^GF\s*=.*',           'GF       = {:.7e}'),
        'M_Z':      (r'^MZ\s*=.*',           'MZ       = {:.5e}'),
        'M_W':      (r'^MW\s*=.*',           'MW       = {:.4e}'),
    }

    if params:
        for key, val in params.items():
            if key in param_map:
                pattern, fmt = param_map[key]
                replacements[pattern] = fmt.format(val)

    for pattern, replacement in replacements.items():
        text = re.sub(pattern, replacement, text, flags=re.MULTILINE)

    out_path = os.path.join(workdir, "hdecay.in")
    with open(out_path, 'w') as f:
        f.write(text)


def run_hdecay(masses, params=None, hdecay_dir=None):
    """Run HDECAY for one or more SM Higgs masses.

    Parameters
    ----------
    masses : float or array-like
        Higgs mass(es) in GeV.
    params : dict, optional
        Override default SM parameters (alpha_s, m_t, etc.)
    hdecay_dir : str or Path, optional
        Path to directory containing the HDECAY executable 'run'.

    Returns
    -------
    list of dict
        Each dict contains branching ratios and total width for one mass.
    """
    if hdecay_dir is None:
        hdecay_dir = _HDECAY_DIR
    hdecay_dir = Path(hdecay_dir)
    exe = hdecay_dir / "run"

    if not exe.exists():
        raise FileNotFoundError(
            f"HDECAY executable not found at {exe}. "
            f"Compile with: cd {hdecay_dir} && "
            "FFLAGS='-fallow-argument-mismatch -O2' make hdecay"
        )

    # Handle single mass or array
    masses = np.atleast_1d(np.asarray(masses, dtype=float))
    mass_beg = float(masses.min())
    mass_end = float(masses.max())
    nma = len(masses) if len(masses) > 1 else 1
    if nma == 1:
        mass_end = mass_beg

    with tempfile.TemporaryDirectory() as tmpdir:
        # Prepare input file (copy original + modify flags)
        _prepare_input(tmpdir, mass_beg, mass_end, nma, params, hdecay_dir)

        # Create dummy auxiliary files needed by 2HDECAY-modified source
        with open(os.path.join(tmpdir, "alphaandbeta.dat"), 'w') as f:
            f.write("         0.00000000000000000E+00\n" * 2)
        with open(os.path.join(tmpdir, "fermionmasses.dat"), 'w') as f:
            pass

        # Copy executable
        exe_dest = os.path.join(tmpdir, "run")
        shutil.copy2(str(exe), exe_dest)
        os.chmod(exe_dest, 0o755)

        # Run HDECAY
        result = subprocess.run(
            ["./run"],
            cwd=tmpdir,
            capture_output=True,
            text=True,
            timeout=120,
        )

        sm1 = os.path.join(tmpdir, "br.sm1")
        if not os.path.exists(sm1):
            raise RuntimeError(
                f"HDECAY did not produce output.\n"
                f"Return code: {result.returncode}\n"
                f"stdout: {result.stdout[:500]}\n"
                f"stderr: {result.stderr[:500]}"
            )

        return _parse_br_files(tmpdir)


def hdecay_total_width(mass, params=None):
    """Get SM Higgs total decay width at a single mass [GeV]."""
    results = run_hdecay(mass, params=params)
    return results[0]['total_width']


def hdecay_branching_ratios(mass, params=None):
    """Get SM Higgs branching ratios at a single mass."""
    results = run_hdecay(mass, params=params)
    return results[0]


# ── Caching for batch efficiency ──
_cache = {}


def hdecay_total_width_cached(mass, params=None):
    """Cached version — avoids re-running HDECAY for repeated calls."""
    key = round(mass, 3)
    if key not in _cache:
        _cache[key] = hdecay_total_width(mass, params=params)
    return _cache[key]


def hdecay_mass_scan(mass_min, mass_max, npoints, params=None):
    """Efficient mass scan using HDECAY's built-in linear scan."""
    masses = np.linspace(mass_min, mass_max, npoints)
    return run_hdecay(masses, params=params)


def build_interpolator(mass_min=10.0, mass_max=1000.0, npoints=200, params=None):
    """Build a fast interpolating function for the SM Higgs total width.

    Runs HDECAY once over a mass grid, then returns a function that
    interpolates log(Gamma) vs log(m) for fast evaluation.
    """
    from scipy.interpolate import CubicSpline

    results = hdecay_mass_scan(mass_min, mass_max, npoints, params=params)

    m_arr = np.array([r['mass'] for r in results])
    w_arr = np.array([r['total_width'] for r in results])

    log_m = np.log(m_arr)
    log_w = np.log(w_arr)

    cs = CubicSpline(log_m, log_w)

    def interpolated_width(mass):
        return float(np.exp(cs(np.log(mass))))

    return interpolated_width
