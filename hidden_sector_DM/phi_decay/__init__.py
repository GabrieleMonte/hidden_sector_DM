"""
phi_decay — scalar mediator decay widths.

Provides Gamma(phi -> SM) / sin^2(eps) via three mass regimes:
    m < 2 GeV  : Winkler dispersive analysis  (scalar_portal)
    2-5 GeV    : perturbative QCD Nf=4        (scalar_portal)
    >= 5 GeV   : NNLO QCD + NLO EW            (HDECAY)
"""

import subprocess
from pathlib import Path

_hdecay_ready = False


def ensure_hdecay_ready():
    """Check for HDECAY binary; run setup script if missing."""
    global _hdecay_ready
    if _hdecay_ready:
        return

    here = Path(__file__).resolve().parent
    exe = here / "2HDECAY" / "HDECAY" / "run"

    if exe.exists():
        _hdecay_ready = True
        return

    script = here / "setup_phi_decay.sh"
    if not script.exists():
        raise FileNotFoundError(
            f"Setup script not found: {script}\n"
            "Cannot auto-install HDECAY.")

    print("=" * 60)
    print(" HDECAY not found -- running first-time setup...")
    print(" (requires gfortran and git)")
    print("=" * 60)

    result = subprocess.run(
        ["bash", str(script)],
        cwd=str(here),
        timeout=300,
    )

    if result.returncode != 0 or not exe.exists():
        raise RuntimeError(
            "HDECAY setup failed. Please run manually:\n"
            f"  cd {here} && bash setup_phi_decay.sh")

    _hdecay_ready = True
    print("HDECAY setup complete.")


from .phi_decay_hybrid import (
    phi_total_width_normalised,
    phi_total_width_normalized,
    phi_total_width,
    phi_lifetime,
    phi_partial_widths,
    phi_branching_ratios,
    width_to_hh,
)
