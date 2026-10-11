#!/bin/bash
# setup_phi_decay.sh
# ==================
# Sets up HDECAY and (optionally) scalar_portal for phi_decay_hybrid.py.
# Installs everything relative to this script's own directory.
#
# Requirements: gfortran, git, Python 3.8+ with numpy

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "============================================="
echo " Setting up phi_decay environment"
echo " Install directory: $SCRIPT_DIR"
echo "============================================="

# ── (1) Check prerequisites ──
echo ""
echo "[1/3] Checking prerequisites..."

if ! command -v gfortran &> /dev/null; then
    echo "ERROR: gfortran not found. Install with:"
    echo "  sudo apt install gfortran      (Debian/Ubuntu)"
    echo "  brew install gcc                (macOS)"
    exit 1
fi
echo "  gfortran: OK"

if ! command -v git &> /dev/null; then
    echo "ERROR: git not found."
    exit 1
fi
echo "  git: OK"

# ── (2) Clone and compile HDECAY ──
echo ""
echo "[2/3] Setting up HDECAY..."

if [ -f "2HDECAY/HDECAY/run" ]; then
    echo "  HDECAY binary already exists, skipping."
else
    if [ ! -d "2HDECAY" ]; then
        echo "  Cloning 2HDECAY repository..."
        git clone --depth 1 https://github.com/marcel-krause/2HDECAY.git
    fi

    cd 2HDECAY/HDECAY

    # Patch: allow pure SM mode (comment out forced 2HDM activation)
    echo "  Patching hdecay.f for standalone SM mode..."
    cp hdecay.f hdecay.f.bak
    python3 -c "
with open('hdecay.f', 'r') as f:
    t = f.read()
old = '''      if(ielw2hdm.eq.0.and.i2hdm.eq.0) then
         print*,''
         print*,\'You chose to calculate the EW corrections to the 2HDM d
     .ecay widths but did not turn on the flag for the 2HDM. This is don
     .e now.\'
         i2hdm = 1
      endif'''
new = '''c      if(ielw2hdm.eq.0.and.i2hdm.eq.0) then
c         i2hdm = 1
c      endif'''
assert old in t, 'Patch target not found in hdecay.f'
with open('hdecay.f', 'w') as f:
    f.write(t.replace(old, new))
print('  Patch applied.')
"

    # Create empty auxiliary files expected by 2HDECAY-modified reader
    touch alphaandbeta.dat fermionmasses.dat

    # Compile
    echo "  Compiling HDECAY..."
    FFLAGS="-fallow-argument-mismatch -O2" make hdecay

    if [ ! -f "run" ]; then
        echo "ERROR: Compilation failed -- 'run' not produced."
        exit 1
    fi
    echo "  HDECAY compiled successfully."

    # Prepare a clean hdecay.in for SM mode
    sed -i 's/^SLHAIN   = .*/SLHAIN   = 0/' hdecay.in
    sed -i 's/^SLHAOUT  = .*/SLHAOUT  = 0/' hdecay.in
    sed -i 's/^COUPVAR  = .*/COUPVAR  = 0/' hdecay.in
    sed -i 's/^HIGGS    = .*/HIGGS    = 0/' hdecay.in
    sed -i 's/^2HDM     = .*/2HDM     = 0/' hdecay.in
    sed -i 's/^OMIT ELW = .*/OMIT ELW = 0/' hdecay.in
    sed -i 's/^OMIT ELW2= .*/OMIT ELW2= 0/' hdecay.in

    cd "$SCRIPT_DIR"
fi

# ── (3) scalar_portal (optional, for m < 5 GeV) ──
echo ""
echo "[3/3] Setting up scalar_portal (optional, for m < 5 GeV)..."

if [ -d "scalar_portal" ]; then
    echo "  scalar_portal already present, skipping."
elif python3 -c "from scalar_portal import Model" 2>/dev/null; then
    echo "  scalar_portal already importable, skipping."
else
    echo "  Attempting to clone scalar_portal..."
    git clone --depth 1 https://github.com/jlp-lu/scalar_portal.git 2>/dev/null && {
        echo "  scalar_portal cloned."
    } || {
        echo "  WARNING: Could not clone scalar_portal."
        echo "  The code will still work for m_phi >= 5 GeV."
    }
fi

echo ""
echo "============================================="
echo " Setup complete."
echo "============================================="
