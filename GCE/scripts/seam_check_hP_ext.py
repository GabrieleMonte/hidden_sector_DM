"""Continuity of the photon-spectrum cache across the mX = 100 GeV seam.

The Higgs-portal extension (make_channel_spectra.extension_grid) continues the main
grid's log-mX spacing past MX_LIM, so the first extension column sits one step of
du = 0.0067 dex above the main grid's last.  Nothing in the driver enforces that the two
sides were generated consistently; this measures it.

Stages:
  1. Read N_gamma per decay along the main grid's last column and the extension's first,
     summing each cached file's own `gamma` histogram.
  2. Interpolate the main column in log mY onto the extension's rows, over the mY range
     the two share.  N_gamma varies steeply with mY, so a nearest-neighbour read would
     show that gradient as a seam discontinuity that is not there.
  3. Report the ratio per channel.  A 1.55% step in mX should move it by far less than
     the Monte-Carlo scatter of an N_EVENTS-event file.

Output: stdout only.
"""

import os
import re
from pathlib import Path

import numpy as np

# ---- configuration --------------------------------------------------------

ROOT = Path(__file__).resolve().parents[2]
CACHE_DIR = Path(os.environ.get('CHANNEL_SPECTRA_DIR', ROOT/'GCE'/'channel_spectra'))
MAIN_MX = '100'               # the main grid's last column, as the filenames spell it
EXT_MX = '101.55'             # the extension's first column
CHANNELS = ('bb', 'cc', 'tautau', 'gg', 'WW')
N_EVENTS = 400_000
SEED = 12345

_NAME = re.compile(rf'^(?P<ch>\w+)_mX(?P<mX>[\d.]+)_mY(?P<mY>[\d.]+)_'
                   rf'N{N_EVENTS}_s{SEED}\.npz$')


def column(channel, mx_tag):
    """(mY [GeV], N_gamma per decay) along one mX column, sorted by mY."""
    rows = []
    for path in CACHE_DIR.glob(f'{channel}_mX{mx_tag}_mY*_N{N_EVENTS}_s{SEED}.npz'):
        name = _NAME.match(path.name)
        if name is None or name.group('mX') != mx_tag:
            continue                      # 'mX100' must not swallow 'mX100.5'
        with np.load(path) as cached:
            rows.append((float(name.group('mY')), float(cached['gamma'].sum())))

    if not rows:
        return np.empty(0), np.empty(0)
    rows.sort()
    return np.array([r[0] for r in rows]), np.array([r[1] for r in rows])


# ---- run ------------------------------------------------------------------

mx_step = 100*(float(EXT_MX)/float(MAIN_MX) - 1)
print(f'cache: {CACHE_DIR}')
print(f'seam: mX = {MAIN_MX} (main) -> {EXT_MX} (extension), '
      f'a {mx_step:.2f}% step in mX\n')
print(f"{'channel':>8} {'rows':>5} {'overlap mY':>18} {'median':>8} {'min':>8} {'max':>8}")

worst = 0.0
for channel in CHANNELS:
    mY_main, n_main = column(channel, MAIN_MX)
    mY_ext, n_ext = column(channel, EXT_MX)
    if mY_main.size < 2 or mY_ext.size == 0:
        print(f'{channel:>8} {"-":>5}  no overlap '
              f'(main {mY_main.size}, extension {mY_ext.size})')
        continue

    shared = (mY_ext >= mY_main.min()) & (mY_ext <= mY_main.max())
    if shared.sum() == 0:
        print(f'{channel:>8} {"-":>5}  extension rows lie outside the main column')
        continue

    reference = np.interp(np.log(mY_ext[shared]), np.log(mY_main), n_main)
    ratio = n_ext[shared]/reference
    worst = max(worst, float(np.abs(ratio - 1).max()))

    print(f'{channel:>8} {shared.sum():>5} {mY_ext[shared].min():8.2f}-'
          f'{mY_ext[shared].max():<9.2f}{np.median(ratio):8.4f} '
          f'{ratio.min():8.4f} {ratio.max():8.4f}')

print(f'\nlargest deviation from 1 across all channels: {100*worst:.2f}%')
