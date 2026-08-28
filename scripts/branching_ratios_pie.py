"""
Pie charts of mediator branching ratios for all portal models.

Reproduces Fig. 1 of arXiv:1912.08821 for comparison.
Light quarks (u, d, s, c) are grouped as q-qbar.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('/home/gab/Desktop/PyCharm_env/mine.mplstyle')

from hidden_sector_DM.model import VectorPortal, BLPortal, LiLjPortal, BaryonPortal, HiggsPortal

# ── Parameters ──────────────────────────────────────────────────────
mY = 50.0   # mediator mass [GeV]

models = {
    r'Hypercharge Portal': VectorPortal(
        mX=100.0, mY=mY, alphaX=1e-2,
        Delta1=1e-2, Delta2=1e-2, Delta3=1e-2),
    r'$B-L$ Portal': BLPortal(
        mX=100.0, mY=mY, alphaX=1e-2,
        Delta1=1e-2, Delta2=1e-2, Delta3=1e-2, g_BL=1e-3),
    r'$L_i - L_j$ Portal': LiLjPortal(
        mX=100.0, mY=mY, alphaX=1e-2,
        Delta1=1e-2, Delta2=1e-2, Delta3=1e-2,
        flavor='e-mu', g_LiLj=1e-3),
    r'Baryon Portal': BaryonPortal(
        mX=100.0, mY=mY, alphaX=1e-2,
        Delta1=1e-2, Delta2=1e-2, Delta3=1e-2, g_B=1e-3),
    r'Higgs Portal': HiggsPortal(
        mX=100.0, mY=mY, lam=1e-2),
}

# ── Channel grouping ────────────────────────────────────────────────

def group_channels(br):
    """Group light quarks -> qq, neutrinos -> nunu, e+mu -> e/mu."""
    out = {}
    qq = 0.0
    nunu = 0.0
    for ch, v in br.items():
        if ch in ('uu', 'dd', 'ss', 'cc'):
            qq += v
        elif ch in ('nunu_e', 'nunu_mu', 'nunu_tau'):
            nunu += v
        elif ch in ('ee', 'mumu'):
            out.setdefault('e/mu', 0.0)
            out['e/mu'] += v
        else:
            out[ch] = v
    if qq > 0:
        out['qq'] = qq
    if nunu > 0:
        out['nunu'] = nunu
    return out

# ── Labels and colors ───────────────────────────────────────────────

channel_labels = {
    'e/mu': r'$e/\mu$',
    'tautau': r'$\tau^+\tau^-$',
    'bb': r'$b\bar{b}$', 'tt': r'$t\bar{t}$',
    'qq': r'$q\bar{q}$',
    'nunu': r'$\nu\bar{\nu}$',
    'gg': r'$gg$', 'gamgam': r'$\gamma\gamma$', 'Zgam': r'$Z\gamma$',
    'WW': r'$WW$', 'ZZ': r'$ZZ$', 'hh': r'$hh$',
}

channel_colors = {
    'e/mu':     '#f58231',  # orange
    'tautau':   '#e6194b',  # red
    'nunu':     '#42d4f4',  # cyan
    'qq':       '#ffe119',  # yellow
    'bb':       '#000075',  # navy
    'tt':       '#800000',  # maroon
    'gg':       '#469990',  # teal
    'gamgam':   '#dcbeff',  # lavender
    'Zgam':     '#9A6324',  # brown
    'WW':       '#fabed4',  # pink
    'ZZ':       '#ffd8b1',  # apricot
    'hh':       '#aaffc3',  # mint
    'other':    '#e0e0e0',  # light grey
}

# ── Compute branching ratios ────────────────────────────────────────

all_brs = {}
for name, model in models.items():
    raw = model.branching_ratios_to_SM()
    all_brs[name] = group_channels(raw)

# ── Plot ────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 5, figsize=(20, 4.5))

for ax, (name, _) in zip(axes, models.items()):
    br = all_brs[name]
    # Keep channels >= 1%, group rest into "other"
    filtered = {ch: val for ch, val in br.items() if val >= 0.01}
    other = sum(val for val in br.values() if 0 < val < 0.01)
    if other > 0:
        filtered['other'] = other

    total = sum(filtered.values())
    sorted_items = sorted(filtered.items(), key=lambda x: -x[1])
    labels = [channel_labels.get(ch, ch) for ch, _ in sorted_items]
    sizes = [v / total for _, v in sorted_items]
    colors = [channel_colors.get(ch, '#e0e0e0') for ch, _ in sorted_items]

    wedges, texts, autotexts = ax.pie(
        sizes, labels=labels, colors=colors,
        autopct=lambda p: r'{:.1f}\%'.format(p) if p >= 1 else '',
        startangle=90, normalize=True,
        pctdistance=0.72,
        textprops={'fontsize': 7},
    )
    for at in autotexts:
        at.set_fontsize(5.5)
    ax.set_title(name, fontsize=10, pad=8)

fig.suptitle(
    r'Mediator branching ratios, $m_Y = {:g}$ GeV'.format(mY),
    fontsize=12, y=1.0)
plt.tight_layout()

outdir = os.path.join(os.path.dirname(__file__), '..', 'figs')
os.makedirs(outdir, exist_ok=True)
outpath = os.path.join(outdir, 'branching_ratios_pie.pdf')
plt.savefig(outpath)
print(f'Saved to {outpath}')
