# Mathematica — Cross-section calculations

Symbolic computation of all 2→2 and 3→2 dark-sector cross sections using
[FeynCalc](https://feyncalc.github.io/).  The results feed into the
analytic formulas in `source/model.py`.

## Notebooks

| Notebook | Process | Model | DM spin |
|----------|---------|-------|---------|
| `XX_YY_hP_and_vP.nb` | X X → Y Y (2→2 annihilation) | Higgs Portal (Majorana X, scalar Y) + Vector Portal (Dirac X, vector Y) | both |
| `YXX_XX_vP.nb` | Y X X → X X (3→2 cannibal) | Vector Portal | Dirac |
| `YXX_XX_hP.nb` | Y X X → X X (3→2 cannibal) | Higgs Portal | Majorana |
| `YYX_YX_vP.nb` | Y Y X → Y X (3→2 cannibal) | Vector Portal | Dirac |
| `YYX_YX_hP.nb` | Y Y X → Y X (3→2 cannibal) | Higgs Portal | Majorana |

Each notebook follows a similar structure: define the Feynman rules, contract
the amplitude using FeynCalc, square and sum/average over spins, integrate
over phase space, and expand in the non-relativistic limit to extract
⟨σv⟩ (2→2) or ⟨σ²v²⟩ (3→2).

## Squared-amplitude term files

For the more complex 3→2 processes, individual |M|² contributions
(diagram × diagram) are exported to `.m` files so the full result can be
assembled and checked term by term.

| Directory | Process | # diagrams |
|-----------|---------|-----------|
| `MsqTerms_YXX_XX_vP/` | Y X X → X X (Vector Portal) | 8 |
| `MsqTerms_YXX_XX_hP/` | Y X X → X X (Higgs Portal) | 16 |
| `MsqTerms_YYX_YX_hP/` | Y Y X → Y X (Higgs Portal) | 6 |

Files are named `diag_i_i.m` (diagonal) and `offdiag_i_j.m` (interference).

## Reference PDFs

- `FeynDiagrams_3to2_processes_vP.pdf` — Feynman diagrams for all 3→2
  processes in the Vector Portal.
- `YYX_to_YX_sigmav2_vP.pdf` — Step-by-step derivation of ⟨σ²v²⟩ for
  Y Y X → Y X in the Vector Portal.
