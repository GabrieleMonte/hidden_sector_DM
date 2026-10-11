# 2HDECAY

A program for the Calculation of Electroweak One-Loop Corrections to Higgs Decays in the Two-Higgs-Doublet Model Including State-of-the-Art QCD Corrections

## Program information

**Program** 2HDECAY 1.1.4

**Date** 04.11.2019

**Authors** [Marcel Krause](mailto:marcel.krause@alumni.kit.edu), [M. Margarete Mühlleitner](mailto:margarete.muehlleitner@kit.edu) and [Michael Spira](mailto:michael.spira@psi.ch)

**Manual** [Comp. Phys. Commun. 246 (2020) 106852](https://www.sciencedirect.com/science/article/pii/S0010465519302292), [arXiv:1810.00768 [hep-ph]](https://arxiv.org/abs/1810.00768)

**Citations** When you use this program, please acknowledge the work of our and other groups by citing the following papers:
- The manual for 2HDECAY:
  - M. Krause, M. Muhlleitner, M. Spira, Comp. Phys. Commun. 246 (2020) 106852, arXiv:1810.00768 (hep-ph)
- The manuals for HDECAY:
  - A. Djouadi, J. Kalinowski, M. Spira, Comp. Phys. Commun. 108 (1998) 56, hep-ph/9704448
  - A. Djouadi, J. Kalinowski, M. Muhlleitner, M. Spira, arXiv:1801.09506 (hep-ph)
- The papers on the electroweak correction to the 2HDM decays:
  - M. Krause, R. Lorenz, M. Muhlleitner, R. Santos, H. Ziesche, JHEP 1609 (2016) 143, arXiv:1605.04853 (hep-ph)
  - M. Krause, M. Muhlleitner, R. Santos, H. Ziesche, Phys.Rev. D95 (2017) no.7, 075019, arXiv:1609.04185 (hep-ph)
  - A. Denner, S. Dittmaier, J.-N. Lang, J. High Energ. Phys. (2018) 2018: 104, arXiv:1808.03466 (hep-ph)
- The publication of LoopTools:
  - T. Hahn, M. Perez-Victoria, Comp. Phys. Commun. 118 (1999) 153-165, hep-ph/9807565

**Abstract** We present the program package 2HDECAY for the calculation of the partial decay widths and branching ratios of the Higgs bosons of a general CP-conserving 2-Higgs doublet model (2HDM). The tool includes the full electroweak one-loop corrections to all two-body onshell Higgs decays in the 2HDM that are not loop-induced. It combines them with the state-of-the-art QCD corrections that are already implemented in the program HDECAY. For the renormalization of the electroweak sector an on-shell scheme is implemented for most of the renormalization parameters. Exceptions are the soft-![](https://latex.codecogs.com/gif.latex?%5Cmathbb%7BZ%7D_2 "\mathbb{Z}_2")-breaking squared mass scale ![](https://latex.codecogs.com/gif.latex?m_%7B12%7D%5E2 "m_{12}^2"), where an ![](https://latex.codecogs.com/gif.latex?%5Coverline%7B%5Ctext%7BMS%7D%7D "\overline{\text{MS}}") condition is applied, as well as the 2HDM mixing angles ![](https://latex.codecogs.com/gif.latex?%5Calpha "\alpha") and ![](https://latex.codecogs.com/gif.latex?%5Cbeta "\beta"), for which several distinct renormalization schemes are implemented. The tool 2HDECAY can be used for phenomenological analyses of the branching ratios of Higgs decays in the 2HDM. Furthermore, the separate output of the electroweak contributions to the tree-level partial decay widths for several different renormalization schemes allows for an efficient analysis of the impact of the electroweak corrections and the remaining theoretical error due to missing higher-order corrections. The latest version of the program package 2HDECAY can be downloaded from the URL https://github.com/marcel-krause/2HDECAY.

**Changelog** For a documentation about the changes made in 2HDECAY, check the [Changelog.md](Changelog.md) file.

**Copyright** Copyright (C) 2018-2019, Marcel Krause, Milada Margarete Mühlleitner and Michael Spira

**License** GNU General Public License (GNU GPL-3.0-or-later). 2HDECAY is released under GNU General Public License (GNU GPL-3.0-or-later). This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You have received a copy ([LICENSE.md](LICENSE.md)) of the GNU General Public License along with this program.

**Contact** For feedback, complaints and bug reports, please send an e-mail to <marcel.krause@alumni.kit.edu>, <margarete.muehlleitner@kit.edu> and <michael.spira@psi.ch>.

## Getting Started

For more detailed information about installing and using the program, please have a look at the [full manual](Documentation/manual.pdf "full manual").

### System requirements

Supported operating systems:
- Windows 7 and Windows 10 (with [Cygwin](https://www.cygwin.com/ "Cygwin") installed)
- Linux
- macOS

The following components must be installed on your system to compile and run 2HDECAY:
- GNU gcc (tested with versions 6.4.0 and 7.3.1)
- GNU g++
- GNU gfortran
- Python 2 or 3 (tested with versions 2.7.14 and 3.5.0)
- find
- cURL

If you install 2HDECAY on Windows, make sure that you install [Cygwin](https://www.cygwin.com/ "Cygwin") first, including the following components:
- GNU gcc
- GNU g++
- GNU gfortran
- find
- cURL

### Installing

For an easy installation, we recommend using the automatic installer of 2HDECAY, which guides you through the installation. Download the latest version of 2HDECAY from https://github.com/marcel-krause/2HDECAY. Open a shell, navigate to the 2HDECAY root folder and execute the following command:
```
python setup.py
```

The installer asks you if you want to download and install LoopTools, which is required for 2HDECAY to run. If you do not have LoopTools installed on your system already type y. The installer downloads the LoopTools version specified in the file [Config.py](Config.py); if you prefer another version, change the corresponding entry in that file first before executing the setup.
If you already have LoopTools on your system, you can choose not to download and install LoopTools. In that case, open the file [Config.py](Config.py) and set
```
useRelativeLoopToolsPath = False
```
and change the variable pathLoopTools to the *absolute* path to the LoopTools main folder (containing the bin, include and library subfolders). Additionally, modify the variables pathLoopToolsLibs and pathLoopToolsExecs to the relative paths to the library and bin subfolders with respect to the LoopTools main directory.

The installer further asks you whether the makefile and electroweakCorrections.F90 file shall be created and whether it should make the program. Type y for all these questions in order to compile 2HDECAY.

### Using 2HDECAY

Using 2HDECAY is simple. Save all input files that you want to use in the 'Input' subdirectory. You can put as many input files in the directory as you wish and the input files can have arbitrary file names. Make sure that the input files have the exact input format required by 2HDECAY, as described in the [full manual](Documentation/manual.pdf "full manual").

To start 2HDECAY, execute the following command in the main folder of 2HDECAY:
```
python 2HDECAY.py
```
The program will iterate over all input files and calculate all higher-order corrections according to the parameter and renormalization scheme choices as given in your input files. The resulting output files are saved with the same file name as the input files in the 'Results' subdirectory. For each input file, two output files are generated: one file containing the branching ratios with and without the electroweak corrections, indicated by a filename suffix '_BR', as well as one file containing the electroweak partial decay widths at tree-level and one-loop order, indicated by the filename suffix '_EW'.