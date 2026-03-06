module constants
    implicit none
    save

	! Copyright (C) 2018-2019, Marcel Krause, Milada Margarete Muehlleitner and Michael Spira
	
	! License: GNU General Public License (GNU GPL-3.0-or-later)

	! 2HDECAY is released under GNU General Public License (GNU GPL-3.0-or-later).
	! This program is free software: you can redistribute it and/or modify it under the terms of the
	! GNU General Public License as published by the Free Software Foundation, either version 3 of
	! the License, or any later version.

	! This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
	! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	! See the GNU General Public License for more details.

	! You have received a copy (LICENSE.md) of the GNU General Public License along with this program
	! in the 2HDECAY root directory.

	! Integer for choosing whether the electroweak corrections should be calculated or not
	integer :: omitELCorr

    ! Imaginary unit for convenience
    double complex, parameter :: I = (0D0, 1D0)

    ! Pi and its square
    double precision, parameter :: PI = 4.D0*atan(1.D0)
    double precision, parameter :: PI2 = PI**2

    ! "Photon mass" for the IR divergence regulation (the final amplitude is independent of IRLambda)
    double precision :: IRLambda

    ! Gauge-fixing parameters
    double precision :: GaugeXiW = 1.0D0
    double precision :: GaugeXiZ = 1.0D0
    double precision :: GaugeXiA = 1.0D0

	! vev for the tree-level mass calculation
	double precision :: vevCalc

	! Fermi constant for consistency cross-checks
	double precision :: GFinput

	! UV scale and regularization scale
	double precision :: UVDelta = 0D0
	double precision :: InputScale, OutputScale, EvalScale
	character(len=50) :: OutputScaleReadIn

	! 2HDM type, input parameter type (1: masses and alpha are given; 2: lambda1 to lambda5 are given), renormalization scheme and reference scheme
	integer :: TypeOf2HDM
	integer :: parameterType
	integer :: RenormScheme
	integer :: RefScheme

    ! Standard Model parameters; the values are stored in the input file in the ./Parameters/input folder (1101.0593 [hep-ph], 1503.07589 [hep-ex],  Chin. Phys. C38 (2014) 090001)
	double precision :: alphaAtMZ
    double precision :: Mh0
	double precision :: MW
	double precision :: MZ
	double precision :: ME
	double precision :: MM
	double precision :: ML
	double precision :: MU
	double precision :: MC
	double precision :: MT
	double precision :: MD
	double precision :: MS
	double precision :: MB
	double precision :: EL
	double precision :: SW
	double precision :: CW
	double precision :: CKM11
	double precision :: CKM12
	double precision :: CKM13
	double precision :: CKM21
	double precision :: CKM22
	double precision :: CKM23
	double precision :: CKM31
	double precision :: CKM32
	double precision :: CKM33
	double precision :: CKMC11
	double precision :: CKMC12
	double precision :: CKMC13
	double precision :: CKMC21
	double precision :: CKMC22
	double precision :: CKMC23
	double precision :: CKMC31
	double precision :: CKMC32
	double precision :: CKMC33

    ! Squared Standard Model parameters
    double precision :: Mh02
	double precision :: MW2
	double precision :: MZ2
	double precision :: ME2
	double precision :: MM2
	double precision :: ML2
	double precision :: MU2
	double precision :: MC2
	double precision :: MT2
	double precision :: MD2
	double precision :: MS2
	double precision :: MB2
	double precision :: EL2
	double precision :: SW2
	double precision :: CW2

    ! 2HDM-specific parameters (these are set in the folder ./Parameters/input by reading the respective input files)
    double precision :: MHH
	double precision :: MA0
	double precision :: MHp
	double precision :: alpha
	double precision :: beta
	double precision :: CA
	double precision :: SA
	double precision :: TA
	double precision :: CB
	double precision :: SB
	double precision :: TB
	double precision :: S2A
	double precision :: C2A
	double precision :: S2B
	double precision :: C2B
	double precision :: CAB
	double precision :: SAB
	double precision :: CBA
	double precision :: SBA
    double precision :: Yuk1
	double precision :: Yuk2
	double precision :: Yuk3
    double precision :: Yuk4
	double precision :: Yuk5
	double precision :: Yuk6
	double precision :: Lambda5
	double precision :: m12squared
	double precision :: hdecayLam1
	double precision :: hdecayLam2
	double precision :: hdecayLam3
	double precision :: hdecayLam4
	double precision :: hdecayLam5

    ! Squared 2HDM-specific parameters
    double precision :: MHH2
	double precision :: MA02
	double precision :: MHp2
    double precision :: CA2
	double precision :: SA2
	double precision :: TA2
	double precision :: TB2
	double precision :: SB2
	double precision :: CB2
	double precision :: C2A2
	double precision :: S2A2
	double precision :: C2B2
	double precision :: S2B2
	double precision :: CAB2
	double precision :: SAB2
	double precision :: CBA2
	double precision :: SBA2

contains
    ! This is the three-point function with C0(0,p,p,0,0,m), appearing in calculations for xi != 1, which diverges individually,
	! but cancels overall in the full amplitude; hence, we set it to zero to avoid the artificial divergence
    double complex function C0Mine(a,b,c,d,e,f)
        implicit none
        double precision, intent(in) :: a,b,c,d,e,f
        C0Mine = (0D0,0D0)
    end function C0Mine

	! This is the four-point function with zero arguments, appearing in calculations for xi != 1, which diverges individually,
	! but cancels overall in the full amplitude; hence, we set it to zero to avoid the artificial divergence
	double complex function D0Mine(a,b,c,d,e,f,g,h,j,k)
        implicit none
        double precision, intent(in) :: a,b,c,d,e,f,g,h,j,k
        D0Mine = (0D0,0D0)
    end function D0Mine

    ! These are the derivatives of the three-point function with C0(0,p,p,0,0,m) w.r.t. p^2, appearing in calculations for xi != 1, 
	! which diverge individually, but which cancel overall in the full amplitude; hence, we set them to zero to avoid these artificial divergence
    double complex function DC01Mine(a,b,c,d,e,f)
        implicit none
        double precision, intent(in) :: a,b,c,d,e,f
        DC01Mine = (0D0,0D0)
    end function DC01Mine
    double complex function DC02Mine(a,b,c,d,e,f)
        implicit none
        double precision, intent(in) :: a,b,c,d,e,f
        DC02Mine = (0D0,0D0)
    end function DC02Mine

    double complex function DiracGamma(a)
        implicit none
        double precision, intent(in) :: a
        DiracGamma = (0D0, 0D0)
    end function DiracGamma

end module constants
