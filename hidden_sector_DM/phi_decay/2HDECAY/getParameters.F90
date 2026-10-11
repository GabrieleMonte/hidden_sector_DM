subroutine getParameters(OStype, param2HDMonly)
    use constants
    implicit none
    character(50) dump, dump2, dump3
    ! double precision dump2
    double precision MHHTemp, MA0Temp, MHpTemp, alphaTemp, TBTemp, m12squaredTemp
    integer TypeOf2HDMTemp    
    integer statOpen, statRead
    integer :: currentLine = 1
    logical fileExists
    character isContinue
    integer, intent(in) :: OStype, param2HDMonly
    character(2) :: OSpathSeparator
    character(300), parameter :: pathToInputFiles = 'HDECAY'
    integer m
    double precision M11SqPot, M22SqPot, M12SqPot, tmpMass

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

    ! Check if we want all parameters or only the 2HDM ones
    if (param2HDMonly .EQ. 0) then
        ! Get the correct OS path separator
        if (OStype .eq. 0) then
            OSpathSeparator = '\\'
        else 
            OSpathSeparator = '/ '
        end if

        ! Check if the input file exists
        inquire(file=trim(trim(pathToInputFiles)//trim(OSpathSeparator))//'hdecay.in', exist=fileExists)
        if (.not. fileExists) then
            do
                print *, "ERROR: Could not find input file!"
                print *, ">>> Do you want to continue with the evaluation of the program? [y/n]"
                read (*,*) isContinue
                if (isContinue == 'n') then
                    print *, "Termination requested by user. 2HDECAY will be terminated now."
                    stop
                else if (isContinue == 'y') then
                    exit
                else
                    print *, "Invalid character. Enter y or n."
                end if
            end do
        end if

        ! Read the input file
        open(unit=42, file=trim(trim(pathToInputFiles)//trim(OSpathSeparator))//'hdecay.in', iostat=statOpen)
            if (statOpen == 0) then

                ! Skip the first five lines
                do m = 1, 5, 1
                    read(42, *)
                end do

                ! Read if the electroweak corrections should be calculated or not
                read(42, *) dump, dump2, omitELCorr

                ! In the input file, the SM block does begin at line 18
                do m = 1, 12, 1
                    read(42, *)
                end do

                ! Read the s quark mass
                read(42, *) dump, dump2, MS

                ! Skip two lines
                read(42, *)
                read(42, *)

                ! Read the c quark mass
                read(42, *) dump, dump2, MC

                ! Read the b quark mass
                read(42, *) dump, dump2, MB

                ! Read the t quark mass
                read(42, *) dump, dump2, MT

                ! Read the tau lepton mass
                read(42, *) dump, dump2, ML

                ! Read the muon mass
                read(42, *) dump, dump2, MM

                ! Skip a line
                read(42, *)
                
                ! Read the inverse fine-structure constant at the Z mass
                read(42, *) dump, dump2, alphaAtMZ

                ! Read the Fermi constant for consistency checks
                read(42, *) dump, dump2, GFinput
                
                ! Skip three lines
                read(42, *)
                read(42, *)
                read(42, *)

                ! Read the Z boson mass
                read(42, *) dump, dump2, MZ

                ! Read the W boson mass
                read(42, *) dump, dump2, MW

                ! Read the CKM elements
                read(42, *) dump, dump2, CKM33
                read(42, *) dump, dump2, CKM32
                read(42, *) dump, dump2, CKM31
                read(42, *) dump, dump2, CKM23
                read(42, *) dump, dump2, CKM22
                read(42, *) dump, dump2, CKM21
                read(42, *) dump, dump2, CKM13
                read(42, *) dump, dump2, CKM12
                read(42, *) dump, dump2, CKM11

                ! Skip 14 lines (after that, the 2HDM input block begins)
                do m = 1, 14, 1
                    read(42, *)
                end do

                ! Read the parameter type (1: masses and alpha are given; 2: lambda1 to lambda5 are given)
                read(42, *) dump, dump2, parameterType

                ! Read the 2HDM type
                read(42, *) dump, dump2, TypeOf2HDM

                ! Read the renormalization scheme
                read(42, *) dump, dump2, RenormScheme

                ! Read the reference scheme
                read(42, *) dump, dump2, RefScheme

                ! Skip a line
                read(42, *)

                ! Read Tan(beta)
                read(42, *) dump, TB
                
                ! Read m12^2
                read(42, *) dump, dump2, m12squared

                ! Read the input and output scales
                read(42, *) dump, dump2, InputScale
                read(42, *) dump, dump2, OutputScaleReadIn

                ! Skip a line
                read(42, *)

                ! Read alpha (relevant if parameter type = 1)
                read(42, *) dump, dump2, alpha

                ! Read Mh0 (relevant if parameter type = 1)
                read(42, *) dump, dump2, Mh0

                ! Read MHH (relevant if parameter type = 1)
                read(42, *) dump, dump2, MHH

                ! Read MA0 (relevant if parameter type = 1)
                read(42, *) dump, dump2, MA0

                ! Read MHp (relevant if parameter type = 1)
                read(42, *) dump, dump2, MHp

                ! Skip a line
                read(42, *) dump

                ! Read lambda1 (relevant if parameter type = 2)
                read(42, *) dump, dump2, hdecayLam1

                ! Read lambda2 (relevant if parameter type = 2)
                read(42, *) dump, dump2, hdecayLam2

                ! Read lambda3 (relevant if parameter type = 2)
                read(42, *) dump, dump2, hdecayLam3

                ! Read lambda4 (relevant if parameter type = 2)
                read(42, *) dump, dump2, hdecayLam4

                ! Read lambda5 (relevant if parameter type = 2)
                read(42, *) dump, dump2, hdecayLam5

            else
                do
                    print *, "ERROR: Generic error when reading input file!"
                    print *, ">>> Do you want to continue with the evaluation of the program? [y/n]"
                    read (*,*) isContinue
                    if (isContinue == 'n') then
                        print *, "Termination requested by user. 2HDECAY will be terminated now."
                        stop
                    else if (isContinue == 'y') then
                        exit
                    else
                        print *, "Invalid character. Enter y or n."
                    end if
                end do
            end if
        close(42)
        
        ! The light quark masses MU and MD and the electron mass ME are set to zero in HDECAY; we use very small values here to avoid numerical instability
        ! IMPORTANT: especially the light quarks u and d should not have exactly the same small mass, because in this case, terms like 1/(MU2 - MD2) become numerically unstable
        MU = 1.0e-5
        MD = 1.5e-5
        ME = 1.8e-5

        ! Calculate the sine and cosine of the Weinberg angle through MW and MZ 
        CW = MW/MZ 
        SW = DSQRT(1D0 - CW**2)

        ! Take the inverse of the read-in fine-structure constant (in the input file, we provide 1/alphaAtMZ, not alphaAtMZ!)
        ! alphaAtMZ = 1D0/alphaAtMZ

        ! Calculate the electromagnetic coupling constant out of the fundamental constants alphaAtMZ, MZ, MW 
        ! (This set is our convention; in HDECAY, we use GF instead of alphaAtMZ. The difference between the two is minimal, 
        ! and we can use the tree-level formula EL = DSQRT(8D0*GFermi*SW**2*MW**2/DSQRT(2D0)) to switch from our scheme to the one used in HDECAY)
        EL = DSQRT(4D0*PI*alphaAtMZ)

        ! Set the complex values of the CKM matrix to the non-complex ones (we do not consider CP violation)
        CKMC11 = CKM11
        CKMC12 = CKM12
        CKMC13 = CKM13
        CKMC21 = CKM21
        CKMC22 = CKM22
        CKMC23 = CKM23
        CKMC31 = CKM31
        CKMC32 = CKM32
        CKMC33 = CKM33

        ! Calculate the current beta with the given tan(beta)
        beta = datan(TB)
        CB = dcos(beta)
        SB = dsin(beta)
        S2B = dsin(2D0*beta)
        C2B = dcos(2D0*beta)

        ! Calculate the vev for the conversion of the masses
        ! vevCalc = 1D0/DSQRT(DSQRT(2D0)*GFinput)
        vevCalc = 2D0*MW*SW/EL

        ! If parameter type 1 is chosen, then alpha and the scalar masses are the relevant parameters and the potential parameters lamdba1 to lambda5 are calculated
        ! If parameter type 2 is chosen, then lambda1 to lambda5 are the relevant parameters and alpha and the scalar masses are calculated through these lambdas
        ! The conversion from alpha and the masses to the lambdas and vice versa is described in hep-ph/0408364
        if (parameterType .eq. 1) then
            CA = dcos(alpha)
            SA = dsin(alpha)
            S2A = dsin(2D0*alpha)
            hdecayLam1 = 1D0/vevCalc**2/CB**2*( -SB**2*(2D0*m12squared/S2B) + SA**2*Mh0**2 + CA**2*MHH**2 )
            hdecayLam2 = 1D0/vevCalc**2/SB**2*( -CB**2*(2D0*m12squared/S2B) + CA**2*Mh0**2 + SA**2*MHH**2 )
            hdecayLam3 = 1D0/vevCalc**2*( -(2D0*m12squared/S2B) + 2D0*MHp**2 + S2A/S2B*(MHH**2 - Mh0**2) )
            hdecayLam4 = 1D0/vevCalc**2*( (2D0*m12squared/S2B) + MA0**2 - 2D0*MHp**2 )
            hdecayLam5 = 1D0/vevCalc**2*( (2D0*m12squared/S2B) - MA0**2 )
        else if (parameterType .eq. 2) then
            ! CP-even mass matrix elements
            M11SqPot = (hdecayLam1*CB**4 + hdecayLam2*SB**4 + 2D0*(hdecayLam3 + hdecayLam4 + hdecayLam5)*CB**2*SB**2)*vevCalc**2
            M12SqPot = ( -hdecayLam1*CB**2 + hdecayLam2*SB**2 + (hdecayLam3 + hdecayLam4 + hdecayLam5)*(CB**2 - SB**2) )*CB*SB*&
                        &vevCalc**2
            M22SqPot = (2D0*m12squared/S2B) + 1D0/8D0*( hdecayLam1 + hdecayLam2 - 2D0*(hdecayLam3 + hdecayLam4 + hdecayLam5) )*&
                        &(1D0 - dcos(4D0*beta))*vevCalc**2

            ! Calculate the mixing angle alpha
            alpha = 1D0/2D0*datan( 2D0*M12SqPot/(M11SqPot - M22SqPot) ) + beta

            ! Calculate the CP-even masses
            Mh0 = DSQRT( (dsin(alpha-beta))**2*M11SqPot - dsin(2D0*(alpha-beta))*M12SqPot + (dcos(alpha-beta))**2*M22SqPot )
            MHH = DSQRT( (dcos(alpha-beta))**2*M11SqPot + dsin(2D0*(alpha-beta))*M12SqPot + (dsin(alpha-beta))**2*M22SqPot )

            ! If the mass hierarchy should be inverted, we switch the mass conventions such that h0 is still the smaller mass
            if(Mh0 .gt. MHH) then
                alpha = alpha - Pi/2D0 
                tmpMass = MHH
                MHH = Mh0 
                Mh0 = tmpMass
            endif

            ! Calculate the charged and CP-odd Higgs masses 
            MHp = DSQRT( (2D0*m12squared/S2B) - 1D0/2D0*(hdecayLam4 + hdecayLam5)*vevCalc**2 )
            MA0 = DSQRT( (2D0*m12squared/S2B) - hdecayLam5*vevCalc**2 )
        else
            write(*,*) "Error: unknown parameter type. Please choose the integers 1 or 2 for the parameter type in the input file!"
        end if

        ! Generate the square of the input parameters
        Mh02 = Mh0**2
        MW2 = MW**2
        MZ2 = MZ**2
        ME2 = ME**2
        MM2 = MM**2
        ML2 = ML**2
        MU2 = MU**2
        MC2 = MC**2
        MT2 = MT**2
        MD2 = MD**2
        MS2 = MS**2
        MB2 = MB**2
        EL2 = EL**2
        SW2 = SW**2
        CW2 = CW**2

        ! Calculate the rest of the 2HDM parameters
        CA = dcos(alpha)
        SA = dsin(alpha)
        TA = dtan(alpha)
        S2A = dsin(2D0*alpha)
        C2A = dcos(2D0*alpha)
        CAB = dcos(alpha + beta)
        SAB = dsin(alpha + beta)
        CBA = dcos(beta - alpha)
        SBA = dsin(beta - alpha)
        Lambda5 = EL2*m12squared/(2D0*SW2*MW2*SB*CB)
        if (TypeOf2HDM == 1) then
            Yuk1 = CA/SB
            Yuk2 = SA/SB
            Yuk3 = - 1D0/TB
            Yuk4 = CA/SB
            Yuk5 = SA/SB
            Yuk6 = - 1D0/TB
        else if (TypeOf2HDM == 2) then
            Yuk1 = - SA/CB
            Yuk2 = CA/CB
            Yuk3 = TB
            Yuk4 = - SA/CB
            Yuk5 = CA/CB
            Yuk6 = TB
        else if (TypeOf2HDM == 3) then
            Yuk1 = CA/SB
            Yuk2 = SA/SB
            Yuk3 = - 1D0/TB
            Yuk4 = - SA/CB
            Yuk5 = CA/CB
            Yuk6 = TB
        else
            Yuk1 = - SA/CB
            Yuk2 = CA/CB
            Yuk3 = TB
            Yuk4 = CA/SB
            Yuk5 = SA/SB
            Yuk6 = - 1D0/TB
        end if

        ! Generate the square of the additional input parameters
        MHH2 = MHH**2
        MA02 = MA0**2
        MHp2 = MHp**2
        CA2 = CA**2
        SA2 = SA**2
        TA2 = TA**2
        TB2 = TB**2
        CB2 = CB**2
        SB2 = SB**2
        S2A2 = S2A**2
        C2A2 = C2A**2
        S2B2 = S2B**2
        C2B2 = C2B**2
        CAB2 = CAB**2
        SAB2 = SAB**2
        CBA2 = CBA**2
        SBA2 = SBA**2
    else
        vevCalc = 2D0*MW*SW/EL

        CB = dcos(beta)
        SB = dsin(beta)
        TB = dtan(beta)
        S2B = dsin(2D0*beta)
        C2B = dcos(2D0*beta)

        CA = dcos(alpha)
        SA = dsin(alpha)
        TA = dtan(alpha)
        S2A = dsin(2D0*alpha)
        C2A = dcos(2D0*alpha)

        CAB = dcos(alpha + beta)
        SAB = dsin(alpha + beta)
        CBA = dcos(beta - alpha)
        SBA = dsin(beta - alpha)

        hdecayLam1 = 1D0/vevCalc**2/CB**2*( -SB**2*(2D0*m12squared/S2B) + SA**2*Mh0**2 + CA**2*MHH**2 )
        hdecayLam2 = 1D0/vevCalc**2/SB**2*( -CB**2*(2D0*m12squared/S2B) + CA**2*Mh0**2 + SA**2*MHH**2 )
        hdecayLam3 = 1D0/vevCalc**2*( -(2D0*m12squared/S2B) + 2D0*MHp**2 + S2A/S2B*(MHH**2 - Mh0**2) )
        hdecayLam4 = 1D0/vevCalc**2*( (2D0*m12squared/S2B) + MA0**2 - 2D0*MHp**2 )
        hdecayLam5 = 1D0/vevCalc**2*( (2D0*m12squared/S2B) - MA0**2 )
        
        Lambda5 = EL2*m12squared/(2D0*SW2*MW2*SB*CB)
        if (TypeOf2HDM == 1) then
            Yuk1 = CA/SB
            Yuk2 = SA/SB
            Yuk3 = - 1D0/TB
            Yuk4 = CA/SB
            Yuk5 = SA/SB
            Yuk6 = - 1D0/TB
        else if (TypeOf2HDM == 2) then
            Yuk1 = - SA/CB
            Yuk2 = CA/CB
            Yuk3 = TB
            Yuk4 = - SA/CB
            Yuk5 = CA/CB
            Yuk6 = TB
        else if (TypeOf2HDM == 3) then
            Yuk1 = CA/SB
            Yuk2 = SA/SB
            Yuk3 = - 1D0/TB
            Yuk4 = - SA/CB
            Yuk5 = CA/CB
            Yuk6 = TB
        else
            Yuk1 = - SA/CB
            Yuk2 = CA/CB
            Yuk3 = TB
            Yuk4 = CA/SB
            Yuk5 = SA/SB
            Yuk6 = - 1D0/TB
        end if

        CA2 = CA**2
        SA2 = SA**2
        TA2 = TA**2
        TB2 = TB**2
        CB2 = CB**2
        SB2 = SB**2
        S2A2 = S2A**2
        C2A2 = C2A**2
        S2B2 = S2B**2
        C2B2 = C2B**2
        CAB2 = CAB**2
        SAB2 = SAB**2
        CBA2 = CBA**2
        SBA2 = SBA**2

    end if

end subroutine getParameters
