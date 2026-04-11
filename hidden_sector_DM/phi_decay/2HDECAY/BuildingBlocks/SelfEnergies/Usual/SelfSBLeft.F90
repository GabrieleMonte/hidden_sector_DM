double complex function SelfSBLeftUsual(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(9)

 amplitudes(1) = (0.015625D0*CKM12*CKMC13*EL2*(MU2*A0(MHp2) - 1.D0*MU2*A0(MU2) - 1.D0*MHp2*MU2*B0(x, MHp2, MU2) + MU2*x*B0(x, MHp&
  &2, MU2) + B0(x, MHp2, MU2)*DBLE(MU**INT(4.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(2) = (0.015625D0*CKM22*CKMC23*EL2*(-1.D0*MC2*A0(MC2) + MC2*A0(MHp2) - 1.D0*MC2*MHp2*B0(x, MC2, MHp2) + MC2*x*B0(x, MC&
  &2, MHp2) + B0(x, MC2, MHp2)*DBLE(MC**INT(4.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(3) = (0.015625D0*CKM32*CKMC33*EL2*(MT2*A0(MHp2) - 1.D0*MT2*A0(MT2) - 1.D0*MHp2*MT2*B0(x, MHp2, MT2) + MT2*x*B0(x, MHp&
  &2, MT2) + B0(x, MHp2, MT2)*DBLE(MT**INT(4.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(4) = (0.015625D0*CKM12*CKMC13*EL2*(-1.D0*MU2*A0(MU2) + MU2*A0(MW2) - 1.D0*MU2*MW2*B0(x, MU2, MW2) + MU2*x*B0(x, MU2, &
  &MW2) + B0(x, MU2, MW2)*DBLE(MU**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(5) = (0.015625D0*CKM22*CKMC23*EL2*(-1.D0*MC2*A0(MC2) + MC2*A0(MW2) - 1.D0*MC2*MW2*B0(x, MC2, MW2) + MC2*x*B0(x, MC2, &
  &MW2) + B0(x, MC2, MW2)*DBLE(MC**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(6) = (0.015625D0*CKM32*CKMC33*EL2*(-1.D0*MT2*A0(MT2) + MT2*A0(MW2) - 1.D0*MT2*MW2*B0(x, MT2, MW2) + MT2*x*B0(x, MT2, &
  &MW2) + B0(x, MT2, MW2)*DBLE(MT**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(7) = (0.03125D0*CKM12*CKMC13*EL2*(-1.D0*x - 1.D0*A0(MU2) + A0(MW2) + MU2*B0(x, MU2, MW2) - 1.D0*MW2*B0(x, MU2, MW2) +&
  & x*B0(x, MU2, MW2)))/(PI2*SW2*x)

 amplitudes(8) = (0.03125D0*CKM22*CKMC23*EL2*(-1.D0*x - 1.D0*A0(MC2) + A0(MW2) + MC2*B0(x, MC2, MW2) - 1.D0*MW2*B0(x, MC2, MW2) +&
  & x*B0(x, MC2, MW2)))/(PI2*SW2*x)

 amplitudes(9) = (0.03125D0*CKM32*CKMC33*EL2*(-1.D0*x - 1.D0*A0(MT2) + A0(MW2) + MT2*B0(x, MT2, MW2) - 1.D0*MW2*B0(x, MT2, MW2) +&
  & x*B0(x, MT2, MW2)))/(PI2*SW2*x)

  totalAmplitude = (0D0,0D0)
 do j=1,9
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 SelfSBLeftUsual = totalAmplitude
end function SelfSBLeftUsual

