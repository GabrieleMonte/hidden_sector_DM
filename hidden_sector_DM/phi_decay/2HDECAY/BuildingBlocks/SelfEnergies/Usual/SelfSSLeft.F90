double complex function SelfSSLeftUsual(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(15)

 amplitudes(1) = (0.0078125D0*EL2*MS2*(A0(Mh02) - 1.D0*A0(MS2) - 1.D0*Mh02*B0(x, Mh02, MS2) + MS2*B0(x, Mh02, MS2) + x*B0(x, Mh02&
  &, MS2))*DBLE(Yuk1**INT(2.D0)))/ (MW2*PI2*SW2*x)

 amplitudes(2) = (0.0078125D0*EL2*MS2*(A0(MHH2) - 1.D0*A0(MS2) - 1.D0*MHH2*B0(x, MHH2, MS2) + MS2*B0(x, MHH2, MS2) + x*B0(x, MHH2&
  &, MS2))*DBLE(Yuk2**INT(2.D0)))/ (MW2*PI2*SW2*x)

 amplitudes(3) = (0.0078125D0*EL2*MS2*(A0(MA02) - 1.D0*A0(MS2) - 1.D0*MA02*B0(x, MA02, MS2) + MS2*B0(x, MA02, MS2) + x*B0(x, MA02&
  &, MS2))*DBLE(Yuk3**INT(2.D0)))/ (MW2*PI2*SW2*x)

 amplitudes(4) = (0.0078125D0*EL2*MS2*(-1.D0*A0(MS2) + A0(MZ2) + MS2*B0(x, MS2, MZ2) - 1.D0*MZ2*B0(x, MS2, MZ2) + x*B0(x, MS2, MZ&
  &2)))/(MW2*PI2*SW2*x)

 amplitudes(5) = (0.015625D0*CKM12*CKMC12*EL2*(MU2*A0(MHp2) - 1.D0*MU2*A0(MU2) - 1.D0*MHp2*MU2*B0(x, MHp2, MU2) + MU2*x*B0(x, MHp&
  &2, MU2) + B0(x, MHp2, MU2)*DBLE(MU**INT(4.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(6) = (0.015625D0*CKM22*CKMC22*EL2*(-1.D0*MC2*A0(MC2) + MC2*A0(MHp2) - 1.D0*MC2*MHp2*B0(x, MC2, MHp2) + MC2*x*B0(x, MC&
  &2, MHp2) + B0(x, MC2, MHp2)*DBLE(MC**INT(4.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(7) = (0.015625D0*CKM32*CKMC32*EL2*(MT2*A0(MHp2) - 1.D0*MT2*A0(MT2) - 1.D0*MHp2*MT2*B0(x, MHp2, MT2) + MT2*x*B0(x, MHp&
  &2, MT2) + B0(x, MHp2, MT2)*DBLE(MT**INT(4.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(8) = (0.015625D0*CKM12*CKMC12*EL2*(-1.D0*MU2*A0(MU2) + MU2*A0(MW2) - 1.D0*MU2*MW2*B0(x, MU2, MW2) + MU2*x*B0(x, MU2, &
  &MW2) + B0(x, MU2, MW2)*DBLE(MU**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(9) = (0.015625D0*CKM22*CKMC22*EL2*(-1.D0*MC2*A0(MC2) + MC2*A0(MW2) - 1.D0*MC2*MW2*B0(x, MC2, MW2) + MC2*x*B0(x, MC2, &
  &MW2) + B0(x, MC2, MW2)*DBLE(MC**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(10) = (0.015625D0*CKM32*CKMC32*EL2*(-1.D0*MT2*A0(MT2) + MT2*A0(MW2) - 1.D0*MT2*MW2*B0(x, MT2, MW2) + MT2*x*B0(x, MT2,&
  & MW2) + B0(x, MT2, MW2)*DBLE(MT**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(11) = (0.006944444444444444D0*EL2*(-1.D0*x - 1.D0*A0(MS2) + MS2*B0(x, 0.D0, MS2) + x*B0(x, 0.D0, MS2)))/(PI2*x)

 amplitudes(12) = (0.001736111111111111D0*EL2*(-9.D0*x + 12.D0*SW2*x - 9.D0*A0(MS2) + 12.D0*SW2*A0(MS2) + 9.D0*A0(MZ2) - 12.D0*SW&
  &2*A0(MZ2) + 9.D0*MS2*B0(x, MS2, MZ2) - 9.D0*MZ2*B0(x, MS2, MZ2) - 12.D0*MS2*SW2*B0(x, MS2, MZ2) + 12.D0*MZ2*SW2*B0(x, MS2, MZ2&
  &) + 9.D0*x*B0(x, MS2, MZ2) - 12.D0*SW2*x*B0(x, MS2, MZ2) - 4.D0*x*DBLE(SW**INT(4.D0)) - 4.D0*A0(MS2)*DBLE(SW**INT(4.D0)) + 4.D&
  &0*A0(MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*MS2*B0(x, MS2, MZ2)*DBLE(SW**INT(4.D0)) - 4.D0*MZ2*B0(x, MS2, MZ2)*DBLE(SW**INT(4.D0)) + &
  &4.D0*x*B0(x, MS2, MZ2)*DBLE(SW**INT(4.D0))))/ (CW2*PI2*SW2*x)

 amplitudes(13) = (0.03125D0*CKM12*CKMC12*EL2*(-1.D0*x - 1.D0*A0(MU2) + A0(MW2) + MU2*B0(x, MU2, MW2) - 1.D0*MW2*B0(x, MU2, MW2) &
  &+ x*B0(x, MU2, MW2)))/(PI2*SW2*x)

 amplitudes(14) = (0.03125D0*CKM22*CKMC22*EL2*(-1.D0*x - 1.D0*A0(MC2) + A0(MW2) + MC2*B0(x, MC2, MW2) - 1.D0*MW2*B0(x, MC2, MW2) &
  &+ x*B0(x, MC2, MW2)))/(PI2*SW2*x)

 amplitudes(15) = (0.03125D0*CKM32*CKMC32*EL2*(-1.D0*x - 1.D0*A0(MT2) + A0(MW2) + MT2*B0(x, MT2, MW2) - 1.D0*MW2*B0(x, MT2, MW2) &
  &+ x*B0(x, MT2, MW2)))/(PI2*SW2*x)

  totalAmplitude = (0D0,0D0)
 do j=1,15
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 SelfSSLeftUsual = totalAmplitude
end function SelfSSLeftUsual

