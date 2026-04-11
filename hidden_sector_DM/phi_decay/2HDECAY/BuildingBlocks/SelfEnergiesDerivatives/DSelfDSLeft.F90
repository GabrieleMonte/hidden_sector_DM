double complex function DSelfDSLeft(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(9)

 amplitudes(1) = (0.015625D0*CKM11*CKMC12*EL2*(MU2*B0(x, MHp2, MU2) - 1.D0*MHp2*MU2*DB0(x, MHp2, MU2) + MU2*x*DB0(x, MHp2, MU2) +&
  & DB0(x, MHp2, MU2)*DBLE(MU**INT(4.D0))))/(MW2*PI2*SW2*TB2*x) - (0.015625D0*CKM11*CKMC12*EL2*(MU2*A0(MHp2) - 1.D0*MU2*A0(MU2) -&
  & 1.D0*MHp2*MU2*B0(x, MHp2, MU2) + MU2*x*B0(x, MHp2, MU2) + B0(x, MHp2, MU2)*DBLE(MU**INT(4.D0)))*DBLE(x**INT(-2.D0)))/(MW2*PI2&
  &*SW2*TB2)

 amplitudes(2) = (0.015625D0*CKM21*CKMC22*EL2*(MC2*B0(x, MC2, MHp2) - 1.D0*MC2*MHp2*DB0(x, MC2, MHp2) + MC2*x*DB0(x, MC2, MHp2) +&
  & DB0(x, MC2, MHp2)*DBLE(MC**INT(4.D0))))/(MW2*PI2*SW2*TB2*x) - (0.015625D0*CKM21*CKMC22*EL2*(-1.D0*MC2*A0(MC2) + MC2*A0(MHp2) &
  &- 1.D0*MC2*MHp2*B0(x, MC2, MHp2) + MC2*x*B0(x, MC2, MHp2) + B0(x, MC2, MHp2)*DBLE(MC**INT(4.D0)))*DBLE(x**INT(-2.D0)))/(MW2*PI&
  &2*SW2*TB2)

 amplitudes(3) = (0.015625D0*CKM31*CKMC32*EL2*(MT2*B0(x, MHp2, MT2) - 1.D0*MHp2*MT2*DB0(x, MHp2, MT2) + MT2*x*DB0(x, MHp2, MT2) +&
  & DB0(x, MHp2, MT2)*DBLE(MT**INT(4.D0))))/(MW2*PI2*SW2*TB2*x) - (0.015625D0*CKM31*CKMC32*EL2*(MT2*A0(MHp2) - 1.D0*MT2*A0(MT2) -&
  & 1.D0*MHp2*MT2*B0(x, MHp2, MT2) + MT2*x*B0(x, MHp2, MT2) + B0(x, MHp2, MT2)*DBLE(MT**INT(4.D0)))*DBLE(x**INT(-2.D0)))/(MW2*PI2&
  &*SW2*TB2)

 amplitudes(4) = (0.015625D0*CKM11*CKMC12*EL2*(MU2*B0(x, MU2, MW2) - 1.D0*MU2*MW2*DB0(x, MU2, MW2) + MU2*x*DB0(x, MU2, MW2) + DB0&
  &(x, MU2, MW2)*DBLE(MU**INT(4.D0))))/(MW2*PI2*SW2*x) - (0.015625D0*CKM11*CKMC12*EL2*(-1.D0*MU2*A0(MU2) + MU2*A0(MW2) - 1.D0*MU2&
  &*MW2*B0(x, MU2, MW2) + MU2*x*B0(x, MU2, MW2) + B0(x, MU2, MW2)*DBLE(MU**INT(4.D0)))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SW2)

 amplitudes(5) = (0.015625D0*CKM21*CKMC22*EL2*(MC2*B0(x, MC2, MW2) - 1.D0*MC2*MW2*DB0(x, MC2, MW2) + MC2*x*DB0(x, MC2, MW2) + DB0&
  &(x, MC2, MW2)*DBLE(MC**INT(4.D0))))/(MW2*PI2*SW2*x) - (0.015625D0*CKM21*CKMC22*EL2*(-1.D0*MC2*A0(MC2) + MC2*A0(MW2) - 1.D0*MC2&
  &*MW2*B0(x, MC2, MW2) + MC2*x*B0(x, MC2, MW2) + B0(x, MC2, MW2)*DBLE(MC**INT(4.D0)))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SW2)

 amplitudes(6) = (0.015625D0*CKM31*CKMC32*EL2*(MT2*B0(x, MT2, MW2) - 1.D0*MT2*MW2*DB0(x, MT2, MW2) + MT2*x*DB0(x, MT2, MW2) + DB0&
  &(x, MT2, MW2)*DBLE(MT**INT(4.D0))))/(MW2*PI2*SW2*x) - (0.015625D0*CKM31*CKMC32*EL2*(-1.D0*MT2*A0(MT2) + MT2*A0(MW2) - 1.D0*MT2&
  &*MW2*B0(x, MT2, MW2) + MT2*x*B0(x, MT2, MW2) + B0(x, MT2, MW2)*DBLE(MT**INT(4.D0)))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SW2)

 amplitudes(7) = (0.03125D0*CKM11*CKMC12*EL2*(-1.D0 + B0(x, MU2, MW2) + MU2*DB0(x, MU2, MW2) - 1.D0*MW2*DB0(x, MU2, MW2) + x*DB0(&
  &x, MU2, MW2)))/(PI2*SW2*x) - (0.03125D0*CKM11*CKMC12*EL2*(-1.D0*x - 1.D0*A0(MU2) + A0(MW2) + MU2*B0(x, MU2, MW2) - 1.D0*MW2*B0&
  &(x, MU2, MW2) + x*B0(x, MU2, MW2))* DBLE(x**INT(-2.D0)))/(PI2*SW2)

 amplitudes(8) = (0.03125D0*CKM21*CKMC22*EL2*(-1.D0 + B0(x, MC2, MW2) + MC2*DB0(x, MC2, MW2) - 1.D0*MW2*DB0(x, MC2, MW2) + x*DB0(&
  &x, MC2, MW2)))/(PI2*SW2*x) - (0.03125D0*CKM21*CKMC22*EL2*(-1.D0*x - 1.D0*A0(MC2) + A0(MW2) + MC2*B0(x, MC2, MW2) - 1.D0*MW2*B0&
  &(x, MC2, MW2) + x*B0(x, MC2, MW2))* DBLE(x**INT(-2.D0)))/(PI2*SW2)

 amplitudes(9) = (0.03125D0*CKM31*CKMC32*EL2*(-1.D0 + B0(x, MT2, MW2) + MT2*DB0(x, MT2, MW2) - 1.D0*MW2*DB0(x, MT2, MW2) + x*DB0(&
  &x, MT2, MW2)))/(PI2*SW2*x) - (0.03125D0*CKM31*CKMC32*EL2*(-1.D0*x - 1.D0*A0(MT2) + A0(MW2) + MT2*B0(x, MT2, MW2) - 1.D0*MW2*B0&
  &(x, MT2, MW2) + x*B0(x, MT2, MW2))* DBLE(x**INT(-2.D0)))/(PI2*SW2)

  totalAmplitude = (0D0,0D0)
 do j=1,9
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 DSelfDSLeft = totalAmplitude
end function DSelfDSLeft

