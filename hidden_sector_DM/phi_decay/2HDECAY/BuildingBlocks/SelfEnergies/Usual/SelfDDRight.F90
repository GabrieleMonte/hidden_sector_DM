double complex function SelfDDRightUsual(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(15)

 amplitudes(1) = (0.0078125D0*EL2*MD2*(-1.D0*A0(MD2) + A0(Mh02) + MD2*B0(x, MD2, Mh02) - 1.D0*Mh02*B0(x, MD2, Mh02) + x*B0(x, MD2&
  &, Mh02))*DBLE(Yuk1**INT(2.D0)))/ (MW2*PI2*SW2*x)

 amplitudes(2) = (0.0078125D0*EL2*MD2*(-1.D0*A0(MD2) + A0(MHH2) + MD2*B0(x, MD2, MHH2) - 1.D0*MHH2*B0(x, MD2, MHH2) + x*B0(x, MD2&
  &, MHH2))*DBLE(Yuk2**INT(2.D0)))/ (MW2*PI2*SW2*x)

 amplitudes(3) = (0.0078125D0*EL2*MD2*(A0(MA02) - 1.D0*A0(MD2) - 1.D0*MA02*B0(x, MA02, MD2) + MD2*B0(x, MA02, MD2) + x*B0(x, MA02&
  &, MD2))*DBLE(Yuk3**INT(2.D0)))/ (MW2*PI2*SW2*x)

 amplitudes(4) = (0.0078125D0*EL2*MD2*(-1.D0*A0(MD2) + A0(MZ2) + MD2*B0(x, MD2, MZ2) - 1.D0*MZ2*B0(x, MD2, MZ2) + x*B0(x, MD2, MZ&
  &2)))/(MW2*PI2*SW2*x)

 amplitudes(5) = (0.015625D0*CKM11*CKMC11*EL2*(MD2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MD2*TB2*A0(MU2)*DBLE(Yuk3**INT(2.D0)&
  &) - 1.D0*MD2*MHp2*TB2*B0(x, MHp2, MU2)*DBLE(Yuk3**INT(2.D0)) + MD2*MU2*TB2*B0(x, MHp2, MU2)*DBLE(Yuk3**INT(2.D0)) + MD2*TB2*x*&
  &B0(x, MHp2, MU2)*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(6) = (0.015625D0*CKM21*CKMC21*EL2*(-1.D0*MD2*TB2*A0(MC2)*DBLE(Yuk3**INT(2.D0)) + MD2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D0&
  &)) + MC2*MD2*TB2*B0(x, MC2, MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MD2*MHp2*TB2*B0(x, MC2, MHp2)*DBLE(Yuk3**INT(2.D0)) + MD2*TB2*x&
  &*B0(x, MC2, MHp2)*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(7) = (0.015625D0*CKM31*CKMC31*EL2*(MD2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MD2*TB2*A0(MT2)*DBLE(Yuk3**INT(2.D0)&
  &) - 1.D0*MD2*MHp2*TB2*B0(x, MHp2, MT2)*DBLE(Yuk3**INT(2.D0)) + MD2*MT2*TB2*B0(x, MHp2, MT2)*DBLE(Yuk3**INT(2.D0)) + MD2*TB2*x*&
  &B0(x, MHp2, MT2)*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(8) = (0.015625D0*CKM11*CKMC11*EL2*(-1.D0*MD2*A0(MU2) + MD2*A0(MW2) + MD2*MU2*B0(x, MU2, MW2) - 1.D0*MD2*MW2*B0(x, MU2&
  &, MW2) + MD2*x*B0(x, MU2, MW2)))/(MW2*PI2*SW2*x)

 amplitudes(9) = (0.015625D0*CKM21*CKMC21*EL2*(-1.D0*MD2*A0(MC2) + MD2*A0(MW2) + MC2*MD2*B0(x, MC2, MW2) - 1.D0*MD2*MW2*B0(x, MC2&
  &, MW2) + MD2*x*B0(x, MC2, MW2)))/(MW2*PI2*SW2*x)

 amplitudes(10) = (0.015625D0*CKM31*CKMC31*EL2*(-1.D0*MD2*A0(MT2) + MD2*A0(MW2) + MD2*MT2*B0(x, MT2, MW2) - 1.D0*MD2*MW2*B0(x, MT&
  &2, MW2) + MD2*x*B0(x, MT2, MW2)))/(MW2*PI2*SW2*x)

 amplitudes(11) = (0.006944444444444444D0*EL2*(-1.D0*x - 1.D0*A0(MD2) + MD2*B0(x, 0.D0, MD2) + x*B0(x, 0.D0, MD2)))/(PI2*x)

 amplitudes(12) = (0.001736111111111111D0*EL2*(-4.D0*x*DBLE(SW**INT(4.D0)) - 4.D0*A0(MD2)*DBLE(SW**INT(4.D0)) + 4.D0*A0(MZ2)*DBLE&
  &(SW**INT(4.D0)) + 4.D0*MD2*B0(x, MD2, MZ2)*DBLE(SW**INT(4.D0)) - 4.D0*MZ2*B0(x, MD2, MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*x*B0(x, M&
  &D2, MZ2)*DBLE(SW**INT(4.D0))))/ (CW2*PI2*SW2*x)

 amplitudes(13) = 0.D0

 amplitudes(14) = 0.D0

 amplitudes(15) = 0.D0

  totalAmplitude = (0D0,0D0)
 do j=1,15
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 SelfDDRightUsual = totalAmplitude
end function SelfDDRightUsual

