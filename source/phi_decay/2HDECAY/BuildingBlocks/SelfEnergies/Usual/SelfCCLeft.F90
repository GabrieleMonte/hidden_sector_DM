double complex function SelfCCLeftUsual(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(15)

 amplitudes(1) = (0.0078125D0*CA2*EL2*MC2*(-1.D0*A0(MC2) + A0(Mh02) + MC2*B0(x, MC2, Mh02) - 1.D0*Mh02*B0(x, MC2, Mh02) + x*B0(x,&
  & MC2, Mh02)))/ (MW2*PI2*SB2*SW2*x)

 amplitudes(2) = (0.0078125D0*EL2*MC2*SA2*(-1.D0*A0(MC2) + A0(MHH2) + MC2*B0(x, MC2, MHH2) - 1.D0*MHH2*B0(x, MC2, MHH2) + x*B0(x,&
  & MC2, MHH2)))/ (MW2*PI2*SB2*SW2*x)

 amplitudes(3) = (0.0078125D0*EL2*MC2*(A0(MA02) - 1.D0*A0(MC2) - 1.D0*MA02*B0(x, MA02, MC2) + MC2*B0(x, MA02, MC2) + x*B0(x, MA02&
  &, MC2)))/(MW2*PI2*SW2*TB2*x)

 amplitudes(4) = (0.0078125D0*EL2*MC2*(-1.D0*A0(MC2) + A0(MZ2) + MC2*B0(x, MC2, MZ2) - 1.D0*MZ2*B0(x, MC2, MZ2) + x*B0(x, MC2, MZ&
  &2)))/(MW2*PI2*SW2*x)

 amplitudes(5) = (0.015625D0*CKM21*CKMC21*EL2*(-1.D0*MD2*TB2*A0(MD2)*DBLE(Yuk3**INT(2.D0)) + MD2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D0&
  &)) - 1.D0*MD2*MHp2*TB2*B0(x, MD2, MHp2)*DBLE(Yuk3**INT(2.D0)) + MD2*TB2*x*B0(x, MD2, MHp2)*DBLE(Yuk3**INT(2.D0)) + TB2*B0(x, M&
  &D2, MHp2)*DBLE(MD**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(6) = (0.015625D0*CKM22*CKMC22*EL2*(MS2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MS2*TB2*A0(MS2)*DBLE(Yuk3**INT(2.D0)&
  &) - 1.D0*MHp2*MS2*TB2*B0(x, MHp2, MS2)*DBLE(Yuk3**INT(2.D0)) + MS2*TB2*x*B0(x, MHp2, MS2)*DBLE(Yuk3**INT(2.D0)) + TB2*B0(x, MH&
  &p2, MS2)*DBLE(MS**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(7) = (0.015625D0*CKM23*CKMC23*EL2*(-1.D0*MB2*TB2*A0(MB2)*DBLE(Yuk3**INT(2.D0)) + MB2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D0&
  &)) - 1.D0*MB2*MHp2*TB2*B0(x, MB2, MHp2)*DBLE(Yuk3**INT(2.D0)) + MB2*TB2*x*B0(x, MB2, MHp2)*DBLE(Yuk3**INT(2.D0)) + TB2*B0(x, M&
  &B2, MHp2)*DBLE(MB**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(8) = (0.015625D0*CKM21*CKMC21*EL2*(-1.D0*MD2*A0(MD2) + MD2*A0(MW2) - 1.D0*MD2*MW2*B0(x, MD2, MW2) + MD2*x*B0(x, MD2, &
  &MW2) + B0(x, MD2, MW2)*DBLE(MD**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(9) = (0.015625D0*CKM22*CKMC22*EL2*(-1.D0*MS2*A0(MS2) + MS2*A0(MW2) - 1.D0*MS2*MW2*B0(x, MS2, MW2) + MS2*x*B0(x, MS2, &
  &MW2) + B0(x, MS2, MW2)*DBLE(MS**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(10) = (0.015625D0*CKM23*CKMC23*EL2*(-1.D0*MB2*A0(MB2) + MB2*A0(MW2) - 1.D0*MB2*MW2*B0(x, MB2, MW2) + MB2*x*B0(x, MB2,&
  & MW2) + B0(x, MB2, MW2)*DBLE(MB**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(11) = (0.027777777777777776D0*EL2*(-1.D0*x - 1.D0*A0(MC2) + MC2*B0(x, 0.D0, MC2) + x*B0(x, 0.D0, MC2)))/(PI2*x)

 amplitudes(12) = (0.001736111111111111D0*EL2*(-9.D0*x + 24.D0*SW2*x - 9.D0*A0(MC2) + 24.D0*SW2*A0(MC2) + 9.D0*A0(MZ2) - 24.D0*SW&
  &2*A0(MZ2) + 9.D0*MC2*B0(x, MC2, MZ2) - 9.D0*MZ2*B0(x, MC2, MZ2) - 24.D0*MC2*SW2*B0(x, MC2, MZ2) + 24.D0*MZ2*SW2*B0(x, MC2, MZ2&
  &) + 9.D0*x*B0(x, MC2, MZ2) - 24.D0*SW2*x*B0(x, MC2, MZ2) - 16.D0*x*DBLE(SW**INT(4.D0)) - 16.D0*A0(MC2)*DBLE(SW**INT(4.D0)) + 1&
  &6.D0*A0(MZ2)*DBLE(SW**INT(4.D0)) + 16.D0*MC2*B0(x, MC2, MZ2)*DBLE(SW**INT(4.D0)) - 16.D0*MZ2*B0(x, MC2, MZ2)*DBLE(SW**INT(4.D0&
  &)) + 16.D0*x*B0(x, MC2, MZ2)*DBLE(SW**INT(4.D0))))/ (CW2*PI2*SW2*x)

 amplitudes(13) = (0.03125D0*CKM21*CKMC21*EL2*(-1.D0*x - 1.D0*A0(MD2) + A0(MW2) + MD2*B0(x, MD2, MW2) - 1.D0*MW2*B0(x, MD2, MW2) &
  &+ x*B0(x, MD2, MW2)))/(PI2*SW2*x)

 amplitudes(14) = (0.03125D0*CKM22*CKMC22*EL2*(-1.D0*x - 1.D0*A0(MS2) + A0(MW2) + MS2*B0(x, MS2, MW2) - 1.D0*MW2*B0(x, MS2, MW2) &
  &+ x*B0(x, MS2, MW2)))/(PI2*SW2*x)

 amplitudes(15) = (0.03125D0*CKM23*CKMC23*EL2*(-1.D0*x - 1.D0*A0(MB2) + A0(MW2) + MB2*B0(x, MB2, MW2) - 1.D0*MW2*B0(x, MB2, MW2) &
  &+ x*B0(x, MB2, MW2)))/(PI2*SW2*x)

  totalAmplitude = (0D0,0D0)
 do j=1,15
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 SelfCCLeftUsual = totalAmplitude
end function SelfCCLeftUsual

