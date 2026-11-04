double complex function DSelfTTLeft(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(15)

 amplitudes(1) = (0.0078125D0*CA2*EL2*MT2*(B0(x, Mh02, MT2) - 1.D0*Mh02*DB0(x, Mh02, MT2) + MT2*DB0(x, Mh02, MT2) + x*DB0(x, Mh02&
  &, MT2)))/ (MW2*PI2*SB2*SW2*x) - (0.0078125D0*CA2*EL2*MT2*(A0(Mh02) - 1.D0*A0(MT2) - 1.D0*Mh02*B0(x, Mh02, MT2) + MT2*B0(x, Mh0&
  &2, MT2) + x*B0(x, Mh02, MT2))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SB2*SW2)

 amplitudes(2) = (0.0078125D0*EL2*MT2*SA2*(B0(x, MHH2, MT2) - 1.D0*MHH2*DB0(x, MHH2, MT2) + MT2*DB0(x, MHH2, MT2) + x*DB0(x, MHH2&
  &, MT2)))/ (MW2*PI2*SB2*SW2*x) - (0.0078125D0*EL2*MT2*SA2*(A0(MHH2) - 1.D0*A0(MT2) - 1.D0*MHH2*B0(x, MHH2, MT2) + MT2*B0(x, MHH&
  &2, MT2) + x*B0(x, MHH2, MT2))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SB2*SW2)

 amplitudes(3) = (0.0078125D0*EL2*MT2*(B0(x, MA02, MT2) - 1.D0*MA02*DB0(x, MA02, MT2) + MT2*DB0(x, MA02, MT2) + x*DB0(x, MA02, MT&
  &2)))/(MW2*PI2*SW2*TB2*x) - (0.0078125D0*EL2*MT2*(A0(MA02) - 1.D0*A0(MT2) - 1.D0*MA02*B0(x, MA02, MT2) + MT2*B0(x, MA02, MT2) +&
  & x*B0(x, MA02, MT2))*DBLE(x**INT(-2.D0)))/ (MW2*PI2*SW2*TB2)

 amplitudes(4) = (0.0078125D0*EL2*MT2*(B0(x, MT2, MZ2) + MT2*DB0(x, MT2, MZ2) - 1.D0*MZ2*DB0(x, MT2, MZ2) + x*DB0(x, MT2, MZ2)))/&
  &(MW2*PI2*SW2*x) - (0.0078125D0*EL2*MT2*(-1.D0*A0(MT2) + A0(MZ2) + MT2*B0(x, MT2, MZ2) - 1.D0*MZ2*B0(x, MT2, MZ2) + x*B0(x, MT2&
  &, MZ2))*DBLE(x**INT(-2.D0)))/ (MW2*PI2*SW2)

 amplitudes(5) = (-0.015625D0*CKM31*CKMC31*EL2*DBLE(x**INT(-2.D0))*(-1.D0*MD2*TB2*A0(MD2)*DBLE(Yuk3**INT(2.D0)) + MD2*TB2*A0(MHp2&
  &)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MD2*MHp2*TB2*B0(x, MD2, MHp2)*DBLE(Yuk3**INT(2.D0)) + MD2*TB2*x*B0(x, MD2, MHp2)*DBLE(Yuk3**INT&
  &(2.D0)) + TB2*B0(x, MD2, MHp2)*DBLE(MD**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2) + (0.015625D0*CKM31*CKMC31*EL2*(M&
  &D2*TB2*B0(x, MD2, MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MD2*MHp2*TB2*DB0(x, MD2, MHp2)*DBLE(Yuk3**INT(2.D0)) + MD2*TB2*x*DB0(x, M&
  &D2, MHp2)*DBLE(Yuk3**INT(2.D0)) + TB2*DB0(x, MD2, MHp2)*DBLE(MD**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(6) = (-0.015625D0*CKM32*CKMC32*EL2*DBLE(x**INT(-2.D0))*(MS2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MS2*TB2*A0(MS2)&
  &*DBLE(Yuk3**INT(2.D0)) - 1.D0*MHp2*MS2*TB2*B0(x, MHp2, MS2)*DBLE(Yuk3**INT(2.D0)) + MS2*TB2*x*B0(x, MHp2, MS2)*DBLE(Yuk3**INT(&
  &2.D0)) + TB2*B0(x, MHp2, MS2)*DBLE(MS**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2) + (0.015625D0*CKM32*CKMC32*EL2*(MS&
  &2*TB2*B0(x, MHp2, MS2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MHp2*MS2*TB2*DB0(x, MHp2, MS2)*DBLE(Yuk3**INT(2.D0)) + MS2*TB2*x*DB0(x, MH&
  &p2, MS2)*DBLE(Yuk3**INT(2.D0)) + TB2*DB0(x, MHp2, MS2)*DBLE(MS**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(7) = (-0.015625D0*CKM33*CKMC33*EL2*DBLE(x**INT(-2.D0))*(-1.D0*MB2*TB2*A0(MB2)*DBLE(Yuk3**INT(2.D0)) + MB2*TB2*A0(MHp2&
  &)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MB2*MHp2*TB2*B0(x, MB2, MHp2)*DBLE(Yuk3**INT(2.D0)) + MB2*TB2*x*B0(x, MB2, MHp2)*DBLE(Yuk3**INT&
  &(2.D0)) + TB2*B0(x, MB2, MHp2)*DBLE(MB**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2) + (0.015625D0*CKM33*CKMC33*EL2*(M&
  &B2*TB2*B0(x, MB2, MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MB2*MHp2*TB2*DB0(x, MB2, MHp2)*DBLE(Yuk3**INT(2.D0)) + MB2*TB2*x*DB0(x, M&
  &B2, MHp2)*DBLE(Yuk3**INT(2.D0)) + TB2*DB0(x, MB2, MHp2)*DBLE(MB**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(8) = (0.015625D0*CKM31*CKMC31*EL2*(MD2*B0(x, MD2, MW2) - 1.D0*MD2*MW2*DB0(x, MD2, MW2) + MD2*x*DB0(x, MD2, MW2) + DB0&
  &(x, MD2, MW2)*DBLE(MD**INT(4.D0))))/(MW2*PI2*SW2*x) - (0.015625D0*CKM31*CKMC31*EL2*(-1.D0*MD2*A0(MD2) + MD2*A0(MW2) - 1.D0*MD2&
  &*MW2*B0(x, MD2, MW2) + MD2*x*B0(x, MD2, MW2) + B0(x, MD2, MW2)*DBLE(MD**INT(4.D0)))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SW2)

 amplitudes(9) = (0.015625D0*CKM32*CKMC32*EL2*(MS2*B0(x, MS2, MW2) - 1.D0*MS2*MW2*DB0(x, MS2, MW2) + MS2*x*DB0(x, MS2, MW2) + DB0&
  &(x, MS2, MW2)*DBLE(MS**INT(4.D0))))/(MW2*PI2*SW2*x) - (0.015625D0*CKM32*CKMC32*EL2*(-1.D0*MS2*A0(MS2) + MS2*A0(MW2) - 1.D0*MS2&
  &*MW2*B0(x, MS2, MW2) + MS2*x*B0(x, MS2, MW2) + B0(x, MS2, MW2)*DBLE(MS**INT(4.D0)))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SW2)

 amplitudes(10) = (0.015625D0*CKM33*CKMC33*EL2*(MB2*B0(x, MB2, MW2) - 1.D0*MB2*MW2*DB0(x, MB2, MW2) + MB2*x*DB0(x, MB2, MW2) + DB&
  &0(x, MB2, MW2)*DBLE(MB**INT(4.D0))))/(MW2*PI2*SW2*x) - (0.015625D0*CKM33*CKMC33*EL2*(-1.D0*MB2*A0(MB2) + MB2*A0(MW2) - 1.D0*MB&
  &2*MW2*B0(x, MB2, MW2) + MB2*x*B0(x, MB2, MW2) + B0(x, MB2, MW2)*DBLE(MB**INT(4.D0)))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SW2)

 amplitudes(11) = (0.027777777777777776D0*EL2*(-1.D0 + B0(x, 0.D0, MT2) + MT2*DB0(x, 0.D0, MT2) + x*DB0(x, 0.D0, MT2)))/(PI2*x) -&
  & (0.027777777777777776D0*EL2*(-1.D0*x - 1.D0*A0(MT2) + MT2*B0(x, 0.D0, MT2) + x*B0(x, 0.D0, MT2))*DBLE(x**INT(-2.D0)))/PI2

 amplitudes(12) = (0.001736111111111111D0*EL2*(-9.D0 + 24.D0*SW2 + 9.D0*B0(x, MT2, MZ2) - 24.D0*SW2*B0(x, MT2, MZ2) + 9.D0*MT2*DB&
  &0(x, MT2, MZ2) - 9.D0*MZ2*DB0(x, MT2, MZ2) - 24.D0*MT2*SW2*DB0(x, MT2, MZ2) + 24.D0*MZ2*SW2*DB0(x, MT2, MZ2) + 9.D0*x*DB0(x, M&
  &T2, MZ2) - 24.D0*SW2*x*DB0(x, MT2, MZ2) - 16.D0*DBLE(SW**INT(4.D0)) + 16.D0*B0(x, MT2, MZ2)*DBLE(SW**INT(4.D0)) + 16.D0*MT2*DB&
  &0(x, MT2, MZ2)*DBLE(SW**INT(4.D0)) - 16.D0*MZ2*DB0(x, MT2, MZ2)*DBLE(SW**INT(4.D0)) + 16.D0*x*DB0(x, MT2, MZ2)*DBLE(SW**INT(4.&
  &D0))))/ (CW2*PI2*SW2*x) - (0.001736111111111111D0*EL2*(-9.D0*x + 24.D0*SW2*x - 9.D0*A0(MT2) + 24.D0*SW2*A0(MT2) + 9.D0*A0(MZ2)&
  & - 24.D0*SW2*A0(MZ2) + 9.D0*MT2*B0(x, MT2, MZ2) - 9.D0*MZ2*B0(x, MT2, MZ2) - 24.D0*MT2*SW2*B0(x, MT2, MZ2) + 24.D0*MZ2*SW2*B0(&
  &x, MT2, MZ2) + 9.D0*x*B0(x, MT2, MZ2) - 24.D0*SW2*x*B0(x, MT2, MZ2) - 16.D0*x*DBLE(SW**INT(4.D0)) - 16.D0*A0(MT2)*DBLE(SW**INT&
  &(4.D0)) + 16.D0*A0(MZ2)*DBLE(SW**INT(4.D0)) + 16.D0*MT2*B0(x, MT2, MZ2)*DBLE(SW**INT(4.D0)) - 16.D0*MZ2*B0(x, MT2, MZ2)*DBLE(S&
  &W**INT(4.D0)) + 16.D0*x*B0(x, MT2, MZ2)*DBLE(SW**INT(4.D0)))* DBLE(x**INT(-2.D0)))/(CW2*PI2*SW2)

 amplitudes(13) = (0.03125D0*CKM31*CKMC31*EL2*(-1.D0 + B0(x, MD2, MW2) + MD2*DB0(x, MD2, MW2) - 1.D0*MW2*DB0(x, MD2, MW2) + x*DB0&
  &(x, MD2, MW2)))/(PI2*SW2*x) - (0.03125D0*CKM31*CKMC31*EL2*(-1.D0*x - 1.D0*A0(MD2) + A0(MW2) + MD2*B0(x, MD2, MW2) - 1.D0*MW2*B&
  &0(x, MD2, MW2) + x*B0(x, MD2, MW2))* DBLE(x**INT(-2.D0)))/(PI2*SW2)

 amplitudes(14) = (0.03125D0*CKM32*CKMC32*EL2*(-1.D0 + B0(x, MS2, MW2) + MS2*DB0(x, MS2, MW2) - 1.D0*MW2*DB0(x, MS2, MW2) + x*DB0&
  &(x, MS2, MW2)))/(PI2*SW2*x) - (0.03125D0*CKM32*CKMC32*EL2*(-1.D0*x - 1.D0*A0(MS2) + A0(MW2) + MS2*B0(x, MS2, MW2) - 1.D0*MW2*B&
  &0(x, MS2, MW2) + x*B0(x, MS2, MW2))* DBLE(x**INT(-2.D0)))/(PI2*SW2)

 amplitudes(15) = (0.03125D0*CKM33*CKMC33*EL2*(-1.D0 + B0(x, MB2, MW2) + MB2*DB0(x, MB2, MW2) - 1.D0*MW2*DB0(x, MB2, MW2) + x*DB0&
  &(x, MB2, MW2)))/(PI2*SW2*x) - (0.03125D0*CKM33*CKMC33*EL2*(-1.D0*x - 1.D0*A0(MB2) + A0(MW2) + MB2*B0(x, MB2, MW2) - 1.D0*MW2*B&
  &0(x, MB2, MW2) + x*B0(x, MB2, MW2))* DBLE(x**INT(-2.D0)))/(PI2*SW2)

  totalAmplitude = (0D0,0D0)
 do j=1,15
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 DSelfTTLeft = totalAmplitude
end function DSelfTTLeft

