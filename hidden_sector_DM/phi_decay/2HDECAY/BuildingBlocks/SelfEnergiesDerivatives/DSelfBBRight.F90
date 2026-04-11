double complex function DSelfBBRight(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(15)

 amplitudes(1) = (0.0078125D0*EL2*MB2*(B0(x, MB2, Mh02) + MB2*DB0(x, MB2, Mh02) - 1.D0*Mh02*DB0(x, MB2, Mh02) + x*DB0(x, MB2, Mh0&
  &2))*DBLE(Yuk1**INT(2.D0)))/ (MW2*PI2*SW2*x) - (0.0078125D0*EL2*MB2*(-1.D0*A0(MB2) + A0(Mh02) + MB2*B0(x, MB2, Mh02) - 1.D0*Mh0&
  &2*B0(x, MB2, Mh02) + x*B0(x, MB2, Mh02))* DBLE(x**INT(-2.D0))*DBLE(Yuk1**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(2) = (0.0078125D0*EL2*MB2*(B0(x, MB2, MHH2) + MB2*DB0(x, MB2, MHH2) - 1.D0*MHH2*DB0(x, MB2, MHH2) + x*DB0(x, MB2, MHH&
  &2))*DBLE(Yuk2**INT(2.D0)))/ (MW2*PI2*SW2*x) - (0.0078125D0*EL2*MB2*(-1.D0*A0(MB2) + A0(MHH2) + MB2*B0(x, MB2, MHH2) - 1.D0*MHH&
  &2*B0(x, MB2, MHH2) + x*B0(x, MB2, MHH2))* DBLE(x**INT(-2.D0))*DBLE(Yuk2**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(3) = (0.0078125D0*EL2*MB2*(B0(x, MA02, MB2) - 1.D0*MA02*DB0(x, MA02, MB2) + MB2*DB0(x, MA02, MB2) + x*DB0(x, MA02, MB&
  &2))*DBLE(Yuk3**INT(2.D0)))/ (MW2*PI2*SW2*x) - (0.0078125D0*EL2*MB2*(A0(MA02) - 1.D0*A0(MB2) - 1.D0*MA02*B0(x, MA02, MB2) + MB2&
  &*B0(x, MA02, MB2) + x*B0(x, MA02, MB2))* DBLE(x**INT(-2.D0))*DBLE(Yuk3**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(4) = (0.0078125D0*EL2*MB2*(B0(x, MB2, MZ2) + MB2*DB0(x, MB2, MZ2) - 1.D0*MZ2*DB0(x, MB2, MZ2) + x*DB0(x, MB2, MZ2)))/&
  &(MW2*PI2*SW2*x) - (0.0078125D0*EL2*MB2*(-1.D0*A0(MB2) + A0(MZ2) + MB2*B0(x, MB2, MZ2) - 1.D0*MZ2*B0(x, MB2, MZ2) + x*B0(x, MB2&
  &, MZ2))*DBLE(x**INT(-2.D0)))/ (MW2*PI2*SW2)

 amplitudes(5) = (-0.015625D0*CKM13*CKMC13*EL2*DBLE(x**INT(-2.D0))*(MB2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MB2*TB2*A0(MU2)&
  &*DBLE(Yuk3**INT(2.D0)) - 1.D0*MB2*MHp2*TB2*B0(x, MHp2, MU2)*DBLE(Yuk3**INT(2.D0)) + MB2*MU2*TB2*B0(x, MHp2, MU2)*DBLE(Yuk3**IN&
  &T(2.D0)) + MB2*TB2*x*B0(x, MHp2, MU2)*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2) + (0.015625D0*CKM13*CKMC13*EL2*(MB2*TB2*B0(x, &
  &MHp2, MU2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MB2*MHp2*TB2*DB0(x, MHp2, MU2)*DBLE(Yuk3**INT(2.D0)) + MB2*MU2*TB2*DB0(x, MHp2, MU2)*D&
  &BLE(Yuk3**INT(2.D0)) + MB2*TB2*x*DB0(x, MHp2, MU2)*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(6) = (-0.015625D0*CKM23*CKMC23*EL2*DBLE(x**INT(-2.D0))*(-1.D0*MB2*TB2*A0(MC2)*DBLE(Yuk3**INT(2.D0)) + MB2*TB2*A0(MHp2&
  &)*DBLE(Yuk3**INT(2.D0)) + MB2*MC2*TB2*B0(x, MC2, MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MB2*MHp2*TB2*B0(x, MC2, MHp2)*DBLE(Yuk3**I&
  &NT(2.D0)) + MB2*TB2*x*B0(x, MC2, MHp2)*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2) + (0.015625D0*CKM23*CKMC23*EL2*(MB2*TB2*B0(x,&
  & MC2, MHp2)*DBLE(Yuk3**INT(2.D0)) + MB2*MC2*TB2*DB0(x, MC2, MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MB2*MHp2*TB2*DB0(x, MC2, MHp2)*&
  &DBLE(Yuk3**INT(2.D0)) + MB2*TB2*x*DB0(x, MC2, MHp2)*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(7) = (-0.015625D0*CKM33*CKMC33*EL2*DBLE(x**INT(-2.D0))*(MB2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MB2*TB2*A0(MT2)&
  &*DBLE(Yuk3**INT(2.D0)) - 1.D0*MB2*MHp2*TB2*B0(x, MHp2, MT2)*DBLE(Yuk3**INT(2.D0)) + MB2*MT2*TB2*B0(x, MHp2, MT2)*DBLE(Yuk3**IN&
  &T(2.D0)) + MB2*TB2*x*B0(x, MHp2, MT2)*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2) + (0.015625D0*CKM33*CKMC33*EL2*(MB2*TB2*B0(x, &
  &MHp2, MT2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MB2*MHp2*TB2*DB0(x, MHp2, MT2)*DBLE(Yuk3**INT(2.D0)) + MB2*MT2*TB2*DB0(x, MHp2, MT2)*D&
  &BLE(Yuk3**INT(2.D0)) + MB2*TB2*x*DB0(x, MHp2, MT2)*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(8) = (0.015625D0*CKM13*CKMC13*EL2*(MB2*B0(x, MU2, MW2) + MB2*MU2*DB0(x, MU2, MW2) - 1.D0*MB2*MW2*DB0(x, MU2, MW2) + M&
  &B2*x*DB0(x, MU2, MW2)))/ (MW2*PI2*SW2*x) - (0.015625D0*CKM13*CKMC13*EL2*(-1.D0*MB2*A0(MU2) + MB2*A0(MW2) + MB2*MU2*B0(x, MU2, &
  &MW2) - 1.D0*MB2*MW2*B0(x, MU2, MW2) + MB2*x*B0(x, MU2, MW2))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SW2)

 amplitudes(9) = (0.015625D0*CKM23*CKMC23*EL2*(MB2*B0(x, MC2, MW2) + MB2*MC2*DB0(x, MC2, MW2) - 1.D0*MB2*MW2*DB0(x, MC2, MW2) + M&
  &B2*x*DB0(x, MC2, MW2)))/ (MW2*PI2*SW2*x) - (0.015625D0*CKM23*CKMC23*EL2*(-1.D0*MB2*A0(MC2) + MB2*A0(MW2) + MB2*MC2*B0(x, MC2, &
  &MW2) - 1.D0*MB2*MW2*B0(x, MC2, MW2) + MB2*x*B0(x, MC2, MW2))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SW2)

 amplitudes(10) = (0.015625D0*CKM33*CKMC33*EL2*(MB2*B0(x, MT2, MW2) + MB2*MT2*DB0(x, MT2, MW2) - 1.D0*MB2*MW2*DB0(x, MT2, MW2) + &
  &MB2*x*DB0(x, MT2, MW2)))/ (MW2*PI2*SW2*x) - (0.015625D0*CKM33*CKMC33*EL2*(-1.D0*MB2*A0(MT2) + MB2*A0(MW2) + MB2*MT2*B0(x, MT2,&
  & MW2) - 1.D0*MB2*MW2*B0(x, MT2, MW2) + MB2*x*B0(x, MT2, MW2))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SW2)

 amplitudes(11) = (0.006944444444444444D0*EL2*(-1.D0 + B0(x, 0.D0, MB2) + MB2*DB0(x, 0.D0, MB2) + x*DB0(x, 0.D0, MB2)))/(PI2*x) -&
  & (0.006944444444444444D0*EL2*(-1.D0*x - 1.D0*A0(MB2) + MB2*B0(x, 0.D0, MB2) + x*B0(x, 0.D0, MB2))*DBLE(x**INT(-2.D0)))/PI2

 amplitudes(12) = (0.001736111111111111D0*EL2*(-4.D0*DBLE(SW**INT(4.D0)) + 4.D0*B0(x, MB2, MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*MB2*DB&
  &0(x, MB2, MZ2)*DBLE(SW**INT(4.D0)) - 4.D0*MZ2*DB0(x, MB2, MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*x*DB0(x, MB2, MZ2)*DBLE(SW**INT(4.D0&
  &))))/(CW2*PI2*SW2*x) - (0.001736111111111111D0*EL2*(-4.D0*x*DBLE(SW**INT(4.D0)) - 4.D0*A0(MB2)*DBLE(SW**INT(4.D0)) + 4.D0*A0(M&
  &Z2)*DBLE(SW**INT(4.D0)) + 4.D0*MB2*B0(x, MB2, MZ2)*DBLE(SW**INT(4.D0)) - 4.D0*MZ2*B0(x, MB2, MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*x&
  &*B0(x, MB2, MZ2)*DBLE(SW**INT(4.D0)))* DBLE(x**INT(-2.D0)))/(CW2*PI2*SW2)

 amplitudes(13) = 0.D0

 amplitudes(14) = 0.D0

 amplitudes(15) = 0.D0

  totalAmplitude = (0D0,0D0)
 do j=1,15
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 DSelfBBRight = totalAmplitude
end function DSelfBBRight

