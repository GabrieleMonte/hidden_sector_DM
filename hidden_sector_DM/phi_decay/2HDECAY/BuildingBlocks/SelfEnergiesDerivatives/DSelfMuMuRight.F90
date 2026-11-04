double complex function DSelfMuMuRight(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(9)

 amplitudes(1) = (0.0078125D0*EL2*MM2*(B0(x, Mh02, MM2) - 1.D0*Mh02*DB0(x, Mh02, MM2) + MM2*DB0(x, Mh02, MM2) + x*DB0(x, Mh02, MM&
  &2))*DBLE(Yuk4**INT(2.D0)))/ (MW2*PI2*SW2*x) - (0.0078125D0*EL2*MM2*(A0(Mh02) - 1.D0*A0(MM2) - 1.D0*Mh02*B0(x, Mh02, MM2) + MM2&
  &*B0(x, Mh02, MM2) + x*B0(x, Mh02, MM2))* DBLE(x**INT(-2.D0))*DBLE(Yuk4**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(2) = (0.0078125D0*EL2*MM2*(B0(x, MHH2, MM2) - 1.D0*MHH2*DB0(x, MHH2, MM2) + MM2*DB0(x, MHH2, MM2) + x*DB0(x, MHH2, MM&
  &2))*DBLE(Yuk5**INT(2.D0)))/ (MW2*PI2*SW2*x) - (0.0078125D0*EL2*MM2*(A0(MHH2) - 1.D0*A0(MM2) - 1.D0*MHH2*B0(x, MHH2, MM2) + MM2&
  &*B0(x, MHH2, MM2) + x*B0(x, MHH2, MM2))* DBLE(x**INT(-2.D0))*DBLE(Yuk5**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(3) = (0.0078125D0*EL2*MM2*(B0(x, MA02, MM2) - 1.D0*MA02*DB0(x, MA02, MM2) + MM2*DB0(x, MA02, MM2) + x*DB0(x, MA02, MM&
  &2))*DBLE(Yuk6**INT(2.D0)))/ (MW2*PI2*SW2*x) - (0.0078125D0*EL2*MM2*(A0(MA02) - 1.D0*A0(MM2) - 1.D0*MA02*B0(x, MA02, MM2) + MM2&
  &*B0(x, MA02, MM2) + x*B0(x, MA02, MM2))* DBLE(x**INT(-2.D0))*DBLE(Yuk6**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(4) = (0.0078125D0*EL2*MM2*(B0(x, MM2, MZ2) + MM2*DB0(x, MM2, MZ2) - 1.D0*MZ2*DB0(x, MM2, MZ2) + x*DB0(x, MM2, MZ2)))/&
  &(MW2*PI2*SW2*x) - (0.0078125D0*EL2*MM2*(-1.D0*A0(MM2) + A0(MZ2) + MM2*B0(x, MM2, MZ2) - 1.D0*MZ2*B0(x, MM2, MZ2) + x*B0(x, MM2&
  &, MZ2))*DBLE(x**INT(-2.D0)))/ (MW2*PI2*SW2)

 amplitudes(5) = (0.015625D0*EL2*MM2*(B0(x, 0.D0, MHp2) + (-1.D0*MHp2 + x)*DB0(x, 0.D0, MHp2))*DBLE(Yuk6**INT(2.D0)))/(MW2*PI2*SW&
  &2*x) - (0.015625D0*EL2*MM2*(A0(MHp2) + (-1.D0*MHp2 + x)*B0(x, 0.D0, MHp2))*DBLE(x**INT(-2.D0))*DBLE(Yuk6**INT(2.D0)))/(MW2*PI2&
  &*SW2)

 amplitudes(6) = (0.015625D0*EL2*MM2*(B0(x, 0.D0, MW2) + (-1.D0*MW2 + x)*DB0(x, 0.D0, MW2)))/(MW2*PI2*SW2*x) - (0.015625D0*EL2*MM&
  &2*(A0(MW2) + (-1.D0*MW2 + x)*B0(x, 0.D0, MW2))*DBLE(x**INT(-2.D0)))/(MW2*PI2*SW2)

 amplitudes(7) = (0.0625D0*EL2*(-1.D0 + B0(x, 0.D0, MM2) + MM2*DB0(x, 0.D0, MM2) + x*DB0(x, 0.D0, MM2)))/(PI2*x) - (0.0625D0*EL2*&
  &(-1.D0*x - 1.D0*A0(MM2) + MM2*B0(x, 0.D0, MM2) + x*B0(x, 0.D0, MM2))*DBLE(x**INT(-2.D0)))/PI2

 amplitudes(8) = (0.015625D0*EL2*(-4.D0*DBLE(SW**INT(4.D0)) + 4.D0*B0(x, MM2, MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*MM2*DB0(x, MM2, MZ2&
  &)*DBLE(SW**INT(4.D0)) - 4.D0*MZ2*DB0(x, MM2, MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*x*DB0(x, MM2, MZ2)*DBLE(SW**INT(4.D0))))/(CW2*PI2&
  &*SW2*x) - (0.015625D0*EL2*(-4.D0*x*DBLE(SW**INT(4.D0)) - 4.D0*A0(MM2)*DBLE(SW**INT(4.D0)) + 4.D0*A0(MZ2)*DBLE(SW**INT(4.D0)) +&
  & 4.D0*MM2*B0(x, MM2, MZ2)*DBLE(SW**INT(4.D0)) - 4.D0*MZ2*B0(x, MM2, MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*x*B0(x, MM2, MZ2)*DBLE(SW*&
  &*INT(4.D0)))* DBLE(x**INT(-2.D0)))/(CW2*PI2*SW2)

 amplitudes(9) = 0.D0

  totalAmplitude = (0D0,0D0)
 do j=1,9
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 DSelfMuMuRight = totalAmplitude
end function DSelfMuMuRight

