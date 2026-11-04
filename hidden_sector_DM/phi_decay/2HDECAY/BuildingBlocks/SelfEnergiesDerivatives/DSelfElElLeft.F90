double complex function DSelfElElLeft(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(9)

 amplitudes(1) = (0.0078125D0*EL2*ME2*(B0(x, ME2, Mh02) + ME2*DB0(x, ME2, Mh02) - 1.D0*Mh02*DB0(x, ME2, Mh02) + x*DB0(x, ME2, Mh0&
  &2))*DBLE(Yuk4**INT(2.D0)))/ (MW2*PI2*SW2*x) - (0.0078125D0*EL2*ME2*(-1.D0*A0(ME2) + A0(Mh02) + ME2*B0(x, ME2, Mh02) - 1.D0*Mh0&
  &2*B0(x, ME2, Mh02) + x*B0(x, ME2, Mh02))* DBLE(x**INT(-2.D0))*DBLE(Yuk4**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(2) = (0.0078125D0*EL2*ME2*(B0(x, ME2, MHH2) + ME2*DB0(x, ME2, MHH2) - 1.D0*MHH2*DB0(x, ME2, MHH2) + x*DB0(x, ME2, MHH&
  &2))*DBLE(Yuk5**INT(2.D0)))/ (MW2*PI2*SW2*x) - (0.0078125D0*EL2*ME2*(-1.D0*A0(ME2) + A0(MHH2) + ME2*B0(x, ME2, MHH2) - 1.D0*MHH&
  &2*B0(x, ME2, MHH2) + x*B0(x, ME2, MHH2))* DBLE(x**INT(-2.D0))*DBLE(Yuk5**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(3) = (0.0078125D0*EL2*ME2*(B0(x, MA02, ME2) - 1.D0*MA02*DB0(x, MA02, ME2) + ME2*DB0(x, MA02, ME2) + x*DB0(x, MA02, ME&
  &2))*DBLE(Yuk6**INT(2.D0)))/ (MW2*PI2*SW2*x) - (0.0078125D0*EL2*ME2*(A0(MA02) - 1.D0*A0(ME2) - 1.D0*MA02*B0(x, MA02, ME2) + ME2&
  &*B0(x, MA02, ME2) + x*B0(x, MA02, ME2))* DBLE(x**INT(-2.D0))*DBLE(Yuk6**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(4) = (0.0078125D0*EL2*ME2*(B0(x, ME2, MZ2) + ME2*DB0(x, ME2, MZ2) - 1.D0*MZ2*DB0(x, ME2, MZ2) + x*DB0(x, ME2, MZ2)))/&
  &(MW2*PI2*SW2*x) - (0.0078125D0*EL2*ME2*(-1.D0*A0(ME2) + A0(MZ2) + ME2*B0(x, ME2, MZ2) - 1.D0*MZ2*B0(x, ME2, MZ2) + x*B0(x, ME2&
  &, MZ2))*DBLE(x**INT(-2.D0)))/ (MW2*PI2*SW2)

 amplitudes(5) = 0.D0

 amplitudes(6) = 0.D0

 amplitudes(7) = (0.0625D0*EL2*(-1.D0 + B0(x, 0.D0, ME2) + ME2*DB0(x, 0.D0, ME2) + x*DB0(x, 0.D0, ME2)))/(PI2*x) - (0.0625D0*EL2*&
  &(-1.D0*x - 1.D0*A0(ME2) + ME2*B0(x, 0.D0, ME2) + x*B0(x, 0.D0, ME2))*DBLE(x**INT(-2.D0)))/PI2

 amplitudes(8) = (0.015625D0*EL2*(-1.D0 + 4.D0*SW2 + B0(x, ME2, MZ2) - 4.D0*SW2*B0(x, ME2, MZ2) + ME2*DB0(x, ME2, MZ2) - 1.D0*MZ2&
  &*DB0(x, ME2, MZ2) - 4.D0*ME2*SW2*DB0(x, ME2, MZ2) + 4.D0*MZ2*SW2*DB0(x, ME2, MZ2) + x*DB0(x, ME2, MZ2) - 4.D0*SW2*x*DB0(x, ME2&
  &, MZ2) - 4.D0*DBLE(SW**INT(4.D0)) + 4.D0*B0(x, ME2, MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*ME2*DB0(x, ME2, MZ2)*DBLE(SW**INT(4.D0)) -&
  & 4.D0*MZ2*DB0(x, ME2, MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*x*DB0(x, ME2, MZ2)*DBLE(SW**INT(4.D0))))/(CW2*PI2*SW2*x) - (0.015625D0*E&
  &L2*(-1.D0*x + 4.D0*SW2*x - 1.D0*A0(ME2) + 4.D0*SW2*A0(ME2) + A0(MZ2) - 4.D0*SW2*A0(MZ2) + ME2*B0(x, ME2, MZ2) - 1.D0*MZ2*B0(x,&
  & ME2, MZ2) - 4.D0*ME2*SW2*B0(x, ME2, MZ2) + 4.D0*MZ2*SW2*B0(x, ME2, MZ2) + x*B0(x, ME2, MZ2) - 4.D0*SW2*x*B0(x, ME2, MZ2) - 4.&
  &D0*x*DBLE(SW**INT(4.D0)) - 4.D0*A0(ME2)*DBLE(SW**INT(4.D0)) + 4.D0*A0(MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*ME2*B0(x, ME2, MZ2)*DBLE&
  &(SW**INT(4.D0)) - 4.D0*MZ2*B0(x, ME2, MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*x*B0(x, ME2, MZ2)*DBLE(SW**INT(4.D0)))*DBLE(x**INT(-2.D0&
  &)))/(CW2*PI2*SW2)

 amplitudes(9) = (0.03125D0*EL2*(-1.D0 + B0(x, 0.D0, MW2) + (-1.D0*MW2 + x)*DB0(x, 0.D0, MW2)))/(PI2*SW2*x) - (0.03125D0*EL2*(-1.&
  &D0*x + A0(MW2) + (-1.D0*MW2 + x)*B0(x, 0.D0, MW2))*DBLE(x**INT(-2.D0)))/(PI2*SW2)

  totalAmplitude = (0D0,0D0)
 do j=1,9
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 DSelfElElLeft = totalAmplitude
end function DSelfElElLeft

