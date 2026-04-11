double complex function SelfMuMuRightUsual(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(9)

 amplitudes(1) = (0.0078125D0*EL2*MM2*(A0(Mh02) - 1.D0*A0(MM2) - 1.D0*Mh02*B0(x, Mh02, MM2) + MM2*B0(x, Mh02, MM2) + x*B0(x, Mh02&
  &, MM2))*DBLE(Yuk4**INT(2.D0)))/ (MW2*PI2*SW2*x)

 amplitudes(2) = (0.0078125D0*EL2*MM2*(A0(MHH2) - 1.D0*A0(MM2) - 1.D0*MHH2*B0(x, MHH2, MM2) + MM2*B0(x, MHH2, MM2) + x*B0(x, MHH2&
  &, MM2))*DBLE(Yuk5**INT(2.D0)))/ (MW2*PI2*SW2*x)

 amplitudes(3) = (0.0078125D0*EL2*MM2*(A0(MA02) - 1.D0*A0(MM2) - 1.D0*MA02*B0(x, MA02, MM2) + MM2*B0(x, MA02, MM2) + x*B0(x, MA02&
  &, MM2))*DBLE(Yuk6**INT(2.D0)))/ (MW2*PI2*SW2*x)

 amplitudes(4) = (0.0078125D0*EL2*MM2*(-1.D0*A0(MM2) + A0(MZ2) + MM2*B0(x, MM2, MZ2) - 1.D0*MZ2*B0(x, MM2, MZ2) + x*B0(x, MM2, MZ&
  &2)))/(MW2*PI2*SW2*x)

 amplitudes(5) = (0.015625D0*EL2*MM2*(A0(MHp2) + (-1.D0*MHp2 + x)*B0(x, 0.D0, MHp2))*DBLE(Yuk6**INT(2.D0)))/(MW2*PI2*SW2*x)

 amplitudes(6) = (0.015625D0*EL2*MM2*(A0(MW2) + (-1.D0*MW2 + x)*B0(x, 0.D0, MW2)))/(MW2*PI2*SW2*x)

 amplitudes(7) = (0.0625D0*EL2*(-1.D0*x - 1.D0*A0(MM2) + MM2*B0(x, 0.D0, MM2) + x*B0(x, 0.D0, MM2)))/(PI2*x)

 amplitudes(8) = (0.015625D0*EL2*(-4.D0*x*DBLE(SW**INT(4.D0)) - 4.D0*A0(MM2)*DBLE(SW**INT(4.D0)) + 4.D0*A0(MZ2)*DBLE(SW**INT(4.D0&
  &)) + 4.D0*MM2*B0(x, MM2, MZ2)*DBLE(SW**INT(4.D0)) - 4.D0*MZ2*B0(x, MM2, MZ2)*DBLE(SW**INT(4.D0)) + 4.D0*x*B0(x, MM2, MZ2)*DBLE&
  &(SW**INT(4.D0))))/ (CW2*PI2*SW2*x)

 amplitudes(9) = 0.D0

  totalAmplitude = (0D0,0D0)
 do j=1,9
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 SelfMuMuRightUsual = totalAmplitude
end function SelfMuMuRightUsual

