double complex function SelfNeuTNeuTLeftUsual(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(4)

 amplitudes(1) = (-0.03125D0*EL2*ML2*B1(x, ML2, MHp2)*DBLE(Yuk6**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(2) = (-0.03125D0*EL2*ML2*B1(x, ML2, MW2))/(MW2*PI2*SW2)

 amplitudes(3) = (-0.015625D0*EL2*(1.D0 + 2.D0*B1(x, 0.D0, MZ2)))/(CW2*PI2*SW2)

 amplitudes(4) = (-0.03125D0*EL2*(1.D0 + 2.D0*B1(x, ML2, MW2)))/(PI2*SW2)

  totalAmplitude = (0D0,0D0)
 do j=1,4
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 SelfNeuTNeuTLeftUsual = totalAmplitude
end function SelfNeuTNeuTLeftUsual

