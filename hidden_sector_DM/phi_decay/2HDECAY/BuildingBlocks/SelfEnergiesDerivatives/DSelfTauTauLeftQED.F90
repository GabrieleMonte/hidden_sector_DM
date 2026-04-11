double complex function DSelfTauTauLeftQED(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(1)

 amplitudes(1) = (0.0625D0*EL2*(-1.D0 + B0(x, 0.D0, ML2) + ML2*DB0(x, 0.D0, ML2) + x*DB0(x, 0.D0, ML2)))/(PI2*x) - (0.0625D0*EL2*&
  &(-1.D0*x - 1.D0*A0(ML2) + ML2*B0(x, 0.D0, ML2) + x*B0(x, 0.D0, ML2))*DBLE(x**INT(-2.D0)))/PI2

  totalAmplitude = (0D0,0D0)
 do j=1,1
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 DSelfTauTauLeftQED = totalAmplitude
end function DSelfTauTauLeftQED

