double complex function SelfAZ0Alter(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(19)

 amplitudes(1) = (0.0625D0*EL2*(CW2 - 1.D0*SW2)*A0(MHp2))/(CW*PI2*SW)

 amplitudes(2) = (0.0625D0*EL2*(CW2 - 1.D0*SW2)*A0(MW2))/(CW*PI2*SW)

 amplitudes(3) = (0.125D0*CW*EL2*(-2.D0*MW2 + 3.D0*A0(MW2)))/(PI2*SW)

 amplitudes(4) = (-0.006944444444444444D0*EL2*(-1.D0 + 4.D0*SW2)*(6.D0*ME2 - 1.D0*x - 6.D0*A0(ME2) + 3.D0*(2.D0*ME2 + x)*B0(x, ME&
  &2, ME2)))/(CW*PI2*SW)

 amplitudes(5) = (-0.006944444444444444D0*EL2*(-1.D0 + 4.D0*SW2)*(6.D0*MM2 - 1.D0*x - 6.D0*A0(MM2) + 3.D0*(2.D0*MM2 + x)*B0(x, MM&
  &2, MM2)))/(CW*PI2*SW)

 amplitudes(6) = (-0.006944444444444444D0*EL2*(-1.D0 + 4.D0*SW2)*(6.D0*ML2 - 1.D0*x - 6.D0*A0(ML2) + 3.D0*(2.D0*ML2 + x)*B0(x, ML&
  &2, ML2)))/(CW*PI2*SW)

 amplitudes(7) = (-0.004629629629629629D0*EL2*(-3.D0 + 8.D0*SW2)*(6.D0*MU2 - 1.D0*x - 6.D0*A0(MU2) + 3.D0*(2.D0*MU2 + x)*B0(x, MU&
  &2, MU2)))/(CW*PI2*SW)

 amplitudes(8) = (-0.004629629629629629D0*EL2*(-3.D0 + 8.D0*SW2)*(6.D0*MC2 - 1.D0*x - 6.D0*A0(MC2) + 3.D0*(2.D0*MC2 + x)*B0(x, MC&
  &2, MC2)))/(CW*PI2*SW)

 amplitudes(9) = (-0.004629629629629629D0*EL2*(-3.D0 + 8.D0*SW2)*(6.D0*MT2 - 1.D0*x - 6.D0*A0(MT2) + 3.D0*(2.D0*MT2 + x)*B0(x, MT&
  &2, MT2)))/(CW*PI2*SW)

 amplitudes(10) = (-0.0023148148148148147D0*EL2*(-3.D0 + 4.D0*SW2)*(6.D0*MD2 - 1.D0*x - 6.D0*A0(MD2) + 3.D0*(2.D0*MD2 + x)*B0(x, &
  &MD2, MD2)))/(CW*PI2*SW)

 amplitudes(11) = (-0.0023148148148148147D0*EL2*(-3.D0 + 4.D0*SW2)*(6.D0*MS2 - 1.D0*x - 6.D0*A0(MS2) + 3.D0*(2.D0*MS2 + x)*B0(x, &
  &MS2, MS2)))/(CW*PI2*SW)

 amplitudes(12) = (-0.0023148148148148147D0*EL2*(-3.D0 + 4.D0*SW2)*(6.D0*MB2 - 1.D0*x - 6.D0*A0(MB2) + 3.D0*(2.D0*MB2 + x)*B0(x, &
  &MB2, MB2)))/(CW*PI2*SW)

 amplitudes(13) = (-0.003472222222222222D0*EL2*(CW2 - 1.D0*SW2)*(12.D0*MHp2 - 2.D0*x + 6.D0*A0(MHp2) + 3.D0*(4.D0*MHp2 - 1.D0*x)*&
  &B0(x, MHp2, MHp2)))/(CW*PI2*SW)

 amplitudes(14) = (-0.003472222222222222D0*EL2*(CW2 - 1.D0*SW2)*(12.D0*MW2 - 2.D0*x + 6.D0*A0(MW2) + 3.D0*(4.D0*MW2 - 1.D0*x)*B0(&
  &x, MW2, MW2)))/(CW*PI2*SW)

 amplitudes(15) = (0.001736111111111111D0*CW*EL2*(12.D0*MW2 - 2.D0*x + 6.D0*A0(MW2) + 3.D0*(4.D0*MW2 - 1.D0*x)*B0(x, MW2, MW2)))/&
  &(PI2*SW)

 amplitudes(16) = (0.001736111111111111D0*CW*EL2*(12.D0*MW2 - 2.D0*x + 6.D0*A0(MW2) + 3.D0*(4.D0*MW2 - 1.D0*x)*B0(x, MW2, MW2)))/&
  &(PI2*SW)

 amplitudes(17) = (-0.003472222222222222D0*CW*EL2*(2.D0*(-6.D0*MW2 + x) + 66.D0*A0(MW2) + (96.D0*MW2 + 57.D0*x)*B0(x, MW2, MW2)))&
  &/(PI2*SW)

 amplitudes(18) = (-0.0625D0*EL2*MW2*SW*B0(x, MW2, MW2))/(CW*PI2)

 amplitudes(19) = (-0.0625D0*EL2*MW2*SW*B0(x, MW2, MW2))/(CW*PI2)

  totalAmplitude = (0D0,0D0)
 do j=1,19
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 SelfAZ0Alter = totalAmplitude
end function SelfAZ0Alter

