double complex function DSelfHHh0(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(36)

 amplitudes(1) = 0.D0

 amplitudes(2) = 0.D0

 amplitudes(3) = 0.D0

 amplitudes(4) = 0.D0

 amplitudes(5) = 0.D0

 amplitudes(6) = 0.D0

 amplitudes(7) = (-0.03125D0*EL2*ME2*Yuk4*Yuk5*(-1.D0*B0(x, ME2, ME2) + (4.D0*ME2 - 1.D0*x)*DB0(x, ME2, ME2)))/(MW2*PI2*SW2)

 amplitudes(8) = (-0.03125D0*EL2*MM2*Yuk4*Yuk5*(-1.D0*B0(x, MM2, MM2) + (4.D0*MM2 - 1.D0*x)*DB0(x, MM2, MM2)))/(MW2*PI2*SW2)

 amplitudes(9) = (-0.03125D0*EL2*ML2*Yuk4*Yuk5*(-1.D0*B0(x, ML2, ML2) + (4.D0*ML2 - 1.D0*x)*DB0(x, ML2, ML2)))/(MW2*PI2*SW2)

 amplitudes(10) = (-0.09375D0*CA*EL2*MU2*SA*(-1.D0*B0(x, MU2, MU2) + (4.D0*MU2 - 1.D0*x)*DB0(x, MU2, MU2)))/(MW2*PI2*SB2*SW2)

 amplitudes(11) = (-0.09375D0*CA*EL2*MC2*SA*(-1.D0*B0(x, MC2, MC2) + (4.D0*MC2 - 1.D0*x)*DB0(x, MC2, MC2)))/(MW2*PI2*SB2*SW2)

 amplitudes(12) = (-0.09375D0*CA*EL2*MT2*SA*(-1.D0*B0(x, MT2, MT2) + (4.D0*MT2 - 1.D0*x)*DB0(x, MT2, MT2)))/(MW2*PI2*SB2*SW2)

 amplitudes(13) = (-0.09375D0*EL2*MD2*Yuk1*Yuk2*(-1.D0*B0(x, MD2, MD2) + (4.D0*MD2 - 1.D0*x)*DB0(x, MD2, MD2)))/(MW2*PI2*SW2)

 amplitudes(14) = (-0.09375D0*EL2*MS2*Yuk1*Yuk2*(-1.D0*B0(x, MS2, MS2) + (4.D0*MS2 - 1.D0*x)*DB0(x, MS2, MS2)))/(MW2*PI2*SW2)

 amplitudes(15) = (-0.09375D0*EL2*MB2*Yuk1*Yuk2*(-1.D0*B0(x, MB2, MB2) + (4.D0*MB2 - 1.D0*x)*DB0(x, MB2, MB2)))/(MW2*PI2*SW2)

 amplitudes(16) = (-0.0234375D0*CBA*(-2.D0*CAB*EL2*Mh02 - 1.D0*EL2*Mh02*S2A*SBA + 4.D0*CAB*CBA2*Lambda5*MW2*SW2)* (EL2*(2.D0*Mh02&
  & + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*DB0(x, Mh02, Mh02))/(EL2*MW2*PI2*S2B2*SW2)

 amplitudes(17) = (0.0234375D0*SBA*(EL2*(Mh02 + 2.D0*MHH2)*S2A - 2.D0*Lambda5*MW2*(3.D0*S2A + S2B)*SW2)*(CBA*EL2*MHH2*S2A - 2.D0*&
  &EL2*MHH2*SAB + 4.D0*Lambda5*MW2*SAB*SBA2*SW2)*DB0(x, MHH2, MHH2))/(EL2*MW2*PI2*S2B2*SW2)

 amplitudes(18) = (-0.015625D0*CBA*SBA*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)* (EL2*(Mh02 + 2.D0*M&
  &HH2)*S2A - 2.D0*Lambda5*MW2*(3.D0*S2A + S2B)*SW2)*DB0(x, Mh02, MHH2))/(EL2*MW2*PI2*S2B2*SW2)

 amplitudes(19) = (0.0078125D0*(EL2*(2.D0*MA02 - 1.D0*Mh02)*S2B*SBA + 2.D0*CAB*(EL2*Mh02 - 2.D0*Lambda5*MW2*SW2))* (CBA*EL2*(2.D0&
  &*MA02 - 1.D0*MHH2)*S2B + 2.D0*SAB*(EL2*MHH2 - 2.D0*Lambda5*MW2*SW2))*DB0(x, MA02, MA02))/(EL2*MW2*PI2*S2B2*SW2)

 amplitudes(20) = (0.0078125D0*CBA*EL2*Mh02*MHH2*SBA*DB0(x, MZ2, MZ2))/(MW2*PI2*SW2)

 amplitudes(21) = (-0.015625D0*CBA*EL2*(MA02 - 1.D0*Mh02)*(MA02 - 1.D0*MHH2)*SBA*DB0(x, MA02, MZ2))/(MW2*PI2*SW2)

 amplitudes(22) = (0.015625D0*(-1.D0*EL2*(Mh02 - 2.D0*MHp2)*S2B*SBA + 2.D0*CAB*(EL2*Mh02 - 2.D0*Lambda5*MW2*SW2))* (-1.D0*CBA*EL2&
  &*(MHH2 - 2.D0*MHp2)*S2B + 2.D0*SAB*(EL2*MHH2 - 2.D0*Lambda5*MW2*SW2))*DB0(x, MHp2, MHp2))/(EL2*MW2*PI2*S2B2*SW2)

 amplitudes(23) = (0.015625D0*CBA*EL2*Mh02*MHH2*SBA*DB0(x, MW2, MW2))/(MW2*PI2*SW2)

 amplitudes(24) = (0.015625D0*CBA*EL2*(Mh02 - 1.D0*MHp2)*(-1.D0*MHH2 + MHp2)*SBA*DB0(x, MHp2, MW2))/(MW2*PI2*SW2)

 amplitudes(25) = (0.015625D0*CBA*EL2*(Mh02 - 1.D0*MHp2)*(-1.D0*MHH2 + MHp2)*SBA*DB0(x, MHp2, MW2))/(MW2*PI2*SW2)

 amplitudes(26) = (-0.015625D0*CBA*EL2*MW2*SBA*DB0(x, MZ2, MZ2)*DBLE(CW**INT(-4.D0)))/(PI2*SW2)

 amplitudes(27) = (-0.015625D0*CBA*EL2*MW2*SBA*DB0(x, MW2, MW2))/(PI2*SW2)

 amplitudes(28) = (-0.015625D0*CBA*EL2*MW2*SBA*DB0(x, MW2, MW2))/(PI2*SW2)

 amplitudes(29) = (0.125D0*CBA*EL2*MW2*SBA*DB0(x, MZ2, MZ2)*DBLE(CW**INT(-4.D0)))/(PI2*SW2)

 amplitudes(30) = (0.25D0*CBA*EL2*MW2*SBA*DB0(x, MW2, MW2))/(PI2*SW2)

 amplitudes(31) = (-0.015625D0*CBA*EL2*SBA*(-2.D0*B0(x, MA02, MZ2) + (-2.D0*MA02 + MZ2 - 2.D0*x)*DB0(x, MA02, MZ2)))/(CW2*PI2*SW2&
  &)

 amplitudes(32) = (-0.015625D0*CBA*EL2*SBA*(2.D0*B0(x, MZ2, MZ2) + (MZ2 + 2.D0*x)*DB0(x, MZ2, MZ2)))/(CW2*PI2*SW2)

 amplitudes(33) = (-0.015625D0*CBA*EL2*SBA*(-2.D0*B0(x, MHp2, MW2) + (-2.D0*MHp2 + MW2 - 2.D0*x)*DB0(x, MHp2, MW2)))/(PI2*SW2)

 amplitudes(34) = (-0.015625D0*CBA*EL2*SBA*(-2.D0*B0(x, MHp2, MW2) + (-2.D0*MHp2 + MW2 - 2.D0*x)*DB0(x, MHp2, MW2)))/(PI2*SW2)

 amplitudes(35) = (-0.015625D0*CBA*EL2*SBA*(2.D0*B0(x, MW2, MW2) + (MW2 + 2.D0*x)*DB0(x, MW2, MW2)))/(PI2*SW2)

 amplitudes(36) = (-0.015625D0*CBA*EL2*SBA*(2.D0*B0(x, MW2, MW2) + (MW2 + 2.D0*x)*DB0(x, MW2, MW2)))/(PI2*SW2)

  totalAmplitude = (0D0,0D0)
 do j=1,36
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 DSelfHHh0 = totalAmplitude
end function DSelfHHh0

