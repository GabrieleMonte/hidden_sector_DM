double complex function Selfh0h0Alter(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(78)

 amplitudes(1) = (0.0234375D0*(4.D0*CAB2*EL2*Mh02 + CBA2*EL2*MHH2*S2A2 + 4.D0*CAB*EL2*Mh02*S2A*SBA + EL2*Mh02*S2A2*SBA2 - 2.D0*C2&
  &A2*Lambda5*MW2*SW2 - 4.D0*C2A*C2B*Lambda5*MW2*SW2 - 2.D0*C2B2*Lambda5*MW2*SW2)*A0(Mh02))/(MW2*PI2*S2B2*SW2)

 amplitudes(2) = (0.0078125D0*(3.D0*CBA2*EL2*Mh02*S2A2 + EL2*S2A*(-1.D0*Mh02*S2B + MHH2*(S2B + 3.D0*S2A*SBA2)) + 2.D0*Lambda5*MW2&
  &*(-3.D0*S2A2 + S2B2)*SW2)*A0(MHH2))/ (MW2*PI2*S2B2*SW2)

 amplitudes(3) = (-0.0078125D0*(-2.D0*CBA*EL2*MHH2*S2A*SAB + 2.D0*CAB*EL2*Mh02*(-1.D0*S2A + S2B)*SBA + EL2*S2B*(Mh02*S2A - 2.D0*M&
  &A02*S2B)*SBA2 + CAB2*(-4.D0*EL2*Mh02 + 4.D0*Lambda5*MW2*SW2) + CBA2*(EL2*MHH2*S2A*S2B + 4.D0*C2B2*Lambda5*MW2*SW2))*A0(MA02))/&
  &(MW2*PI2*S2B2*SW2)

 amplitudes(4) = (-0.0078125D0*(-1.D0*EL2*Mh02*S2B + CBA2*(EL2*(Mh02*S2A - 1.D0*MHH2*S2A - 2.D0*MA02*S2B) + 4.D0*Lambda5*MW2*S2B*&
  &SW2))*A0(MZ2))/(MW2*PI2*S2B*SW2)

 amplitudes(5) = (-0.015625D0*(-2.D0*CBA*EL2*MHH2*S2A*SAB + 2.D0*CAB*EL2*Mh02*(-1.D0*S2A + S2B)*SBA + EL2*S2B*(Mh02*S2A - 2.D0*MH&
  &p2*S2B)*SBA2 + CAB2*(-4.D0*EL2*Mh02 + 4.D0*Lambda5*MW2*SW2) + CBA2*(EL2*MHH2*S2A*S2B + 4.D0*C2B2*Lambda5*MW2*SW2))*A0(MHp2))/(&
  &MW2*PI2*S2B2*SW2)

 amplitudes(6) = (-0.015625D0*(-1.D0*EL2*Mh02*S2B + CBA2*(EL2*(Mh02*S2A - 1.D0*MHH2*S2A - 2.D0*MHp2*S2B) + 4.D0*Lambda5*MW2*S2B*S&
  &W2))*A0(MW2))/(MW2*PI2*S2B*SW2)

 amplitudes(7) = (-0.03125D0*EL2*(MZ2 - 2.D0*A0(MZ2)))/(CW2*PI2*SW2)

 amplitudes(8) = (-0.0625D0*EL2*(MW2 - 2.D0*A0(MW2)))/(PI2*SW2)

 amplitudes(9) = (0.1875D0*ME2*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*Yuk4*A0(ME2))/(Mh02*MW2*PI2&
  &*S2B*SW2)

 amplitudes(10) = (0.1875D0*MM2*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*Yuk4*A0(MM2))/(Mh02*MW2*PI&
  &2*S2B*SW2)

 amplitudes(11) = (0.1875D0*ML2*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*Yuk4*A0(ML2))/(Mh02*MW2*PI&
  &2*S2B*SW2)

 amplitudes(12) = (0.5625D0*CA*MU2*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*A0(MU2))/(Mh02*MW2*PI2*&
  &S2B*SB*SW2)

 amplitudes(13) = (0.5625D0*CA*MC2*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*A0(MC2))/(Mh02*MW2*PI2*&
  &S2B*SB*SW2)

 amplitudes(14) = (0.5625D0*CA*MT2*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*A0(MT2))/(Mh02*MW2*PI2*&
  &S2B*SB*SW2)

 amplitudes(15) = (0.5625D0*MD2*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*Yuk1*A0(MD2))/(Mh02*MW2*PI&
  &2*S2B*SW2)

 amplitudes(16) = (0.5625D0*MS2*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*Yuk1*A0(MS2))/(Mh02*MW2*PI&
  &2*S2B*SW2)

 amplitudes(17) = (0.5625D0*MB2*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*Yuk1*A0(MB2))/(Mh02*MW2*PI&
  &2*S2B*SW2)

 amplitudes(18) = (0.0625D0*CBA*ME2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*Yuk5*A0(ME2))/(MHH2*MW2&
  &*PI2*S2B*SW2)

 amplitudes(19) = (0.0625D0*CBA*MM2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*Yuk5*A0(MM2))/(MHH2*MW2&
  &*PI2*S2B*SW2)

 amplitudes(20) = (0.0625D0*CBA*ML2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*Yuk5*A0(ML2))/(MHH2*MW2&
  &*PI2*S2B*SW2)

 amplitudes(21) = (0.1875D0*CBA*MU2*SA*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*A0(MU2))/(MHH2*MW2*P&
  &I2*S2B*SB*SW2)

 amplitudes(22) = (0.1875D0*CBA*MC2*SA*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*A0(MC2))/(MHH2*MW2*P&
  &I2*S2B*SB*SW2)

 amplitudes(23) = (0.1875D0*CBA*MT2*SA*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*A0(MT2))/(MHH2*MW2*P&
  &I2*S2B*SB*SW2)

 amplitudes(24) = (0.1875D0*CBA*MD2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*Yuk2*A0(MD2))/(MHH2*MW2&
  &*PI2*S2B*SW2)

 amplitudes(25) = (0.1875D0*CBA*MS2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*Yuk2*A0(MS2))/(MHH2*MW2&
  &*PI2*S2B*SW2)

 amplitudes(26) = (0.1875D0*CBA*MB2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*Yuk2*A0(MB2))/(MHH2*MW2&
  &*PI2*S2B*SW2)

 amplitudes(27) = (-0.0703125D0*A0(Mh02)*DBLE((2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)**INT(2.D0)))&
  &/(EL2*Mh02*MW2*PI2*S2B2*SW2)

 amplitudes(28) = (0.0234375D0*SBA*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)* (EL2*(Mh02 + 2.D0*MHH2&
  &)*S2A - 2.D0*Lambda5*MW2*(3.D0*S2A + S2B)*SW2)*A0(MHH2))/(EL2*Mh02*MW2*PI2*S2B2*SW2)

 amplitudes(29) = (0.0234375D0*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)* (EL2*(-2.D0*MA02 + Mh02)*S&
  &2B*SBA + CAB*(-2.D0*EL2*Mh02 + 4.D0*Lambda5*MW2*SW2))*A0(MA02))/(EL2*Mh02*MW2*PI2*S2B2*SW2)

 amplitudes(30) = (-0.0234375D0*SBA*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*A0(MZ2))/(MW2*PI2*S2B*&
  &SW2)

 amplitudes(31) = (0.046875D0*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)* (EL2*(Mh02 - 2.D0*MHp2)*S2B&
  &*SBA + CAB*(-2.D0*EL2*Mh02 + 4.D0*Lambda5*MW2*SW2))*A0(MHp2))/(EL2*Mh02*MW2*PI2*S2B2*SW2)

 amplitudes(32) = (-0.046875D0*SBA*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*A0(MW2))/(MW2*PI2*S2B*S&
  &W2)

 amplitudes(33) = (-0.0078125D0*CBA2*A0(Mh02)*DBLE((EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)**INT(2.D&
  &0)))/(EL2*MHH2*MW2*PI2*S2B2*SW2)

 amplitudes(34) = (0.0234375D0*CBA*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*(CBA*EL2*MHH2*S2A - 2.D0&
  &*EL2*MHH2*SAB + 4.D0*Lambda5*MW2*SAB*SBA2*SW2)*A0(MHH2))/(EL2*MHH2*MW2*PI2*S2B2*SW2)

 amplitudes(35) = (0.0078125D0*CBA*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*(CBA*EL2*(-2.D0*MA02 + M&
  &HH2)*S2B - 2.D0*EL2*MHH2*SAB + 4.D0*Lambda5*MW2*SAB*SW2)*A0(MA02))/(EL2*MHH2*MW2*PI2*S2B2*SW2)

 amplitudes(36) = (-0.0078125D0*CBA2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*A0(MZ2))/(MW2*PI2*S2B*&
  &SW2)

 amplitudes(37) = (0.015625D0*CBA*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*(CBA*EL2*(MHH2 - 2.D0*MHp&
  &2)*S2B - 2.D0*EL2*MHH2*SAB + 4.D0*Lambda5*MW2*SAB*SW2)*A0(MHp2))/(EL2*MHH2*MW2*PI2*S2B2*SW2)

 amplitudes(38) = (-0.015625D0*CBA2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*A0(MW2))/(MW2*PI2*S2B*S&
  &W2)

 amplitudes(39) = (0.046875D0*SBA*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*A0(MZ2))/(CW2*Mh02*PI2*S&
  &2B*SW2)

 amplitudes(40) = (0.046875D0*SBA*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*A0(MW2))/(Mh02*PI2*S2B*S&
  &W2)

 amplitudes(41) = (0.046875D0*SBA*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*A0(MW2))/(Mh02*PI2*S2B*S&
  &W2)

 amplitudes(42) = (0.015625D0*CBA2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*A0(MZ2))/(CW2*MHH2*PI2*S&
  &2B*SW2)

 amplitudes(43) = (0.015625D0*CBA2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*A0(MW2))/(MHH2*PI2*S2B*S&
  &W2)

 amplitudes(44) = (0.015625D0*CBA2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*A0(MW2))/(MHH2*PI2*S2B*S&
  &W2)

 amplitudes(45) = (0.09375D0*SBA*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*(MZ2 - 2.D0*A0(MZ2)))/(CW&
  &2*Mh02*PI2*S2B*SW2)

 amplitudes(46) = (0.1875D0*SBA*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*(MW2 - 2.D0*A0(MW2)))/(Mh0&
  &2*PI2*S2B*SW2)

 amplitudes(47) = (0.03125D0*CBA2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*(MZ2 - 2.D0*A0(MZ2)))/(CW&
  &2*MHH2*PI2*S2B*SW2)

 amplitudes(48) = (0.0625D0*CBA2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*(MW2 - 2.D0*A0(MW2)))/(MHH&
  &2*PI2*S2B*SW2)

 amplitudes(49) = (-0.03125D0*EL2*ME2*(2.D0*A0(ME2) + (4.D0*ME2 - 1.D0*x)*B0(x, ME2, ME2))*DBLE(Yuk4**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(50) = (-0.03125D0*EL2*MM2*(2.D0*A0(MM2) + (4.D0*MM2 - 1.D0*x)*B0(x, MM2, MM2))*DBLE(Yuk4**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(51) = (-0.03125D0*EL2*ML2*(2.D0*A0(ML2) + (4.D0*ML2 - 1.D0*x)*B0(x, ML2, ML2))*DBLE(Yuk4**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(52) = (-0.09375D0*CA2*EL2*MU2*(2.D0*A0(MU2) + (4.D0*MU2 - 1.D0*x)*B0(x, MU2, MU2)))/(MW2*PI2*SB2*SW2)

 amplitudes(53) = (-0.09375D0*CA2*EL2*MC2*(2.D0*A0(MC2) + (4.D0*MC2 - 1.D0*x)*B0(x, MC2, MC2)))/(MW2*PI2*SB2*SW2)

 amplitudes(54) = (-0.09375D0*CA2*EL2*MT2*(2.D0*A0(MT2) + (4.D0*MT2 - 1.D0*x)*B0(x, MT2, MT2)))/(MW2*PI2*SB2*SW2)

 amplitudes(55) = (-0.09375D0*EL2*MD2*(2.D0*A0(MD2) + (4.D0*MD2 - 1.D0*x)*B0(x, MD2, MD2))*DBLE(Yuk1**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(56) = (-0.09375D0*EL2*MS2*(2.D0*A0(MS2) + (4.D0*MS2 - 1.D0*x)*B0(x, MS2, MS2))*DBLE(Yuk1**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(57) = (-0.09375D0*EL2*MB2*(2.D0*A0(MB2) + (4.D0*MB2 - 1.D0*x)*B0(x, MB2, MB2))*DBLE(Yuk1**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(58) = (0.0703125D0*B0(x, Mh02, Mh02)*DBLE((2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)**INT&
  &(2.D0)))/(EL2*MW2*PI2*S2B2*SW2)

 amplitudes(59) = (0.0078125D0*SBA2*B0(x, MHH2, MHH2)*DBLE((EL2*(Mh02 + 2.D0*MHH2)*S2A - 2.D0*Lambda5*MW2*(3.D0*S2A + S2B)*SW2)**&
  &INT(2.D0)))/(EL2*MW2*PI2*S2B2*SW2)

 amplitudes(60) = (0.015625D0*CBA2*B0(x, Mh02, MHH2)*DBLE((EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)**&
  &INT(2.D0)))/(EL2*MW2*PI2*S2B2*SW2)

 amplitudes(61) = (0.0078125D0*B0(x, MA02, MA02)*DBLE((EL2*(-2.D0*MA02 + Mh02)*S2B*SBA + CAB*(-2.D0*EL2*Mh02 + 4.D0*Lambda5*MW2*S&
  &W2))**INT(2.D0)))/ (EL2*MW2*PI2*S2B2*SW2)

 amplitudes(62) = (0.0078125D0*EL2*SBA2*B0(x, MZ2, MZ2)*DBLE(Mh0**INT(4.D0)))/(MW2*PI2*SW2)

 amplitudes(63) = (0.015625D0*CBA2*EL2*B0(x, MA02, MZ2)*DBLE((MA02 - 1.D0*Mh02)**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(64) = (0.015625D0*B0(x, MHp2, MHp2)*DBLE((EL2*(Mh02 - 2.D0*MHp2)*S2B*SBA + CAB*(-2.D0*EL2*Mh02 + 4.D0*Lambda5*MW2*SW2&
  &))**INT(2.D0)))/ (EL2*MW2*PI2*S2B2*SW2)

 amplitudes(65) = (0.015625D0*EL2*SBA2*B0(x, MW2, MW2)*DBLE(Mh0**INT(4.D0)))/(MW2*PI2*SW2)

 amplitudes(66) = (0.015625D0*CBA2*EL2*B0(x, MHp2, MW2)*DBLE((Mh02 - 1.D0*MHp2)**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(67) = (0.015625D0*CBA2*EL2*B0(x, MHp2, MW2)*DBLE((Mh02 - 1.D0*MHp2)**INT(2.D0)))/(MW2*PI2*SW2)

 amplitudes(68) = (-0.015625D0*EL2*MW2*SBA2*B0(x, MZ2, MZ2)*DBLE(CW**INT(-4.D0)))/(PI2*SW2)

 amplitudes(69) = (-0.015625D0*EL2*MW2*SBA2*B0(x, MW2, MW2))/(PI2*SW2)

 amplitudes(70) = (-0.015625D0*EL2*MW2*SBA2*B0(x, MW2, MW2))/(PI2*SW2)

 amplitudes(71) = (0.0625D0*EL2*MW2*SBA2*(-1.D0 + 2.D0*B0(x, MZ2, MZ2))*DBLE(CW**INT(-4.D0)))/(PI2*SW2)

 amplitudes(72) = (0.125D0*EL2*MW2*SBA2*(-1.D0 + 2.D0*B0(x, MW2, MW2)))/(PI2*SW2)

 amplitudes(73) = (0.015625D0*CBA2*EL2*(A0(MA02) - 2.D0*A0(MZ2) + (-2.D0*MA02 + MZ2 - 2.D0*x)*B0(x, MA02, MZ2)))/(CW2*PI2*SW2)

 amplitudes(74) = (-0.015625D0*EL2*SBA2*(A0(MZ2) + (MZ2 + 2.D0*x)*B0(x, MZ2, MZ2)))/(CW2*PI2*SW2)

 amplitudes(75) = (0.015625D0*CBA2*EL2*(A0(MHp2) - 2.D0*A0(MW2) + (-2.D0*MHp2 + MW2 - 2.D0*x)*B0(x, MHp2, MW2)))/(PI2*SW2)

 amplitudes(76) = (0.015625D0*CBA2*EL2*(A0(MHp2) - 2.D0*A0(MW2) + (-2.D0*MHp2 + MW2 - 2.D0*x)*B0(x, MHp2, MW2)))/(PI2*SW2)

 amplitudes(77) = (-0.015625D0*EL2*SBA2*(A0(MW2) + (MW2 + 2.D0*x)*B0(x, MW2, MW2)))/(PI2*SW2)

 amplitudes(78) = (-0.015625D0*EL2*SBA2*(A0(MW2) + (MW2 + 2.D0*x)*B0(x, MW2, MW2)))/(PI2*SW2)

  totalAmplitude = (0D0,0D0)
 do j=1,78
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 Selfh0h0Alter = totalAmplitude
end function Selfh0h0Alter

