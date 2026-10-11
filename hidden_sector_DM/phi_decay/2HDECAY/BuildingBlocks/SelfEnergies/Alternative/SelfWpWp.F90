double complex function SelfWpWpAlter(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(108)

 amplitudes(1) = (0.015625D0*EL2*A0(Mh02))/(PI2*SW2)

 amplitudes(2) = (0.015625D0*EL2*A0(MHH2))/(PI2*SW2)

 amplitudes(3) = (0.015625D0*EL2*A0(MA02))/(PI2*SW2)

 amplitudes(4) = (0.015625D0*EL2*A0(MZ2))/(PI2*SW2)

 amplitudes(5) = (0.03125D0*EL2*A0(MHp2))/(PI2*SW2)

 amplitudes(6) = (0.03125D0*EL2*A0(MW2))/(PI2*SW2)

 amplitudes(7) = 0.D0

 amplitudes(8) = (0.0625D0*CW2*EL2*(-2.D0*MZ2 + 3.D0*A0(MZ2)))/(PI2*SW2)

 amplitudes(9) = (0.0625D0*EL2*(-2.D0*MW2 + 3.D0*A0(MW2)))/(PI2*SW2)

 amplitudes(10) = (0.125D0*EL2*ME2*SBA*Yuk4*A0(ME2))/(Mh02*PI2*SW2)

 amplitudes(11) = (0.125D0*EL2*MM2*SBA*Yuk4*A0(MM2))/(Mh02*PI2*SW2)

 amplitudes(12) = (0.125D0*EL2*ML2*SBA*Yuk4*A0(ML2))/(Mh02*PI2*SW2)

 amplitudes(13) = (0.375D0*CA*EL2*MU2*SBA*A0(MU2))/(Mh02*PI2*SB*SW2)

 amplitudes(14) = (0.375D0*CA*EL2*MC2*SBA*A0(MC2))/(Mh02*PI2*SB*SW2)

 amplitudes(15) = (0.375D0*CA*EL2*MT2*SBA*A0(MT2))/(Mh02*PI2*SB*SW2)

 amplitudes(16) = (0.375D0*EL2*MD2*SBA*Yuk1*A0(MD2))/(Mh02*PI2*SW2)

 amplitudes(17) = (0.375D0*EL2*MS2*SBA*Yuk1*A0(MS2))/(Mh02*PI2*SW2)

 amplitudes(18) = (0.375D0*EL2*MB2*SBA*Yuk1*A0(MB2))/(Mh02*PI2*SW2)

 amplitudes(19) = (0.125D0*CBA*EL2*ME2*Yuk5*A0(ME2))/(MHH2*PI2*SW2)

 amplitudes(20) = (0.125D0*CBA*EL2*MM2*Yuk5*A0(MM2))/(MHH2*PI2*SW2)

 amplitudes(21) = (0.125D0*CBA*EL2*ML2*Yuk5*A0(ML2))/(MHH2*PI2*SW2)

 amplitudes(22) = (0.375D0*CBA*EL2*MU2*SA*A0(MU2))/(MHH2*PI2*SB*SW2)

 amplitudes(23) = (0.375D0*CBA*EL2*MC2*SA*A0(MC2))/(MHH2*PI2*SB*SW2)

 amplitudes(24) = (0.375D0*CBA*EL2*MT2*SA*A0(MT2))/(MHH2*PI2*SB*SW2)

 amplitudes(25) = (0.375D0*CBA*EL2*MD2*Yuk2*A0(MD2))/(MHH2*PI2*SW2)

 amplitudes(26) = (0.375D0*CBA*EL2*MS2*Yuk2*A0(MS2))/(MHH2*PI2*SW2)

 amplitudes(27) = (0.375D0*CBA*EL2*MB2*Yuk2*A0(MB2))/(MHH2*PI2*SW2)

 amplitudes(28) = (-0.046875D0*SBA*(2.D0*CAB*EL2*Mh02 + EL2*Mh02*S2A*SBA - 4.D0*CAB*CBA2*Lambda5*MW2*SW2)*A0(Mh02))/(Mh02*PI2*S2B&
  &*SW2)

 amplitudes(29) = (0.015625D0*SBA2*(EL2*(Mh02 + 2.D0*MHH2)*S2A - 2.D0*Lambda5*MW2*(3.D0*S2A + S2B)*SW2)*A0(MHH2))/(Mh02*PI2*S2B*S&
  &W2)

 amplitudes(30) = (0.015625D0*SBA*(EL2*(-2.D0*MA02 + Mh02)*S2B*SBA + CAB*(-2.D0*EL2*Mh02 + 4.D0*Lambda5*MW2*SW2))*A0(MA02))/(Mh02&
  &*PI2*S2B*SW2)

 amplitudes(31) = (-0.015625D0*EL2*SBA2*A0(MZ2))/(PI2*SW2)

 amplitudes(32) = (0.03125D0*SBA*(EL2*(Mh02 - 2.D0*MHp2)*S2B*SBA + CAB*(-2.D0*EL2*Mh02 + 4.D0*Lambda5*MW2*SW2))*A0(MHp2))/(Mh02*P&
  &I2*S2B*SW2)

 amplitudes(33) = (-0.03125D0*EL2*SBA2*A0(MW2))/(PI2*SW2)

 amplitudes(34) = (-0.015625D0*CBA2*(EL2*(2.D0*Mh02 + MHH2)*S2A + 2.D0*Lambda5*MW2*(-3.D0*S2A + S2B)*SW2)*A0(Mh02))/(MHH2*PI2*S2B&
  &*SW2)

 amplitudes(35) = (0.046875D0*CBA*(CBA*EL2*MHH2*S2A - 2.D0*EL2*MHH2*SAB + 4.D0*Lambda5*MW2*SAB*SBA2*SW2)*A0(MHH2))/(MHH2*PI2*S2B*&
  &SW2)

 amplitudes(36) = (0.015625D0*CBA*(CBA*EL2*(-2.D0*MA02 + MHH2)*S2B - 2.D0*EL2*MHH2*SAB + 4.D0*Lambda5*MW2*SAB*SW2)*A0(MA02))/(MHH&
  &2*PI2*S2B*SW2)

 amplitudes(37) = (-0.015625D0*CBA2*EL2*A0(MZ2))/(PI2*SW2)

 amplitudes(38) = (0.03125D0*CBA*(CBA*EL2*(MHH2 - 2.D0*MHp2)*S2B - 2.D0*EL2*MHH2*SAB + 4.D0*Lambda5*MW2*SAB*SW2)*A0(MHp2))/(MHH2*&
  &PI2*S2B*SW2)

 amplitudes(39) = (-0.03125D0*CBA2*EL2*A0(MW2))/(PI2*SW2)

 amplitudes(40) = (0.03125D0*EL2*MW2*SBA2*A0(MZ2))/(CW2*Mh02*PI2*SW2)

 amplitudes(41) = (0.03125D0*EL2*MW2*SBA2*A0(MW2))/(Mh02*PI2*SW2)

 amplitudes(42) = (0.03125D0*EL2*MW2*SBA2*A0(MW2))/(Mh02*PI2*SW2)

 amplitudes(43) = (0.03125D0*CBA2*EL2*MW2*A0(MZ2))/(CW2*MHH2*PI2*SW2)

 amplitudes(44) = (0.03125D0*CBA2*EL2*MW2*A0(MW2))/(MHH2*PI2*SW2)

 amplitudes(45) = (0.03125D0*CBA2*EL2*MW2*A0(MW2))/(MHH2*PI2*SW2)

 amplitudes(46) = (0.0625D0*EL2*MW2*SBA2*(MZ2 - 2.D0*A0(MZ2)))/(CW2*Mh02*PI2*SW2)

 amplitudes(47) = (0.125D0*EL2*MW2*SBA2*(MW2 - 2.D0*A0(MW2)))/(Mh02*PI2*SW2)

 amplitudes(48) = (0.0625D0*CBA2*EL2*MW2*(MZ2 - 2.D0*A0(MZ2)))/(CW2*MHH2*PI2*SW2)

 amplitudes(49) = (0.125D0*CBA2*EL2*MW2*(MW2 - 2.D0*A0(MW2)))/(MHH2*PI2*SW2)

 amplitudes(50) = 0.D0

 amplitudes(51) = 0.D0

 amplitudes(52) = 0.D0

 amplitudes(53) = 0.D0

 amplitudes(54) = 0.D0

 amplitudes(55) = 0.D0

 amplitudes(56) = 0.D0

 amplitudes(57) = 0.D0

 amplitudes(58) = 0.D0

 amplitudes(59) = 0.D0

 amplitudes(60) = 0.D0

 amplitudes(61) = 0.D0

 amplitudes(62) = 0.D0

 amplitudes(63) = 0.D0

 amplitudes(64) = 0.D0

 amplitudes(65) = 0.D0

 amplitudes(66) = 0.D0

 amplitudes(67) = 0.D0

 amplitudes(68) = 0.D0

 amplitudes(69) = 0.D0

 amplitudes(70) = 0.D0

 amplitudes(71) = 0.D0

 amplitudes(72) = 0.D0

 amplitudes(73) = 0.D0

 amplitudes(74) = 0.D0

 amplitudes(75) = 0.D0

 amplitudes(76) = 0.D0

 amplitudes(77) = 0.D0

 amplitudes(78) = 0.D0

 amplitudes(79) = 0.D0

 amplitudes(80) = 0.D0

 amplitudes(81) = (0.003472222222222222D0*EL2*(2.D0*(3.D0*ME2 - 1.D0*x)*x + 3.D0*(ME2 - 2.D0*x)*A0(ME2) - 3.D0*B0(x, 0.D0, ME2)*(&
  &ME2*x + DBLE(ME**INT(4.D0)) - 2.D0*DBLE(x**INT(2.D0)))))/(PI2*SW2*x)

 amplitudes(82) = (0.003472222222222222D0*EL2*(2.D0*(3.D0*MM2 - 1.D0*x)*x + 3.D0*(MM2 - 2.D0*x)*A0(MM2) - 3.D0*B0(x, 0.D0, MM2)*(&
  &MM2*x + DBLE(MM**INT(4.D0)) - 2.D0*DBLE(x**INT(2.D0)))))/(PI2*SW2*x)

 amplitudes(83) = (0.003472222222222222D0*EL2*(2.D0*(3.D0*ML2 - 1.D0*x)*x + 3.D0*(ML2 - 2.D0*x)*A0(ML2) - 3.D0*B0(x, 0.D0, ML2)*(&
  &ML2*x + DBLE(ML**INT(4.D0)) - 2.D0*DBLE(x**INT(2.D0)))))/(PI2*SW2*x)

 amplitudes(84) = (-0.010416666666666666D0*CKM11*CKMC11*EL2*(-6.D0*MD2*x - 6.D0*MU2*x + (-3.D0*MD2 + 3.D0*MU2 + 6.D0*x)*A0(MD2) +&
  & 3.D0*(MD2 - 1.D0*MU2 + 2.D0*x)*A0(MU2) - 6.D0*MD2*MU2*B0(x, MD2, MU2) + 3.D0*MD2*x*B0(x, MD2, MU2) + 3.D0*MU2*x*B0(x, MD2, MU&
  &2) + 3.D0*B0(x, MD2, MU2)*DBLE(MD**INT(4.D0)) + 3.D0*B0(x, MD2, MU2)*DBLE(MU**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) - 6.D0*B0(x&
  &, MD2, MU2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(85) = (-0.010416666666666666D0*CKM21*CKMC21*EL2*(-6.D0*MC2*x - 6.D0*MD2*x + (-3.D0*MC2 + 3.D0*MD2 + 6.D0*x)*A0(MC2) +&
  & 3.D0*(MC2 - 1.D0*MD2 + 2.D0*x)*A0(MD2) - 6.D0*MC2*MD2*B0(x, MC2, MD2) + 3.D0*MC2*x*B0(x, MC2, MD2) + 3.D0*MD2*x*B0(x, MC2, MD&
  &2) + 3.D0*B0(x, MC2, MD2)*DBLE(MC**INT(4.D0)) + 3.D0*B0(x, MC2, MD2)*DBLE(MD**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) - 6.D0*B0(x&
  &, MC2, MD2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(86) = (-0.010416666666666666D0*CKM31*CKMC31*EL2*(-6.D0*MD2*x - 6.D0*MT2*x + (-3.D0*MD2 + 3.D0*MT2 + 6.D0*x)*A0(MD2) +&
  & 3.D0*(MD2 - 1.D0*MT2 + 2.D0*x)*A0(MT2) - 6.D0*MD2*MT2*B0(x, MD2, MT2) + 3.D0*MD2*x*B0(x, MD2, MT2) + 3.D0*MT2*x*B0(x, MD2, MT&
  &2) + 3.D0*B0(x, MD2, MT2)*DBLE(MD**INT(4.D0)) + 3.D0*B0(x, MD2, MT2)*DBLE(MT**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) - 6.D0*B0(x&
  &, MD2, MT2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(87) = (-0.010416666666666666D0*CKM12*CKMC12*EL2*(-6.D0*MS2*x - 6.D0*MU2*x + (-3.D0*MS2 + 3.D0*MU2 + 6.D0*x)*A0(MS2) +&
  & 3.D0*(MS2 - 1.D0*MU2 + 2.D0*x)*A0(MU2) - 6.D0*MS2*MU2*B0(x, MS2, MU2) + 3.D0*MS2*x*B0(x, MS2, MU2) + 3.D0*MU2*x*B0(x, MS2, MU&
  &2) + 3.D0*B0(x, MS2, MU2)*DBLE(MS**INT(4.D0)) + 3.D0*B0(x, MS2, MU2)*DBLE(MU**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) - 6.D0*B0(x&
  &, MS2, MU2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(88) = (-0.010416666666666666D0*CKM22*CKMC22*EL2*(-6.D0*MC2*x - 6.D0*MS2*x + (-3.D0*MC2 + 3.D0*MS2 + 6.D0*x)*A0(MC2) +&
  & 3.D0*(MC2 - 1.D0*MS2 + 2.D0*x)*A0(MS2) - 6.D0*MC2*MS2*B0(x, MC2, MS2) + 3.D0*MC2*x*B0(x, MC2, MS2) + 3.D0*MS2*x*B0(x, MC2, MS&
  &2) + 3.D0*B0(x, MC2, MS2)*DBLE(MC**INT(4.D0)) + 3.D0*B0(x, MC2, MS2)*DBLE(MS**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) - 6.D0*B0(x&
  &, MC2, MS2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(89) = (-0.010416666666666666D0*CKM32*CKMC32*EL2*(-6.D0*MS2*x - 6.D0*MT2*x + (-3.D0*MS2 + 3.D0*MT2 + 6.D0*x)*A0(MS2) +&
  & 3.D0*(MS2 - 1.D0*MT2 + 2.D0*x)*A0(MT2) - 6.D0*MS2*MT2*B0(x, MS2, MT2) + 3.D0*MS2*x*B0(x, MS2, MT2) + 3.D0*MT2*x*B0(x, MS2, MT&
  &2) + 3.D0*B0(x, MS2, MT2)*DBLE(MS**INT(4.D0)) + 3.D0*B0(x, MS2, MT2)*DBLE(MT**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) - 6.D0*B0(x&
  &, MS2, MT2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(90) = (-0.010416666666666666D0*CKM13*CKMC13*EL2*(-6.D0*MB2*x - 6.D0*MU2*x + (-3.D0*MB2 + 3.D0*MU2 + 6.D0*x)*A0(MB2) +&
  & 3.D0*(MB2 - 1.D0*MU2 + 2.D0*x)*A0(MU2) - 6.D0*MB2*MU2*B0(x, MB2, MU2) + 3.D0*MB2*x*B0(x, MB2, MU2) + 3.D0*MU2*x*B0(x, MB2, MU&
  &2) + 3.D0*B0(x, MB2, MU2)*DBLE(MB**INT(4.D0)) + 3.D0*B0(x, MB2, MU2)*DBLE(MU**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) - 6.D0*B0(x&
  &, MB2, MU2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(91) = (-0.010416666666666666D0*CKM23*CKMC23*EL2*(-6.D0*MB2*x - 6.D0*MC2*x + (-3.D0*MB2 + 3.D0*MC2 + 6.D0*x)*A0(MB2) +&
  & 3.D0*(MB2 - 1.D0*MC2 + 2.D0*x)*A0(MC2) - 6.D0*MB2*MC2*B0(x, MB2, MC2) + 3.D0*MB2*x*B0(x, MB2, MC2) + 3.D0*MC2*x*B0(x, MB2, MC&
  &2) + 3.D0*B0(x, MB2, MC2)*DBLE(MB**INT(4.D0)) + 3.D0*B0(x, MB2, MC2)*DBLE(MC**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) - 6.D0*B0(x&
  &, MB2, MC2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(92) = (-0.010416666666666666D0*CKM33*CKMC33*EL2*(-6.D0*MB2*x - 6.D0*MT2*x + (-3.D0*MB2 + 3.D0*MT2 + 6.D0*x)*A0(MB2) +&
  & 3.D0*(MB2 - 1.D0*MT2 + 2.D0*x)*A0(MT2) - 6.D0*MB2*MT2*B0(x, MB2, MT2) + 3.D0*MB2*x*B0(x, MB2, MT2) + 3.D0*MT2*x*B0(x, MB2, MT&
  &2) + 3.D0*B0(x, MB2, MT2)*DBLE(MB**INT(4.D0)) + 3.D0*B0(x, MB2, MT2)*DBLE(MT**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) - 6.D0*B0(x&
  &, MB2, MT2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(93) = (0.001736111111111111D0*CBA2*EL2*(-6.D0*Mh02*x - 6.D0*MHp2*x - 3.D0*(Mh02 - 1.D0*MHp2 + x)*A0(Mh02) + 3.D0*(Mh0&
  &2 - 1.D0*MHp2 - 1.D0*x)*A0(MHp2) - 6.D0*Mh02*MHp2*B0(x, Mh02, MHp2) - 6.D0*Mh02*x*B0(x, Mh02, MHp2) - 6.D0*MHp2*x*B0(x, Mh02, &
  &MHp2) + 3.D0*B0(x, Mh02, MHp2)*DBLE(Mh0**INT(4.D0)) + 3.D0*B0(x, Mh02, MHp2)*DBLE(MHp**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) + &
  &3.D0*B0(x, Mh02, MHp2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(94) = (0.001736111111111111D0*EL2*SBA2*(-6.D0*MHH2*x - 6.D0*MHp2*x - 3.D0*(MHH2 - 1.D0*MHp2 + x)*A0(MHH2) + 3.D0*(MHH&
  &2 - 1.D0*MHp2 - 1.D0*x)*A0(MHp2) - 6.D0*MHH2*MHp2*B0(x, MHH2, MHp2) - 6.D0*MHH2*x*B0(x, MHH2, MHp2) - 6.D0*MHp2*x*B0(x, MHH2, &
  &MHp2) + 3.D0*B0(x, MHH2, MHp2)*DBLE(MHH**INT(4.D0)) + 3.D0*B0(x, MHH2, MHp2)*DBLE(MHp**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) + &
  &3.D0*B0(x, MHH2, MHp2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(95) = (0.001736111111111111D0*EL2*(-6.D0*MA02*x - 6.D0*MHp2*x - 3.D0*(MA02 - 1.D0*MHp2 + x)*A0(MA02) + 3.D0*(MA02 - 1&
  &.D0*MHp2 - 1.D0*x)*A0(MHp2) - 6.D0*MA02*MHp2*B0(x, MA02, MHp2) - 6.D0*MA02*x*B0(x, MA02, MHp2) - 6.D0*MHp2*x*B0(x, MA02, MHp2)&
  & + 3.D0*B0(x, MA02, MHp2)*DBLE(MA0**INT(4.D0)) + 3.D0*B0(x, MA02, MHp2)*DBLE(MHp**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) + 3.D0*&
  &B0(x, MA02, MHp2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(96) = (0.001736111111111111D0*EL2*SBA2*(-6.D0*Mh02*x - 6.D0*MW2*x - 3.D0*(Mh02 - 1.D0*MW2 + x)*A0(Mh02) + 3.D0*(Mh02 &
  &- 1.D0*MW2 - 1.D0*x)*A0(MW2) - 6.D0*Mh02*MW2*B0(x, Mh02, MW2) - 6.D0*Mh02*x*B0(x, Mh02, MW2) - 6.D0*MW2*x*B0(x, Mh02, MW2) + 3&
  &.D0*B0(x, Mh02, MW2)*DBLE(Mh0**INT(4.D0)) + 3.D0*B0(x, Mh02, MW2)*DBLE(MW**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) + 3.D0*B0(x, M&
  &h02, MW2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(97) = (0.001736111111111111D0*CBA2*EL2*(-6.D0*MHH2*x - 6.D0*MW2*x - 3.D0*(MHH2 - 1.D0*MW2 + x)*A0(MHH2) + 3.D0*(MHH2 &
  &- 1.D0*MW2 - 1.D0*x)*A0(MW2) - 6.D0*MHH2*MW2*B0(x, MHH2, MW2) - 6.D0*MHH2*x*B0(x, MHH2, MW2) - 6.D0*MW2*x*B0(x, MHH2, MW2) + 3&
  &.D0*B0(x, MHH2, MW2)*DBLE(MHH**INT(4.D0)) + 3.D0*B0(x, MHH2, MW2)*DBLE(MW**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) + 3.D0*B0(x, M&
  &HH2, MW2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(98) = (0.001736111111111111D0*EL2*(-6.D0*MW2*x - 6.D0*MZ2*x - 3.D0*(MW2 - 1.D0*MZ2 + x)*A0(MW2) + 3.D0*(MW2 - 1.D0*MZ&
  &2 - 1.D0*x)*A0(MZ2) - 6.D0*MW2*MZ2*B0(x, MW2, MZ2) - 6.D0*MW2*x*B0(x, MW2, MZ2) - 6.D0*MZ2*x*B0(x, MW2, MZ2) + 3.D0*B0(x, MW2,&
  & MZ2)*DBLE(MW**INT(4.D0)) + 3.D0*B0(x, MW2, MZ2)*DBLE(MZ**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) + 3.D0*B0(x, MW2, MZ2)*DBLE(x**&
  &INT(2.D0))))/(PI2*SW2*x)

 amplitudes(99) = (-0.001736111111111111D0*EL2*(2.D0*x*(-3.D0*MW2 + x) - 3.D0*(MW2 + x)*A0(MW2) + 3.D0*B0(x, 0.D0, MW2)*DBLE((MW2&
  & - 1.D0*x)**INT(2.D0))))/(PI2*x)

 amplitudes(100) = (0.001736111111111111D0*CW2*EL2*(6.D0*MW2*x + 6.D0*MZ2*x + 3.D0*(MW2 - 1.D0*MZ2 + x)*A0(MW2) + 3.D0*(-1.D0*MW2&
  & + MZ2 + x)*A0(MZ2) + 6.D0*MW2*MZ2*B0(x, MW2, MZ2) + 6.D0*MW2*x*B0(x, MW2, MZ2) + 6.D0*MZ2*x*B0(x, MW2, MZ2) - 3.D0*B0(x, MW2,&
  & MZ2)*DBLE(MW**INT(4.D0)) - 3.D0*B0(x, MW2, MZ2)*DBLE(MZ**INT(4.D0)) - 2.D0*DBLE(x**INT(2.D0)) - 3.D0*B0(x, MW2, MZ2)*DBLE(x**&
  &INT(2.D0))))/(PI2*SW2*x)

 amplitudes(101) = (-0.001736111111111111D0*EL2*(2.D0*x*(-3.D0*MW2 + x) - 3.D0*(MW2 + x)*A0(MW2) + 3.D0*B0(x, 0.D0, MW2)*DBLE((MW&
  &2 - 1.D0*x)**INT(2.D0))))/(PI2*x)

 amplitudes(102) = (0.001736111111111111D0*CW2*EL2*(6.D0*MW2*x + 6.D0*MZ2*x + 3.D0*(MW2 - 1.D0*MZ2 + x)*A0(MW2) + 3.D0*(-1.D0*MW2&
  & + MZ2 + x)*A0(MZ2) + 6.D0*MW2*MZ2*B0(x, MW2, MZ2) + 6.D0*MW2*x*B0(x, MW2, MZ2) + 6.D0*MZ2*x*B0(x, MW2, MZ2) - 3.D0*B0(x, MW2,&
  & MZ2)*DBLE(MW**INT(4.D0)) - 3.D0*B0(x, MW2, MZ2)*DBLE(MZ**INT(4.D0)) - 2.D0*DBLE(x**INT(2.D0)) - 3.D0*B0(x, MW2, MZ2)*DBLE(x**&
  &INT(2.D0))))/(PI2*SW2*x)

 amplitudes(103) = (-0.003472222222222222D0*EL2*(2.D0*x*(-3.D0*MW2 + x) + 3.D0*(5.D0*MW2 + 11.D0*x)*A0(MW2) + B0(x, 0.D0, MW2)*(4&
  &8.D0*MW2*x - 15.D0*DBLE(MW**INT(4.D0)) + 57.D0*DBLE(x**INT(2.D0)))))/(PI2*x)

 amplitudes(104) = (-0.003472222222222222D0*CW2*EL2*(-6.D0*MW2*x - 6.D0*MZ2*x + 3.D0*(5.D0*MW2 - 5.D0*MZ2 + 11.D0*x)*A0(MW2) + (-&
  &15.D0*MW2 + 15.D0*MZ2 + 33.D0*x)*A0(MZ2) + 30.D0*MW2*MZ2*B0(x, MW2, MZ2) + 48.D0*MW2*x*B0(x, MW2, MZ2) + 48.D0*MZ2*x*B0(x, MW2&
  &, MZ2) - 15.D0*B0(x, MW2, MZ2)*DBLE(MW**INT(4.D0)) - 15.D0*B0(x, MW2, MZ2)*DBLE(MZ**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) + 57.&
  &D0*B0(x, MW2, MZ2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(105) = (0.0625D0*EL2*MW2*B0(x, 0.D0, MW2))/PI2

 amplitudes(106) = (0.0625D0*EL2*MW2*SW2*B0(x, MW2, MZ2))/(CW2*PI2)

 amplitudes(107) = (0.0625D0*EL2*MW2*SBA2*B0(x, Mh02, MW2))/(PI2*SW2)

 amplitudes(108) = (0.0625D0*CBA2*EL2*MW2*B0(x, MHH2, MW2))/(PI2*SW2)

  totalAmplitude = (0D0,0D0)
 do j=1,108
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 SelfWpWpAlter = totalAmplitude
end function SelfWpWpAlter

