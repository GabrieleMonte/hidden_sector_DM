double complex function SelfCCLeftAlter(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(106)

 amplitudes(1) = 0.D0

 amplitudes(2) = 0.D0

 amplitudes(3) = 0.D0

 amplitudes(4) = 0.D0

 amplitudes(5) = 0.D0

 amplitudes(6) = 0.D0

 amplitudes(7) = 0.D0

 amplitudes(8) = 0.D0

 amplitudes(9) = 0.D0

 amplitudes(10) = 0.D0

 amplitudes(11) = 0.D0

 amplitudes(12) = 0.D0

 amplitudes(13) = 0.D0

 amplitudes(14) = 0.D0

 amplitudes(15) = 0.D0

 amplitudes(16) = 0.D0

 amplitudes(17) = 0.D0

 amplitudes(18) = 0.D0

 amplitudes(19) = 0D0

 amplitudes(20) = 0D0

 amplitudes(21) = 0D0

 amplitudes(22) = 0D0

 amplitudes(23) = 0D0

 amplitudes(24) = 0D0

 amplitudes(25) = 0D0

 amplitudes(26) = 0D0

 amplitudes(27) = 0D0

 amplitudes(28) = 0D0

 amplitudes(29) = 0D0

 amplitudes(30) = 0D0

 amplitudes(31) = 0D0

 amplitudes(32) = 0D0

 amplitudes(33) = 0D0

 amplitudes(34) = 0D0

 amplitudes(35) = 0D0

 amplitudes(36) = 0D0

 amplitudes(37) = 0.D0

 amplitudes(38) = 0.D0

 amplitudes(39) = 0.D0

 amplitudes(40) = 0.D0

 amplitudes(41) = 0.D0

 amplitudes(42) = 0.D0

 amplitudes(43) = 0.D0

 amplitudes(44) = 0.D0

 amplitudes(45) = 0.D0

 amplitudes(46) = 0.D0

 amplitudes(47) = 0.D0

 amplitudes(48) = 0.D0

 amplitudes(49) = 0.D0

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

 amplitudes(81) = 0.D0

 amplitudes(82) = 0.D0

 amplitudes(83) = 0.D0

 amplitudes(84) = 0.D0

 amplitudes(85) = 0.D0

 amplitudes(86) = 0.D0

 amplitudes(87) = 0.D0

 amplitudes(88) = 0.D0

 amplitudes(89) = 0.D0

 amplitudes(90) = 0.D0

 amplitudes(91) = 0.D0

 amplitudes(92) = (0.0078125D0*CA2*EL2*MC2*(-1.D0*A0(MC2) + A0(Mh02) + MC2*B0(x, MC2, Mh02) - 1.D0*Mh02*B0(x, MC2, Mh02) + x*B0(x&
  &, MC2, Mh02)))/ (MW2*PI2*SB2*SW2*x)

 amplitudes(93) = (0.0078125D0*EL2*MC2*SA2*(-1.D0*A0(MC2) + A0(MHH2) + MC2*B0(x, MC2, MHH2) - 1.D0*MHH2*B0(x, MC2, MHH2) + x*B0(x&
  &, MC2, MHH2)))/ (MW2*PI2*SB2*SW2*x)

 amplitudes(94) = (0.0078125D0*EL2*MC2*(A0(MA02) - 1.D0*A0(MC2) - 1.D0*MA02*B0(x, MA02, MC2) + MC2*B0(x, MA02, MC2) + x*B0(x, MA0&
  &2, MC2)))/(MW2*PI2*SW2*TB2*x)

 amplitudes(95) = (0.0078125D0*EL2*MC2*(-1.D0*A0(MC2) + A0(MZ2) + MC2*B0(x, MC2, MZ2) - 1.D0*MZ2*B0(x, MC2, MZ2) + x*B0(x, MC2, M&
  &Z2)))/(MW2*PI2*SW2*x)

 amplitudes(96) = (0.015625D0*CKM21*CKMC21*EL2*(-1.D0*MD2*TB2*A0(MD2)*DBLE(Yuk3**INT(2.D0)) + MD2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D&
  &0)) - 1.D0*MD2*MHp2*TB2*B0(x, MD2, MHp2)*DBLE(Yuk3**INT(2.D0)) + MD2*TB2*x*B0(x, MD2, MHp2)*DBLE(Yuk3**INT(2.D0)) + TB2*B0(x, &
  &MD2, MHp2)*DBLE(MD**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(97) = (0.015625D0*CKM22*CKMC22*EL2*(MS2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D0)) - 1.D0*MS2*TB2*A0(MS2)*DBLE(Yuk3**INT(2.D0&
  &)) - 1.D0*MHp2*MS2*TB2*B0(x, MHp2, MS2)*DBLE(Yuk3**INT(2.D0)) + MS2*TB2*x*B0(x, MHp2, MS2)*DBLE(Yuk3**INT(2.D0)) + TB2*B0(x, M&
  &Hp2, MS2)*DBLE(MS**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(98) = (0.015625D0*CKM23*CKMC23*EL2*(-1.D0*MB2*TB2*A0(MB2)*DBLE(Yuk3**INT(2.D0)) + MB2*TB2*A0(MHp2)*DBLE(Yuk3**INT(2.D&
  &0)) - 1.D0*MB2*MHp2*TB2*B0(x, MB2, MHp2)*DBLE(Yuk3**INT(2.D0)) + MB2*TB2*x*B0(x, MB2, MHp2)*DBLE(Yuk3**INT(2.D0)) + TB2*B0(x, &
  &MB2, MHp2)*DBLE(MB**INT(4.D0))*DBLE(Yuk3**INT(2.D0))))/(MW2*PI2*SW2*TB2*x)

 amplitudes(99) = (0.015625D0*CKM21*CKMC21*EL2*(-1.D0*MD2*A0(MD2) + MD2*A0(MW2) - 1.D0*MD2*MW2*B0(x, MD2, MW2) + MD2*x*B0(x, MD2,&
  & MW2) + B0(x, MD2, MW2)*DBLE(MD**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(100) = (0.015625D0*CKM22*CKMC22*EL2*(-1.D0*MS2*A0(MS2) + MS2*A0(MW2) - 1.D0*MS2*MW2*B0(x, MS2, MW2) + MS2*x*B0(x, MS2&
  &, MW2) + B0(x, MS2, MW2)*DBLE(MS**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(101) = (0.015625D0*CKM23*CKMC23*EL2*(-1.D0*MB2*A0(MB2) + MB2*A0(MW2) - 1.D0*MB2*MW2*B0(x, MB2, MW2) + MB2*x*B0(x, MB2&
  &, MW2) + B0(x, MB2, MW2)*DBLE(MB**INT(4.D0))))/(MW2*PI2*SW2*x)

 amplitudes(102) = (0.027777777777777776D0*EL2*(-1.D0*x - 1.D0*A0(MC2) + MC2*B0(x, 0.D0, MC2) + x*B0(x, 0.D0, MC2)))/(PI2*x)

 amplitudes(103) = (0.001736111111111111D0*EL2*(-9.D0*x + 24.D0*SW2*x - 9.D0*A0(MC2) + 24.D0*SW2*A0(MC2) + 9.D0*A0(MZ2) - 24.D0*S&
  &W2*A0(MZ2) + 9.D0*MC2*B0(x, MC2, MZ2) - 9.D0*MZ2*B0(x, MC2, MZ2) - 24.D0*MC2*SW2*B0(x, MC2, MZ2) + 24.D0*MZ2*SW2*B0(x, MC2, MZ&
  &2) + 9.D0*x*B0(x, MC2, MZ2) - 24.D0*SW2*x*B0(x, MC2, MZ2) - 16.D0*x*DBLE(SW**INT(4.D0)) - 16.D0*A0(MC2)*DBLE(SW**INT(4.D0)) + &
  &16.D0*A0(MZ2)*DBLE(SW**INT(4.D0)) + 16.D0*MC2*B0(x, MC2, MZ2)*DBLE(SW**INT(4.D0)) - 16.D0*MZ2*B0(x, MC2, MZ2)*DBLE(SW**INT(4.D&
  &0)) + 16.D0*x*B0(x, MC2, MZ2)*DBLE(SW**INT(4.D0))))/ (CW2*PI2*SW2*x)

 amplitudes(104) = (0.03125D0*CKM21*CKMC21*EL2*(-1.D0*x - 1.D0*A0(MD2) + A0(MW2) + MD2*B0(x, MD2, MW2) - 1.D0*MW2*B0(x, MD2, MW2)&
  & + x*B0(x, MD2, MW2)))/(PI2*SW2*x)

 amplitudes(105) = (0.03125D0*CKM22*CKMC22*EL2*(-1.D0*x - 1.D0*A0(MS2) + A0(MW2) + MS2*B0(x, MS2, MW2) - 1.D0*MW2*B0(x, MS2, MW2)&
  & + x*B0(x, MS2, MW2)))/(PI2*SW2*x)

 amplitudes(106) = (0.03125D0*CKM23*CKMC23*EL2*(-1.D0*x - 1.D0*A0(MB2) + A0(MW2) + MB2*B0(x, MB2, MW2) - 1.D0*MW2*B0(x, MB2, MW2)&
  & + x*B0(x, MB2, MW2)))/(PI2*SW2*x)

  totalAmplitude = (0D0,0D0)
 do j=1,106
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 SelfCCLeftAlter = totalAmplitude
end function SelfCCLeftAlter

