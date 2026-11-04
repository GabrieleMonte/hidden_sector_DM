double complex function DSelfWpWp(x)
 use constants
 implicit none
#include "looptools.h"
 double precision, intent(in) :: x
 integer :: j
 double complex :: totalAmplitude
 double complex :: amplitudes(37)

 amplitudes(1) = 0.D0

 amplitudes(2) = 0.D0

 amplitudes(3) = 0.D0

 amplitudes(4) = 0.D0

 amplitudes(5) = 0.D0

 amplitudes(6) = 0.D0

 amplitudes(7) = 0.D0

 amplitudes(8) = 0.D0

 amplitudes(9) = 0.D0

 amplitudes(10) = (-0.003472222222222222D0*EL2*DBLE(x**INT(-2.D0))*(2.D0*(3.D0*ME2 - 1.D0*x)*x + 3.D0*(ME2 - 2.D0*x)*A0(ME2) - 3.&
  &D0*B0(x, 0.D0, ME2)*(ME2*x + DBLE(ME**INT(4.D0)) - 2.D0*DBLE(x**INT(2.D0)))))/(PI2*SW2) + (0.003472222222222222D0*EL2*(2.D0*(3&
  &.D0*ME2 - 1.D0*x) - 2.D0*x - 6.D0*A0(ME2) - 3.D0*(ME2 - 4.D0*x)*B0(x, 0.D0, ME2) - 3.D0*DB0(x, 0.D0, ME2)*(ME2*x + DBLE(ME**IN&
  &T(4.D0)) - 2.D0*DBLE(x**INT(2.D0)))))/(PI2*SW2*x)

 amplitudes(11) = (-0.003472222222222222D0*EL2*DBLE(x**INT(-2.D0))*(2.D0*(3.D0*MM2 - 1.D0*x)*x + 3.D0*(MM2 - 2.D0*x)*A0(MM2) - 3.&
  &D0*B0(x, 0.D0, MM2)*(MM2*x + DBLE(MM**INT(4.D0)) - 2.D0*DBLE(x**INT(2.D0)))))/(PI2*SW2) + (0.003472222222222222D0*EL2*(2.D0*(3&
  &.D0*MM2 - 1.D0*x) - 2.D0*x - 6.D0*A0(MM2) - 3.D0*(MM2 - 4.D0*x)*B0(x, 0.D0, MM2) - 3.D0*DB0(x, 0.D0, MM2)*(MM2*x + DBLE(MM**IN&
  &T(4.D0)) - 2.D0*DBLE(x**INT(2.D0)))))/(PI2*SW2*x)

 amplitudes(12) = (-0.003472222222222222D0*EL2*DBLE(x**INT(-2.D0))*(2.D0*(3.D0*ML2 - 1.D0*x)*x + 3.D0*(ML2 - 2.D0*x)*A0(ML2) - 3.&
  &D0*B0(x, 0.D0, ML2)*(ML2*x + DBLE(ML**INT(4.D0)) - 2.D0*DBLE(x**INT(2.D0)))))/(PI2*SW2) + (0.003472222222222222D0*EL2*(2.D0*(3&
  &.D0*ML2 - 1.D0*x) - 2.D0*x - 6.D0*A0(ML2) - 3.D0*(ML2 - 4.D0*x)*B0(x, 0.D0, ML2) - 3.D0*DB0(x, 0.D0, ML2)*(ML2*x + DBLE(ML**IN&
  &T(4.D0)) - 2.D0*DBLE(x**INT(2.D0)))))/(PI2*SW2*x)

 amplitudes(13) = (0.010416666666666666D0*CKM11*CKMC11*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MD2*x - 6.D0*MU2*x + (-3.D0*MD2 + 3.D0*MU2 &
  &+ 6.D0*x)*A0(MD2) + 3.D0*(MD2 - 1.D0*MU2 + 2.D0*x)*A0(MU2) - 6.D0*MD2*MU2*B0(x, MD2, MU2) + 3.D0*MD2*x*B0(x, MD2, MU2) + 3.D0*&
  &MU2*x*B0(x, MD2, MU2) + 3.D0*B0(x, MD2, MU2)*DBLE(MD**INT(4.D0)) + 3.D0*B0(x, MD2, MU2)*DBLE(MU**INT(4.D0)) + 2.D0*DBLE(x**INT&
  &(2.D0)) - 6.D0*B0(x, MD2, MU2)*DBLE(x**INT(2.D0))))/(PI2*SW2) - (0.010416666666666666D0*CKM11*CKMC11*EL2*(-6.D0*MD2 - 6.D0*MU2&
  & + 4.D0*x + 6.D0*A0(MD2) + 6.D0*A0(MU2) + 3.D0*MD2*B0(x, MD2, MU2) + 3.D0*MU2*B0(x, MD2, MU2) - 12.D0*x*B0(x, MD2, MU2) - 6.D0&
  &*MD2*MU2*DB0(x, MD2, MU2) + 3.D0*MD2*x*DB0(x, MD2, MU2) + 3.D0*MU2*x*DB0(x, MD2, MU2) + 3.D0*DB0(x, MD2, MU2)*DBLE(MD**INT(4.D&
  &0)) + 3.D0*DB0(x, MD2, MU2)*DBLE(MU**INT(4.D0)) - 6.D0*DB0(x, MD2, MU2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(14) = (0.010416666666666666D0*CKM21*CKMC21*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MC2*x - 6.D0*MD2*x + (-3.D0*MC2 + 3.D0*MD2 &
  &+ 6.D0*x)*A0(MC2) + 3.D0*(MC2 - 1.D0*MD2 + 2.D0*x)*A0(MD2) - 6.D0*MC2*MD2*B0(x, MC2, MD2) + 3.D0*MC2*x*B0(x, MC2, MD2) + 3.D0*&
  &MD2*x*B0(x, MC2, MD2) + 3.D0*B0(x, MC2, MD2)*DBLE(MC**INT(4.D0)) + 3.D0*B0(x, MC2, MD2)*DBLE(MD**INT(4.D0)) + 2.D0*DBLE(x**INT&
  &(2.D0)) - 6.D0*B0(x, MC2, MD2)*DBLE(x**INT(2.D0))))/(PI2*SW2) - (0.010416666666666666D0*CKM21*CKMC21*EL2*(-6.D0*MC2 - 6.D0*MD2&
  & + 4.D0*x + 6.D0*A0(MC2) + 6.D0*A0(MD2) + 3.D0*MC2*B0(x, MC2, MD2) + 3.D0*MD2*B0(x, MC2, MD2) - 12.D0*x*B0(x, MC2, MD2) - 6.D0&
  &*MC2*MD2*DB0(x, MC2, MD2) + 3.D0*MC2*x*DB0(x, MC2, MD2) + 3.D0*MD2*x*DB0(x, MC2, MD2) + 3.D0*DB0(x, MC2, MD2)*DBLE(MC**INT(4.D&
  &0)) + 3.D0*DB0(x, MC2, MD2)*DBLE(MD**INT(4.D0)) - 6.D0*DB0(x, MC2, MD2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(15) = (0.010416666666666666D0*CKM31*CKMC31*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MD2*x - 6.D0*MT2*x + (-3.D0*MD2 + 3.D0*MT2 &
  &+ 6.D0*x)*A0(MD2) + 3.D0*(MD2 - 1.D0*MT2 + 2.D0*x)*A0(MT2) - 6.D0*MD2*MT2*B0(x, MD2, MT2) + 3.D0*MD2*x*B0(x, MD2, MT2) + 3.D0*&
  &MT2*x*B0(x, MD2, MT2) + 3.D0*B0(x, MD2, MT2)*DBLE(MD**INT(4.D0)) + 3.D0*B0(x, MD2, MT2)*DBLE(MT**INT(4.D0)) + 2.D0*DBLE(x**INT&
  &(2.D0)) - 6.D0*B0(x, MD2, MT2)*DBLE(x**INT(2.D0))))/(PI2*SW2) - (0.010416666666666666D0*CKM31*CKMC31*EL2*(-6.D0*MD2 - 6.D0*MT2&
  & + 4.D0*x + 6.D0*A0(MD2) + 6.D0*A0(MT2) + 3.D0*MD2*B0(x, MD2, MT2) + 3.D0*MT2*B0(x, MD2, MT2) - 12.D0*x*B0(x, MD2, MT2) - 6.D0&
  &*MD2*MT2*DB0(x, MD2, MT2) + 3.D0*MD2*x*DB0(x, MD2, MT2) + 3.D0*MT2*x*DB0(x, MD2, MT2) + 3.D0*DB0(x, MD2, MT2)*DBLE(MD**INT(4.D&
  &0)) + 3.D0*DB0(x, MD2, MT2)*DBLE(MT**INT(4.D0)) - 6.D0*DB0(x, MD2, MT2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(16) = (0.010416666666666666D0*CKM12*CKMC12*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MS2*x - 6.D0*MU2*x + (-3.D0*MS2 + 3.D0*MU2 &
  &+ 6.D0*x)*A0(MS2) + 3.D0*(MS2 - 1.D0*MU2 + 2.D0*x)*A0(MU2) - 6.D0*MS2*MU2*B0(x, MS2, MU2) + 3.D0*MS2*x*B0(x, MS2, MU2) + 3.D0*&
  &MU2*x*B0(x, MS2, MU2) + 3.D0*B0(x, MS2, MU2)*DBLE(MS**INT(4.D0)) + 3.D0*B0(x, MS2, MU2)*DBLE(MU**INT(4.D0)) + 2.D0*DBLE(x**INT&
  &(2.D0)) - 6.D0*B0(x, MS2, MU2)*DBLE(x**INT(2.D0))))/(PI2*SW2) - (0.010416666666666666D0*CKM12*CKMC12*EL2*(-6.D0*MS2 - 6.D0*MU2&
  & + 4.D0*x + 6.D0*A0(MS2) + 6.D0*A0(MU2) + 3.D0*MS2*B0(x, MS2, MU2) + 3.D0*MU2*B0(x, MS2, MU2) - 12.D0*x*B0(x, MS2, MU2) - 6.D0&
  &*MS2*MU2*DB0(x, MS2, MU2) + 3.D0*MS2*x*DB0(x, MS2, MU2) + 3.D0*MU2*x*DB0(x, MS2, MU2) + 3.D0*DB0(x, MS2, MU2)*DBLE(MS**INT(4.D&
  &0)) + 3.D0*DB0(x, MS2, MU2)*DBLE(MU**INT(4.D0)) - 6.D0*DB0(x, MS2, MU2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(17) = (0.010416666666666666D0*CKM22*CKMC22*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MC2*x - 6.D0*MS2*x + (-3.D0*MC2 + 3.D0*MS2 &
  &+ 6.D0*x)*A0(MC2) + 3.D0*(MC2 - 1.D0*MS2 + 2.D0*x)*A0(MS2) - 6.D0*MC2*MS2*B0(x, MC2, MS2) + 3.D0*MC2*x*B0(x, MC2, MS2) + 3.D0*&
  &MS2*x*B0(x, MC2, MS2) + 3.D0*B0(x, MC2, MS2)*DBLE(MC**INT(4.D0)) + 3.D0*B0(x, MC2, MS2)*DBLE(MS**INT(4.D0)) + 2.D0*DBLE(x**INT&
  &(2.D0)) - 6.D0*B0(x, MC2, MS2)*DBLE(x**INT(2.D0))))/(PI2*SW2) - (0.010416666666666666D0*CKM22*CKMC22*EL2*(-6.D0*MC2 - 6.D0*MS2&
  & + 4.D0*x + 6.D0*A0(MC2) + 6.D0*A0(MS2) + 3.D0*MC2*B0(x, MC2, MS2) + 3.D0*MS2*B0(x, MC2, MS2) - 12.D0*x*B0(x, MC2, MS2) - 6.D0&
  &*MC2*MS2*DB0(x, MC2, MS2) + 3.D0*MC2*x*DB0(x, MC2, MS2) + 3.D0*MS2*x*DB0(x, MC2, MS2) + 3.D0*DB0(x, MC2, MS2)*DBLE(MC**INT(4.D&
  &0)) + 3.D0*DB0(x, MC2, MS2)*DBLE(MS**INT(4.D0)) - 6.D0*DB0(x, MC2, MS2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(18) = (0.010416666666666666D0*CKM32*CKMC32*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MS2*x - 6.D0*MT2*x + (-3.D0*MS2 + 3.D0*MT2 &
  &+ 6.D0*x)*A0(MS2) + 3.D0*(MS2 - 1.D0*MT2 + 2.D0*x)*A0(MT2) - 6.D0*MS2*MT2*B0(x, MS2, MT2) + 3.D0*MS2*x*B0(x, MS2, MT2) + 3.D0*&
  &MT2*x*B0(x, MS2, MT2) + 3.D0*B0(x, MS2, MT2)*DBLE(MS**INT(4.D0)) + 3.D0*B0(x, MS2, MT2)*DBLE(MT**INT(4.D0)) + 2.D0*DBLE(x**INT&
  &(2.D0)) - 6.D0*B0(x, MS2, MT2)*DBLE(x**INT(2.D0))))/(PI2*SW2) - (0.010416666666666666D0*CKM32*CKMC32*EL2*(-6.D0*MS2 - 6.D0*MT2&
  & + 4.D0*x + 6.D0*A0(MS2) + 6.D0*A0(MT2) + 3.D0*MS2*B0(x, MS2, MT2) + 3.D0*MT2*B0(x, MS2, MT2) - 12.D0*x*B0(x, MS2, MT2) - 6.D0&
  &*MS2*MT2*DB0(x, MS2, MT2) + 3.D0*MS2*x*DB0(x, MS2, MT2) + 3.D0*MT2*x*DB0(x, MS2, MT2) + 3.D0*DB0(x, MS2, MT2)*DBLE(MS**INT(4.D&
  &0)) + 3.D0*DB0(x, MS2, MT2)*DBLE(MT**INT(4.D0)) - 6.D0*DB0(x, MS2, MT2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(19) = (0.010416666666666666D0*CKM13*CKMC13*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MB2*x - 6.D0*MU2*x + (-3.D0*MB2 + 3.D0*MU2 &
  &+ 6.D0*x)*A0(MB2) + 3.D0*(MB2 - 1.D0*MU2 + 2.D0*x)*A0(MU2) - 6.D0*MB2*MU2*B0(x, MB2, MU2) + 3.D0*MB2*x*B0(x, MB2, MU2) + 3.D0*&
  &MU2*x*B0(x, MB2, MU2) + 3.D0*B0(x, MB2, MU2)*DBLE(MB**INT(4.D0)) + 3.D0*B0(x, MB2, MU2)*DBLE(MU**INT(4.D0)) + 2.D0*DBLE(x**INT&
  &(2.D0)) - 6.D0*B0(x, MB2, MU2)*DBLE(x**INT(2.D0))))/(PI2*SW2) - (0.010416666666666666D0*CKM13*CKMC13*EL2*(-6.D0*MB2 - 6.D0*MU2&
  & + 4.D0*x + 6.D0*A0(MB2) + 6.D0*A0(MU2) + 3.D0*MB2*B0(x, MB2, MU2) + 3.D0*MU2*B0(x, MB2, MU2) - 12.D0*x*B0(x, MB2, MU2) - 6.D0&
  &*MB2*MU2*DB0(x, MB2, MU2) + 3.D0*MB2*x*DB0(x, MB2, MU2) + 3.D0*MU2*x*DB0(x, MB2, MU2) + 3.D0*DB0(x, MB2, MU2)*DBLE(MB**INT(4.D&
  &0)) + 3.D0*DB0(x, MB2, MU2)*DBLE(MU**INT(4.D0)) - 6.D0*DB0(x, MB2, MU2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(20) = (0.010416666666666666D0*CKM23*CKMC23*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MB2*x - 6.D0*MC2*x + (-3.D0*MB2 + 3.D0*MC2 &
  &+ 6.D0*x)*A0(MB2) + 3.D0*(MB2 - 1.D0*MC2 + 2.D0*x)*A0(MC2) - 6.D0*MB2*MC2*B0(x, MB2, MC2) + 3.D0*MB2*x*B0(x, MB2, MC2) + 3.D0*&
  &MC2*x*B0(x, MB2, MC2) + 3.D0*B0(x, MB2, MC2)*DBLE(MB**INT(4.D0)) + 3.D0*B0(x, MB2, MC2)*DBLE(MC**INT(4.D0)) + 2.D0*DBLE(x**INT&
  &(2.D0)) - 6.D0*B0(x, MB2, MC2)*DBLE(x**INT(2.D0))))/(PI2*SW2) - (0.010416666666666666D0*CKM23*CKMC23*EL2*(-6.D0*MB2 - 6.D0*MC2&
  & + 4.D0*x + 6.D0*A0(MB2) + 6.D0*A0(MC2) + 3.D0*MB2*B0(x, MB2, MC2) + 3.D0*MC2*B0(x, MB2, MC2) - 12.D0*x*B0(x, MB2, MC2) - 6.D0&
  &*MB2*MC2*DB0(x, MB2, MC2) + 3.D0*MB2*x*DB0(x, MB2, MC2) + 3.D0*MC2*x*DB0(x, MB2, MC2) + 3.D0*DB0(x, MB2, MC2)*DBLE(MB**INT(4.D&
  &0)) + 3.D0*DB0(x, MB2, MC2)*DBLE(MC**INT(4.D0)) - 6.D0*DB0(x, MB2, MC2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(21) = (0.010416666666666666D0*CKM33*CKMC33*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MB2*x - 6.D0*MT2*x + (-3.D0*MB2 + 3.D0*MT2 &
  &+ 6.D0*x)*A0(MB2) + 3.D0*(MB2 - 1.D0*MT2 + 2.D0*x)*A0(MT2) - 6.D0*MB2*MT2*B0(x, MB2, MT2) + 3.D0*MB2*x*B0(x, MB2, MT2) + 3.D0*&
  &MT2*x*B0(x, MB2, MT2) + 3.D0*B0(x, MB2, MT2)*DBLE(MB**INT(4.D0)) + 3.D0*B0(x, MB2, MT2)*DBLE(MT**INT(4.D0)) + 2.D0*DBLE(x**INT&
  &(2.D0)) - 6.D0*B0(x, MB2, MT2)*DBLE(x**INT(2.D0))))/(PI2*SW2) - (0.010416666666666666D0*CKM33*CKMC33*EL2*(-6.D0*MB2 - 6.D0*MT2&
  & + 4.D0*x + 6.D0*A0(MB2) + 6.D0*A0(MT2) + 3.D0*MB2*B0(x, MB2, MT2) + 3.D0*MT2*B0(x, MB2, MT2) - 12.D0*x*B0(x, MB2, MT2) - 6.D0&
  &*MB2*MT2*DB0(x, MB2, MT2) + 3.D0*MB2*x*DB0(x, MB2, MT2) + 3.D0*MT2*x*DB0(x, MB2, MT2) + 3.D0*DB0(x, MB2, MT2)*DBLE(MB**INT(4.D&
  &0)) + 3.D0*DB0(x, MB2, MT2)*DBLE(MT**INT(4.D0)) - 6.D0*DB0(x, MB2, MT2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(22) = (-0.001736111111111111D0*CBA2*EL2*DBLE(x**INT(-2.D0))*(-6.D0*Mh02*x - 6.D0*MHp2*x - 3.D0*(Mh02 - 1.D0*MHp2 + x)&
  &*A0(Mh02) + 3.D0*(Mh02 - 1.D0*MHp2 - 1.D0*x)*A0(MHp2) - 6.D0*Mh02*MHp2*B0(x, Mh02, MHp2) - 6.D0*Mh02*x*B0(x, Mh02, MHp2) - 6.D&
  &0*MHp2*x*B0(x, Mh02, MHp2) + 3.D0*B0(x, Mh02, MHp2)*DBLE(Mh0**INT(4.D0)) + 3.D0*B0(x, Mh02, MHp2)*DBLE(MHp**INT(4.D0)) + 2.D0*&
  &DBLE(x**INT(2.D0)) + 3.D0*B0(x, Mh02, MHp2)*DBLE(x**INT(2.D0))))/(PI2*SW2) + (0.001736111111111111D0*CBA2*EL2*(-6.D0*Mh02 - 6.&
  &D0*MHp2 + 4.D0*x - 3.D0*A0(Mh02) - 3.D0*A0(MHp2) - 6.D0*Mh02*B0(x, Mh02, MHp2) - 6.D0*MHp2*B0(x, Mh02, MHp2) + 6.D0*x*B0(x, Mh&
  &02, MHp2) - 6.D0*Mh02*MHp2*DB0(x, Mh02, MHp2) - 6.D0*Mh02*x*DB0(x, Mh02, MHp2) - 6.D0*MHp2*x*DB0(x, Mh02, MHp2) + 3.D0*DB0(x, &
  &Mh02, MHp2)*DBLE(Mh0**INT(4.D0)) + 3.D0*DB0(x, Mh02, MHp2)*DBLE(MHp**INT(4.D0)) + 3.D0*DB0(x, Mh02, MHp2)*DBLE(x**INT(2.D0))))&
  &/(PI2*SW2*x)

 amplitudes(23) = (-0.001736111111111111D0*EL2*SBA2*DBLE(x**INT(-2.D0))*(-6.D0*MHH2*x - 6.D0*MHp2*x - 3.D0*(MHH2 - 1.D0*MHp2 + x)&
  &*A0(MHH2) + 3.D0*(MHH2 - 1.D0*MHp2 - 1.D0*x)*A0(MHp2) - 6.D0*MHH2*MHp2*B0(x, MHH2, MHp2) - 6.D0*MHH2*x*B0(x, MHH2, MHp2) - 6.D&
  &0*MHp2*x*B0(x, MHH2, MHp2) + 3.D0*B0(x, MHH2, MHp2)*DBLE(MHH**INT(4.D0)) + 3.D0*B0(x, MHH2, MHp2)*DBLE(MHp**INT(4.D0)) + 2.D0*&
  &DBLE(x**INT(2.D0)) + 3.D0*B0(x, MHH2, MHp2)*DBLE(x**INT(2.D0))))/(PI2*SW2) + (0.001736111111111111D0*EL2*SBA2*(-6.D0*MHH2 - 6.&
  &D0*MHp2 + 4.D0*x - 3.D0*A0(MHH2) - 3.D0*A0(MHp2) - 6.D0*MHH2*B0(x, MHH2, MHp2) - 6.D0*MHp2*B0(x, MHH2, MHp2) + 6.D0*x*B0(x, MH&
  &H2, MHp2) - 6.D0*MHH2*MHp2*DB0(x, MHH2, MHp2) - 6.D0*MHH2*x*DB0(x, MHH2, MHp2) - 6.D0*MHp2*x*DB0(x, MHH2, MHp2) + 3.D0*DB0(x, &
  &MHH2, MHp2)*DBLE(MHH**INT(4.D0)) + 3.D0*DB0(x, MHH2, MHp2)*DBLE(MHp**INT(4.D0)) + 3.D0*DB0(x, MHH2, MHp2)*DBLE(x**INT(2.D0))))&
  &/(PI2*SW2*x)

 amplitudes(24) = (-0.001736111111111111D0*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MA02*x - 6.D0*MHp2*x - 3.D0*(MA02 - 1.D0*MHp2 + x)*A0(M&
  &A02) + 3.D0*(MA02 - 1.D0*MHp2 - 1.D0*x)*A0(MHp2) - 6.D0*MA02*MHp2*B0(x, MA02, MHp2) - 6.D0*MA02*x*B0(x, MA02, MHp2) - 6.D0*MHp&
  &2*x*B0(x, MA02, MHp2) + 3.D0*B0(x, MA02, MHp2)*DBLE(MA0**INT(4.D0)) + 3.D0*B0(x, MA02, MHp2)*DBLE(MHp**INT(4.D0)) + 2.D0*DBLE(&
  &x**INT(2.D0)) + 3.D0*B0(x, MA02, MHp2)*DBLE(x**INT(2.D0))))/(PI2*SW2) + (0.001736111111111111D0*EL2*(-6.D0*MA02 - 6.D0*MHp2 + &
  &4.D0*x - 3.D0*A0(MA02) - 3.D0*A0(MHp2) - 6.D0*MA02*B0(x, MA02, MHp2) - 6.D0*MHp2*B0(x, MA02, MHp2) + 6.D0*x*B0(x, MA02, MHp2) &
  &- 6.D0*MA02*MHp2*DB0(x, MA02, MHp2) - 6.D0*MA02*x*DB0(x, MA02, MHp2) - 6.D0*MHp2*x*DB0(x, MA02, MHp2) + 3.D0*DB0(x, MA02, MHp2&
  &)*DBLE(MA0**INT(4.D0)) + 3.D0*DB0(x, MA02, MHp2)*DBLE(MHp**INT(4.D0)) + 3.D0*DB0(x, MA02, MHp2)*DBLE(x**INT(2.D0))))/(PI2*SW2*&
  &x)

 amplitudes(25) = (-0.001736111111111111D0*EL2*SBA2*DBLE(x**INT(-2.D0))*(-6.D0*Mh02*x - 6.D0*MW2*x - 3.D0*(Mh02 - 1.D0*MW2 + x)*A&
  &0(Mh02) + 3.D0*(Mh02 - 1.D0*MW2 - 1.D0*x)*A0(MW2) - 6.D0*Mh02*MW2*B0(x, Mh02, MW2) - 6.D0*Mh02*x*B0(x, Mh02, MW2) - 6.D0*MW2*x&
  &*B0(x, Mh02, MW2) + 3.D0*B0(x, Mh02, MW2)*DBLE(Mh0**INT(4.D0)) + 3.D0*B0(x, Mh02, MW2)*DBLE(MW**INT(4.D0)) + 2.D0*DBLE(x**INT(&
  &2.D0)) + 3.D0*B0(x, Mh02, MW2)*DBLE(x**INT(2.D0))))/(PI2*SW2) + (0.001736111111111111D0*EL2*SBA2*(-6.D0*Mh02 - 6.D0*MW2 + 4.D0&
  &*x - 3.D0*A0(Mh02) - 3.D0*A0(MW2) - 6.D0*Mh02*B0(x, Mh02, MW2) - 6.D0*MW2*B0(x, Mh02, MW2) + 6.D0*x*B0(x, Mh02, MW2) - 6.D0*Mh&
  &02*MW2*DB0(x, Mh02, MW2) - 6.D0*Mh02*x*DB0(x, Mh02, MW2) - 6.D0*MW2*x*DB0(x, Mh02, MW2) + 3.D0*DB0(x, Mh02, MW2)*DBLE(Mh0**INT&
  &(4.D0)) + 3.D0*DB0(x, Mh02, MW2)*DBLE(MW**INT(4.D0)) + 3.D0*DB0(x, Mh02, MW2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(26) = (-0.001736111111111111D0*CBA2*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MHH2*x - 6.D0*MW2*x - 3.D0*(MHH2 - 1.D0*MW2 + x)*A&
  &0(MHH2) + 3.D0*(MHH2 - 1.D0*MW2 - 1.D0*x)*A0(MW2) - 6.D0*MHH2*MW2*B0(x, MHH2, MW2) - 6.D0*MHH2*x*B0(x, MHH2, MW2) - 6.D0*MW2*x&
  &*B0(x, MHH2, MW2) + 3.D0*B0(x, MHH2, MW2)*DBLE(MHH**INT(4.D0)) + 3.D0*B0(x, MHH2, MW2)*DBLE(MW**INT(4.D0)) + 2.D0*DBLE(x**INT(&
  &2.D0)) + 3.D0*B0(x, MHH2, MW2)*DBLE(x**INT(2.D0))))/(PI2*SW2) + (0.001736111111111111D0*CBA2*EL2*(-6.D0*MHH2 - 6.D0*MW2 + 4.D0&
  &*x - 3.D0*A0(MHH2) - 3.D0*A0(MW2) - 6.D0*MHH2*B0(x, MHH2, MW2) - 6.D0*MW2*B0(x, MHH2, MW2) + 6.D0*x*B0(x, MHH2, MW2) - 6.D0*MH&
  &H2*MW2*DB0(x, MHH2, MW2) - 6.D0*MHH2*x*DB0(x, MHH2, MW2) - 6.D0*MW2*x*DB0(x, MHH2, MW2) + 3.D0*DB0(x, MHH2, MW2)*DBLE(MHH**INT&
  &(4.D0)) + 3.D0*DB0(x, MHH2, MW2)*DBLE(MW**INT(4.D0)) + 3.D0*DB0(x, MHH2, MW2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(27) = (-0.001736111111111111D0*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MW2*x - 6.D0*MZ2*x - 3.D0*(MW2 - 1.D0*MZ2 + x)*A0(MW2) &
  &+ 3.D0*(MW2 - 1.D0*MZ2 - 1.D0*x)*A0(MZ2) - 6.D0*MW2*MZ2*B0(x, MW2, MZ2) - 6.D0*MW2*x*B0(x, MW2, MZ2) - 6.D0*MZ2*x*B0(x, MW2, M&
  &Z2) + 3.D0*B0(x, MW2, MZ2)*DBLE(MW**INT(4.D0)) + 3.D0*B0(x, MW2, MZ2)*DBLE(MZ**INT(4.D0)) + 2.D0*DBLE(x**INT(2.D0)) + 3.D0*B0(&
  &x, MW2, MZ2)*DBLE(x**INT(2.D0))))/(PI2*SW2) + (0.001736111111111111D0*EL2*(-6.D0*MW2 - 6.D0*MZ2 + 4.D0*x - 3.D0*A0(MW2) - 3.D0&
  &*A0(MZ2) - 6.D0*MW2*B0(x, MW2, MZ2) - 6.D0*MZ2*B0(x, MW2, MZ2) + 6.D0*x*B0(x, MW2, MZ2) - 6.D0*MW2*MZ2*DB0(x, MW2, MZ2) - 6.D0&
  &*MW2*x*DB0(x, MW2, MZ2) - 6.D0*MZ2*x*DB0(x, MW2, MZ2) + 3.D0*DB0(x, MW2, MZ2)*DBLE(MW**INT(4.D0)) + 3.D0*DB0(x, MW2, MZ2)*DBLE&
  &(MZ**INT(4.D0)) + 3.D0*DB0(x, MW2, MZ2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(28) = (-0.001736111111111111D0*EL2*(2.D0*x + 2.D0*(-3.D0*MW2 + x) - 3.D0*A0(MW2) - 6.D0*(MW2 - 1.D0*x)*B0(x, 0.D0, MW&
  &2) + 3.D0*DB0(x, 0.D0, MW2)*DBLE((MW2 - 1.D0*x)**INT(2.D0))))/(PI2*x) + (0.001736111111111111D0*EL2*(2.D0*x*(-3.D0*MW2 + x) - &
  &3.D0*(MW2 + x)*A0(MW2) + 3.D0*B0(x, 0.D0, MW2)*DBLE((MW2 - 1.D0*x)**INT(2.D0)))* DBLE(x**INT(-2.D0)))/PI2

 amplitudes(29) = (-0.001736111111111111D0*CW2*EL2*DBLE(x**INT(-2.D0))*(6.D0*MW2*x + 6.D0*MZ2*x + 3.D0*(MW2 - 1.D0*MZ2 + x)*A0(MW&
  &2) + 3.D0*(-1.D0*MW2 + MZ2 + x)*A0(MZ2) + 6.D0*MW2*MZ2*B0(x, MW2, MZ2) + 6.D0*MW2*x*B0(x, MW2, MZ2) + 6.D0*MZ2*x*B0(x, MW2, MZ&
  &2) - 3.D0*B0(x, MW2, MZ2)*DBLE(MW**INT(4.D0)) - 3.D0*B0(x, MW2, MZ2)*DBLE(MZ**INT(4.D0)) - 2.D0*DBLE(x**INT(2.D0)) - 3.D0*B0(x&
  &, MW2, MZ2)*DBLE(x**INT(2.D0))))/(PI2*SW2) + (0.001736111111111111D0*CW2*EL2*(6.D0*MW2 + 6.D0*MZ2 - 4.D0*x + 3.D0*A0(MW2) + 3.&
  &D0*A0(MZ2) + 6.D0*MW2*B0(x, MW2, MZ2) + 6.D0*MZ2*B0(x, MW2, MZ2) - 6.D0*x*B0(x, MW2, MZ2) + 6.D0*MW2*MZ2*DB0(x, MW2, MZ2) + 6.&
  &D0*MW2*x*DB0(x, MW2, MZ2) + 6.D0*MZ2*x*DB0(x, MW2, MZ2) - 3.D0*DB0(x, MW2, MZ2)*DBLE(MW**INT(4.D0)) - 3.D0*DB0(x, MW2, MZ2)*DB&
  &LE(MZ**INT(4.D0)) - 3.D0*DB0(x, MW2, MZ2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(30) = (-0.001736111111111111D0*EL2*(2.D0*x + 2.D0*(-3.D0*MW2 + x) - 3.D0*A0(MW2) - 6.D0*(MW2 - 1.D0*x)*B0(x, 0.D0, MW&
  &2) + 3.D0*DB0(x, 0.D0, MW2)*DBLE((MW2 - 1.D0*x)**INT(2.D0))))/(PI2*x) + (0.001736111111111111D0*EL2*(2.D0*x*(-3.D0*MW2 + x) - &
  &3.D0*(MW2 + x)*A0(MW2) + 3.D0*B0(x, 0.D0, MW2)*DBLE((MW2 - 1.D0*x)**INT(2.D0)))* DBLE(x**INT(-2.D0)))/PI2

 amplitudes(31) = (-0.001736111111111111D0*CW2*EL2*DBLE(x**INT(-2.D0))*(6.D0*MW2*x + 6.D0*MZ2*x + 3.D0*(MW2 - 1.D0*MZ2 + x)*A0(MW&
  &2) + 3.D0*(-1.D0*MW2 + MZ2 + x)*A0(MZ2) + 6.D0*MW2*MZ2*B0(x, MW2, MZ2) + 6.D0*MW2*x*B0(x, MW2, MZ2) + 6.D0*MZ2*x*B0(x, MW2, MZ&
  &2) - 3.D0*B0(x, MW2, MZ2)*DBLE(MW**INT(4.D0)) - 3.D0*B0(x, MW2, MZ2)*DBLE(MZ**INT(4.D0)) - 2.D0*DBLE(x**INT(2.D0)) - 3.D0*B0(x&
  &, MW2, MZ2)*DBLE(x**INT(2.D0))))/(PI2*SW2) + (0.001736111111111111D0*CW2*EL2*(6.D0*MW2 + 6.D0*MZ2 - 4.D0*x + 3.D0*A0(MW2) + 3.&
  &D0*A0(MZ2) + 6.D0*MW2*B0(x, MW2, MZ2) + 6.D0*MZ2*B0(x, MW2, MZ2) - 6.D0*x*B0(x, MW2, MZ2) + 6.D0*MW2*MZ2*DB0(x, MW2, MZ2) + 6.&
  &D0*MW2*x*DB0(x, MW2, MZ2) + 6.D0*MZ2*x*DB0(x, MW2, MZ2) - 3.D0*DB0(x, MW2, MZ2)*DBLE(MW**INT(4.D0)) - 3.D0*DB0(x, MW2, MZ2)*DB&
  &LE(MZ**INT(4.D0)) - 3.D0*DB0(x, MW2, MZ2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(32) = (0.003472222222222222D0*EL2*DBLE(x**INT(-2.D0))*(2.D0*x*(-3.D0*MW2 + x) + 3.D0*(5.D0*MW2 + 11.D0*x)*A0(MW2) + B&
  &0(x, 0.D0, MW2)*(48.D0*MW2*x - 15.D0*DBLE(MW**INT(4.D0)) + 57.D0*DBLE(x**INT(2.D0)))))/PI2 - (0.003472222222222222D0*EL2*(2.D0&
  &*x + 2.D0*(-3.D0*MW2 + x) + 33.D0*A0(MW2) + (48.D0*MW2 + 114.D0*x)*B0(x, 0.D0, MW2) + DB0(x, 0.D0, MW2)*(48.D0*MW2*x - 15.D0*D&
  &BLE(MW**INT(4.D0)) + 57.D0*DBLE(x**INT(2.D0)))))/(PI2*x)

 amplitudes(33) = (0.003472222222222222D0*CW2*EL2*DBLE(x**INT(-2.D0))*(-6.D0*MW2*x - 6.D0*MZ2*x + 3.D0*(5.D0*MW2 - 5.D0*MZ2 + 11.&
  &D0*x)*A0(MW2) + (-15.D0*MW2 + 15.D0*MZ2 + 33.D0*x)*A0(MZ2) + 30.D0*MW2*MZ2*B0(x, MW2, MZ2) + 48.D0*MW2*x*B0(x, MW2, MZ2) + 48.&
  &D0*MZ2*x*B0(x, MW2, MZ2) - 15.D0*B0(x, MW2, MZ2)*DBLE(MW**INT(4.D0)) - 15.D0*B0(x, MW2, MZ2)*DBLE(MZ**INT(4.D0)) + 2.D0*DBLE(x&
  &**INT(2.D0)) + 57.D0*B0(x, MW2, MZ2)*DBLE(x**INT(2.D0))))/(PI2*SW2) - (0.003472222222222222D0*CW2*EL2*(-6.D0*MW2 - 6.D0*MZ2 + &
  &4.D0*x + 33.D0*A0(MW2) + 33.D0*A0(MZ2) + 48.D0*MW2*B0(x, MW2, MZ2) + 48.D0*MZ2*B0(x, MW2, MZ2) + 114.D0*x*B0(x, MW2, MZ2) + 30&
  &.D0*MW2*MZ2*DB0(x, MW2, MZ2) + 48.D0*MW2*x*DB0(x, MW2, MZ2) + 48.D0*MZ2*x*DB0(x, MW2, MZ2) - 15.D0*DB0(x, MW2, MZ2)*DBLE(MW**I&
  &NT(4.D0)) - 15.D0*DB0(x, MW2, MZ2)*DBLE(MZ**INT(4.D0)) + 57.D0*DB0(x, MW2, MZ2)*DBLE(x**INT(2.D0))))/(PI2*SW2*x)

 amplitudes(34) = (0.0625D0*EL2*MW2*DB0(x, 0.D0, MW2))/PI2

 amplitudes(35) = (0.0625D0*EL2*MW2*SW2*DB0(x, MW2, MZ2))/(CW2*PI2)

 amplitudes(36) = (0.0625D0*EL2*MW2*SBA2*DB0(x, Mh02, MW2))/(PI2*SW2)

 amplitudes(37) = (0.0625D0*CBA2*EL2*MW2*DB0(x, MHH2, MW2))/(PI2*SW2)

  totalAmplitude = (0D0,0D0)
 do j=1,37
  totalAmplitude = totalAmplitude + amplitudes(j)
 end do
 DSelfWpWp = totalAmplitude
end function DSelfWpWp

