!=========================================================================
! This file is part of MOLGW.
! Authors: Fabien Bruneval
!
! This module contains
!   evaluation of the solid harmonics for a contracted gaussian basis function at point Δr
!
! NOTES:
!
! 1. Mathematica code to generate the solid harmonics
!
! Gamma[l_, m_, k_] :=
!   (-1)^k 2^(-l) Binomial[l, k] Binomial[2 l - 2 k, l] (l - 2 k)!/(l - 2 k - m)!;
! Pi[l_, m_] :=
!   Sum[gamma[l, m, k] (x^2 + y^2 + z^2)^k z^(l - 2 k - m), {k, 0, Floor[(l - m)/2]}];
! A[m_] :=
!   Sum[Binomial[m, p] x^p y^(m - p) Cos[(m - p) Pi/2], {p, 0, m}];
! B[m_] :=
!   Sum[Binomial[m, p] x^p y^(m - p) Sin[(m - p) Pi/2], {p, 0, m}];
! C[l_, m_] :=
!   Sqrt[(2 - KroneckerDelta[m, 0]) (l - m)!/(l + m)!] pi[l, m] a[m];
! S[l_, m_] :=
!   Sqrt[2 (l - m)!/(l + m)!] pi[l, m] b[m];
! Genfun[l_, m_] := If[m >= 0, If[m <= l, c[l, m], 0], If[m >= -l, s[l, -m], 0]] // Simplify;
! Result =
!  Table[Table[genfun[l, m], {m, {0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5}}], {l, 0, 5}] // TableForm
!
! 2. shells are ordered as [-m, ... m]
!
!=========================================================================
#include "molgw.h"
#if !defined(NO_LIBINT)
#include<libint2/libint2_params.h>
#endif




module m_solid_harmonics
  use m_definitions
  implicit none


contains


!=========================================================================
pure function eval_solid_harmonics(dr, am, alpha, c) result(shell)
  integer, intent(in) :: am
  real(dp), intent(in) :: dr(3), alpha(:), c(:)
  real(dp), allocatable :: shell(:)
  !=====
  !=====

  allocate(shell(2 * am + 1))

  select case(am)
  case(0)
    shell(:) = eval_l0(dr, alpha, c)
  case(1)
    shell(:) = eval_l1(dr, alpha, c)
  case(2)
    shell(:) = eval_l2(dr, alpha, c)
  case(3)
    shell(:) = eval_l3(dr, alpha, c)
  case(4)
    shell(:) = eval_l4(dr, alpha, c)
  case(5)
    shell(:) = eval_l5(dr, alpha, c)
  case default
    !eval_solid_harmonics: angular momentum l > 5 not implemented'
    shell(:) = 0.0_dp
  end select


end function eval_solid_harmonics


!=========================================================================
pure function eval_l0(dr, alpha, c) result(shell)
  real(dp), intent(in) :: dr(3), alpha(:), c(:)
  real(dp) :: shell(1)
  !=====
  real(dp) :: x, y, z
  real(dp) :: r2, expo
  integer :: i, n
  !=====

  n = SIZE(alpha)
  x = dr(1)
  y = dr(2)
  z = dr(3)
  r2 = x*x + y*y + z*z
  expo = 0.0_dp
  do i=1, n
    expo = expo + EXP(-alpha(i) * r2) * c(i)
  enddo
  shell(1) = expo

end function eval_l0


!=========================================================================
pure function eval_l1(dr, alpha, c) result(shell)
  real(dp), intent(in) :: dr(3), alpha(:), c(:)
  real(dp) :: shell(3)
  !=====
  real(dp) :: x, y, z
  real(dp) :: r2, expo
  integer :: i, n
  !=====

  n = SIZE(alpha)
  x = dr(1)
  y = dr(2)
  z = dr(3)
  r2 = x*x + y*y + z*z
  expo = 0.0_dp
  do i=1, n
    expo = expo + EXP(-alpha(i) * r2) * c(i)
  enddo

  shell(1) = y * expo
  shell(2) = z * expo
  shell(3) = x * expo

end function eval_l1


!=========================================================================
pure function eval_l2(dr, alpha, c) result(shell)
  real(dp), intent(in) :: dr(3), alpha(:), c(:)
  real(dp) :: shell(5)
  !=====
  real(dp) :: x, y, z
  real(dp) :: r2, expo
  real(dp) :: a0
  integer :: i, n
  !=====

  n = SIZE(alpha)
  x = dr(1)
  y = dr(2)
  z = dr(3)
  r2 = x*x + y*y + z*z
  expo = 0.0_dp
  do i=1, n
    expo = expo + EXP(-alpha(i) * r2) * c(i)
  enddo

  a0  = SQRT(1.0_dp/12.0_dp)

  shell(1) =         x*y         * expo
  shell(2) =         y*z         * expo
  shell(3) = a0     * (3*z*z - r2) * expo
  shell(4) =         x*z         * expo
  shell(5) = 0.5_dp * (x*x - y*y) * expo

end function eval_l2


!=========================================================================
pure function eval_l3(dr, alpha, c) result(shell)
  real(dp), intent(in) :: dr(3), alpha(:), c(:)
  real(dp) :: shell(7)
  !=====
  real(dp) :: x, y, z
  real(dp) :: r2, expo, N3
  real(dp) :: a0, a1, a2, a3
  integer  :: i, n
  !=====

  n = SIZE(alpha)
  x = dr(1)
  y = dr(2)
  z = dr(3)
  r2 = x*x + y*y + z*z
  expo = 0.0_dp
  do i=1, n
    expo = expo + EXP(-alpha(i) * r2) * c(i)
  enddo

  a0  = SQRT(1.0_dp/60.0_dp)
  a1  = SQRT(1.0_dp/40.0_dp)
  a2  = 0.5_dp
  a3  = SQRT(1.0_dp/24.0_dp)

  shell(1) = a3  * (3*x*x*y - y*y*y)       * expo
  shell(2) = a2  * 2*z*x*y                 * expo
  shell(3) = a1  * y*(5*z*z - r2)          * expo
  shell(4) = a0  * (5*z*z*z - 3*z*r2)      * expo
  shell(5) = a1  * x*(5*z*z - r2)          * expo
  shell(6) = a2  * z*(x*x - y*y)           * expo
  shell(7) = a3  * (x*x*x - 3*x*y*y)       * expo

end function eval_l3


!=========================================================================
pure function eval_l4(dr, alpha, c) result(shell)
  real(dp), intent(in) :: dr(3), alpha(:), c(:)
  real(8) :: shell(9)
  !=====
  real(dp) :: x, y, z
  real(8) :: r2, expo
  real(8) :: a0, a1, a2, a3, a4
  integer  :: i, n
  !=====

  n = SIZE(alpha)
  x = dr(1)
  y = dr(2)
  z = dr(3)
  r2 = x*x + y*y + z*z
  expo = 0.0_dp
  do i=1, n
    expo = expo + EXP(-alpha(i) * r2) * c(i)
  enddo

  a0 = SQRT(1.0_dp/6720)
  a1 = SQRT(1.0_dp/168)
  a2 = SQRT(1.0_dp/336)
  a3 = SQRT(1.0_dp/24)
  a4 = SQRT(1.0_dp/192)


  shell(1) = a4 * (4*x*x*x*y - 4*x*y*y*y) * expo
  shell(2) = a3 * z*(3*x*x*y - y*y*y) * expo
  shell(3) = a2 * (2*x*y) *(7*z*z - r2) * expo
  shell(4) = a1 * z*(7*z*z - 3*r2)*y * expo
  shell(5) = a0 * (35*z**4 - 30*z*z*r2 + 3*r2*r2) * expo
  shell(6) = a1 * z*(7*z*z - 3*r2)*x * expo
  shell(7) = a2 * (x*x - y*y)*(7*z*z - r2) * expo
  shell(8) = a3 * z*(x*x*x - 3*x*y*y) * expo
  shell(9) = a4 * (x**4 - 6*x*x*y*y + y**4) * expo

end function eval_l4


!=========================================================================
pure function eval_l5(dr, alpha, c) result(shell)
  real(dp), intent(in) :: dr(3), alpha(:), c(:)
  real(8) :: shell(11)
  !=====
  real(dp) :: x, y, z
  real(8) :: r2, expo
  real(8) :: a0, a1, a2, a3, a4, a5
  integer  :: i, n
  !=====

  n = SIZE(alpha)
  x = dr(1)
  y = dr(2)
  z = dr(3)
  r2 = x*x + y*y + z*z
  expo = 0.0_dp
  do i=1, n
    expo = expo + EXP(-alpha(i) * r2) * c(i)
  enddo

  a0 = SQRT(1.0_dp/60480)
  a1 = SQRT(1.0_dp/36288)
  a2 = SQRT(1.0_dp/144)
  a3 = SQRT(1.0_dp/3456)
  a4 = SQRT(1.0_dp/192)
  a5 = SQRT(1.0_dp/1920)

  shell(1)  = a5 * (5*x*x*x*x*y - 10*x*x*y*y*y + y**5) * expo
  shell(2)  = a4 * z*(4*x*x*x*y - 4*x*y*y*y) * expo
  shell(3)  = a3 * (3*x*x*y - y*y*y)*(9*z*z - r2) * expo
  shell(4)  = a2 * z*(2*x*y) *(3*z*z - r2) * expo
  shell(5)  = a1 * y*(63*z**4 - 42*z*z*r2 + 3*r2*r2) * expo
  shell(6)  = a0 * z*(63*z**4 - 70*z*z*r2 + 15*r2*r2) * expo
  shell(7)  = a1 * x*(63*z**4 - 42*z*z*r2 + 3*r2*r2) * expo
  shell(8)  = a2 * z*(x*x - y*y)*(3*z*z - r2) * expo
  shell(9)  = a3 * (x*x*x - 3*x*y*y)*(9*z*z - r2) * expo
  shell(10) = a4 * z*(x**4 - 6*x*x*y*y + y**4) * expo
  shell(11) = a5 * (x**5 - 10*x*x*x*y*y + 5*x*y*y*y*y) * expo

end function eval_l5


end module m_solid_harmonics


!=========================================================================
