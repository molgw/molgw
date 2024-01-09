!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval (freely adapted from numerical recipe)
!
! This file contains
!  the evaluation of the real spherical harmonics S_lm(cos theta, phi)
!
!=========================================================================
#include "molgw.h"
function real_spherical_harmonics(ll,mm,cos_theta,phi) result(slm)
  use m_definitions
  implicit none
  integer,intent(in)  :: ll,mm
  real(dp),intent(in) :: cos_theta,phi
  real(dp)            :: slm
  !=====
  integer  :: ilm
  integer  :: abs_mm
  real(dp) :: factor,factorial
  real(dp),external :: legendre_polynomial
  !=====

  abs_mm = ABS(mm)

  factorial = 1.0_dp
  if( abs_mm > 0 ) then
    do ilm = ll-abs_mm+1,ll+abs_mm
      factorial = factorial * REAL(ilm,dp)
    enddo
  endif

  factor = SQRT( (2.0_dp * ll + 1.0_dp) / (2.0_dp * pi) / factorial )

  if( mm == 0 ) then
    slm = 1.0_dp / SQRT(2.0_dp)
  else if( mm > 0 ) then
    slm = COS( abs_mm * phi )
  else
    slm = SIN( abs_mm * phi )
  endif

  slm = slm * factor * legendre_polynomial(ll,abs_mm,cos_theta)

end function real_spherical_harmonics


!=========================================================================
function legendre_polynomial(ll,mm,xx) result(plgndr)
  use m_definitions
  implicit none
  integer,intent(in)  :: ll,mm
  real(dp),intent(in) :: xx
  real(dp)            :: plgndr
  !=====
  ! Computes the associated Legendre polynomial Pml(xx).
  ! Here mm and ll are integers satisfying
  ! 0 <= mm <= ll, while xx lies in the range −1 ≤ xx ≤ 1.
  integer  :: i,ill
  real(dp) ::  fact,pll,pmm,pmmp1,somx2
  !=====

  if( mm < 0 .OR. mm > ll .OR. ABS(xx) > 1.0 ) stop 'bad arguments in plgndr'

  pmm = 1.0
  ! Compute P mm .
  if( mm > 0) then
    somx2 = SQRT((1.0 - xx) * (1.0 + xx))

    fact  = 1.0
    do i=1,mm
      pmm=-pmm*fact*somx2
      fact=fact+2.
    enddo
  endif

  if( ll == mm ) then
    plgndr = pmm
  else
    pmmp1 = xx * (2*mm+1) * pmm ! Compute P mm mm+1.
    if(ll == mm+1) then
      plgndr=pmmp1
    else
      ! Compute P mm ll , ll>mm + 1.
      do ill=mm+2,ll
        pll = (xx*(2*ill-1) * pmmp1 - (ill+mm-1)*pmm) / ( ill - mm )
        pmm = pmmp1
        pmmp1 = pll
      enddo
      plgndr = pll
    endif
  endif

end function  legendre_polynomial


!=========================================================================
