!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval (Fortran version of the LIBINT example)
!
! This file contains
!  the evaluation of the Boys function used for Coulomb integrals
!
!=========================================================================
subroutine boys_function_c(fnt,nn,tt)  BIND(C)
  use,intrinsic :: ISO_C_BINDING,only: C_INT,C_DOUBLE
  implicit none
  integer(C_INT),value       :: nn
  real(C_DOUBLE),value       :: tt
  real(C_DOUBLE),intent(out) :: fnt(0:nn)
  !=====
  integer(C_INT),parameter  :: maxfac=100
  real(C_DOUBLE),parameter  :: eps=1.0e-17_C_DOUBLE
  integer(C_INT)            :: ii,mm
  integer(C_INT)            :: m2
  real(C_DOUBLE)            :: t2,num,sum,term1,et,t_sqrt
  real(C_DOUBLE),parameter  :: kk = 0.8862269254527579_C_DOUBLE  ! 0.5 * sqrt(pi)
  real(C_DOUBLE),save       :: df(2*maxfac)=0.0_C_DOUBLE
  !=====
 
  if( ABS(df(1)) < 1.0e-10_C_DOUBLE ) then
    df(1:3) = 1.0_C_DOUBLE
    do ii=4,2*maxfac
      df(ii) = (ii-2) * df(ii-2)
    enddo
  endif
 
  et = EXP(-tt)
  t2 = 2.0_C_DOUBLE * tt
 
  if( tt > 20.0_C_DOUBLE ) then ! For big tt's do upward recursion
 
    t_sqrt  = SQRT(tt)
    fnt(0) = kk * ERF(t_sqrt) / t_sqrt
    do mm=0,nn-1
      fnt(mm+1) = ( (2*mm+1) * fnt(mm) - et ) / t2
    enddo
 
  else
    !   For smaller tt's compute F with highest nn using
    !   asymptotic series (see I. Shavitt in
    !   Methods in Computational Physics, ed. B. Alder eta l,
    !   vol 2, 1963, page 8)
 
    m2 = 2 * nn
    num = df(m2+1)
    sum = 1.0_C_DOUBLE / ( m2 + 1 )
    do ii=1,maxfac-1
      num = num * t2
      term1 = num / df(m2 + 2*ii + 3)
      sum = sum + term1
      if( ABS(term1) < eps ) exit
    enddo
    fnt(nn) = sum * et
    !
    ! And then do downward recursion
    do mm=nn-1,0,-1
      fnt(mm)= ( t2 * fnt(mm+1) + et ) / ( 2 * mm + 1 )
    enddo
 
  endif


end subroutine boys_function_c


!=========================================================================
