!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval (Fortran version of the LIBINT example)
!
! This file contains
!  the evaluation of the Boys function used for Coulomb integrals
!
!=========================================================================
subroutine boys_function_c(fnt,n,t)  BIND(C)
 use,intrinsic :: ISO_C_BINDING,only: C_INT,C_DOUBLE
 implicit none
 integer(C_INT),value       :: n
 real(C_DOUBLE),value       :: t
 real(C_DOUBLE),intent(out) :: fnt(0:n)
!=====
 integer(C_INT),parameter  :: maxfac=100
 real(C_DOUBLE),parameter  :: eps=1.0e-17_C_DOUBLE
 integer(C_INT)            :: i,m,k
 integer(C_INT)            :: m2
 real(C_DOUBLE)            :: t2,num,sum,term1,term2,et,tt
 real(C_DOUBLE),parameter  :: kk = 0.8862269254527579_C_DOUBLE  ! 0.5 * sqrt(pi)
 real(C_DOUBLE),save       :: df(2*maxfac)=0.0_C_DOUBLE
!=====

 if( ABS(df(1)) < 1.0e-10_C_DOUBLE ) then
   df(1:3) = 1.0_C_DOUBLE
   do i=4,2*maxfac
     df(i) = (i-2) * df(i-2)
   enddo
 endif

 if( t > 20.0_C_DOUBLE ) then ! For big t's do upward recursion 
   t2 = 2.0_C_DOUBLE * t
   et = EXP(-t)
   tt  = SQRT(t)
   fnt(0) = kk * ERF(tt) / tt
   do m=0,n-1
     fnt(m+1) = ( (2*m+1) * fnt(m) - et ) / t2
   enddo

 else
   !   For smaller t's compute F with highest n using
   !   asymptotic series (see I. Shavitt in
   !   Methods in Computational Physics, ed. B. Alder eta l,
   !   vol 2, 1963, page 8)

   et = EXP(-t)
   t2 = 2.0_C_DOUBLE * t
   m2 = 2 * n
   num = df(m2+1)
   sum = 1.0_C_DOUBLE / ( m2 + 1 )
   do i=1,maxfac-1
     num = num * t2
     term1 = num / df(m2 + 2*i + 3)
     sum = sum + term1
     if( ABS(term1) < eps ) exit
   enddo
   fnt(n) = sum * et 
   !
   ! And then do downward recursion 
   do m=n-1,0,-1
     fnt(m)= ( t2 * fnt(m+1) + et ) / ( 2 * m + 1 )
   enddo

 endif


end subroutine boys_function_c


!=========================================================================
