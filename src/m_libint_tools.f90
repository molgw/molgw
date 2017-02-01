!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to interface with the C++ library LIBINT
!
!=========================================================================
module m_libint_tools
 use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
 use m_definitions
 use m_basis_set



!=========================================================================
! interfaces to libint wrappers in C/C++ with iso_c_binding
!=========================================================================

 interface

#ifdef HAVE_LIBINT_ONEBODY
   subroutine libint_overlap(amA,contrdepthA,A,alphaA,cA, &
                             amB,contrdepthB,B,alphaB,cB, &
                             overlapAB) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value  :: amA,contrdepthA
     real(C_DOUBLE),intent(in) :: A(*)
     real(C_DOUBLE),intent(in) :: alphaA(*)
     real(C_DOUBLE),intent(in) :: cA(*)
     integer(C_INT),value  :: amB,contrdepthB
     real(C_DOUBLE),intent(in) :: B(*)
     real(C_DOUBLE),intent(in) :: alphaB(*)
     real(C_DOUBLE),intent(in) :: cB(*)
     real(C_DOUBLE),intent(out) :: overlapAB(*)
     
   end subroutine libint_overlap

   subroutine libint_overlap_grad(amA,contrdepthA,A,alphaA,cA, &
                             amB,contrdepthB,B,alphaB,cB, &
                             overlapABx,overlapABy,overlapABz) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value  :: amA,contrdepthA
     real(C_DOUBLE),intent(in) :: A(*)
     real(C_DOUBLE),intent(in) :: alphaA(*)
     real(C_DOUBLE),intent(in) :: cA(*)
     integer(C_INT),value  :: amB,contrdepthB
     real(C_DOUBLE),intent(in) :: B(*)
     real(C_DOUBLE),intent(in) :: alphaB(*)
     real(C_DOUBLE),intent(in) :: cB(*)
     real(C_DOUBLE),intent(out) :: overlapABx(*)
     real(C_DOUBLE),intent(out) :: overlapABy(*)
     real(C_DOUBLE),intent(out) :: overlapABz(*)
     
   end subroutine libint_overlap_grad

   subroutine libint_kinetic(amA,contrdepthA,A,alphaA,cA, &
                             amB,contrdepthB,B,alphaB,cB, &
                             kineticAB) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value       :: amA,contrdepthA
     real(C_DOUBLE),intent(in)  :: A(*)
     real(C_DOUBLE),intent(in)  :: alphaA(*)
     real(C_DOUBLE),intent(in)  :: cA(*)
     integer(C_INT),value       :: amB,contrdepthB
     real(C_DOUBLE),intent(in)  :: B(*)
     real(C_DOUBLE),intent(in)  :: alphaB(*)
     real(C_DOUBLE),intent(in)  :: cB(*)
     real(C_DOUBLE),intent(out) :: kineticAB(*)

   end subroutine libint_kinetic

   subroutine libint_kinetic_grad(amA,contrdepthA,A,alphaA,cA, &
                                  amB,contrdepthB,B,alphaB,cB, &
                                  kineticABx,kineticABy,kineticABz) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value       :: amA,contrdepthA
     real(C_DOUBLE),intent(in)  :: A(*)
     real(C_DOUBLE),intent(in)  :: alphaA(*)
     real(C_DOUBLE),intent(in)  :: cA(*)
     integer(C_INT),value       :: amB,contrdepthB
     real(C_DOUBLE),intent(in)  :: B(*)
     real(C_DOUBLE),intent(in)  :: alphaB(*)
     real(C_DOUBLE),intent(in)  :: cB(*)
     real(C_DOUBLE),intent(out) :: kineticABx(*)
     real(C_DOUBLE),intent(out) :: kineticABy(*)
     real(C_DOUBLE),intent(out) :: kineticABz(*)

   end subroutine libint_kinetic_grad

   subroutine libint_elecpot(amA,contrdepthA,A,alphaA,cA, &
                             amB,contrdepthB,B,alphaB,cB, &
                             C,elecpotAB) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value         :: amA,contrdepthA
     real(C_DOUBLE),intent(in)    :: A(*)
     real(C_DOUBLE),intent(in)    :: alphaA(*)
     real(C_DOUBLE),intent(in)    :: cA(*)
     integer(C_INT),value         :: amB,contrdepthB
     real(C_DOUBLE),intent(in)    :: B(*)
     real(C_DOUBLE),intent(in)    :: alphaB(*)
     real(C_DOUBLE),intent(in)    :: cB(*)
     real(C_DOUBLE),intent(in)    :: C(*)
     real(C_DOUBLE),intent(inout) :: elecpotAB(*)
     
   end subroutine libint_elecpot

   subroutine libint_elecpot_grad(amA,contrdepthA,A,alphaA,cA, &
                                  amB,contrdepthB,B,alphaB,cB, &
                                  C,elecpotAx,elecpotAy,elecpotAz, &
                                    elecpotBx,elecpotBy,elecpotBz) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value         :: amA,contrdepthA
     real(C_DOUBLE),intent(in)    :: A(*)
     real(C_DOUBLE),intent(in)    :: alphaA(*)
     real(C_DOUBLE),intent(in)    :: cA(*)
     integer(C_INT),value         :: amB,contrdepthB
     real(C_DOUBLE),intent(in)    :: B(*)
     real(C_DOUBLE),intent(in)    :: alphaB(*)
     real(C_DOUBLE),intent(in)    :: cB(*)
     real(C_DOUBLE),intent(in)    :: C(*)
     real(C_DOUBLE),intent(inout) :: elecpotAx(*)
     real(C_DOUBLE),intent(inout) :: elecpotAy(*)
     real(C_DOUBLE),intent(inout) :: elecpotAz(*)
     real(C_DOUBLE),intent(inout) :: elecpotBx(*)
     real(C_DOUBLE),intent(inout) :: elecpotBy(*)
     real(C_DOUBLE),intent(inout) :: elecpotBz(*)
     
   end subroutine libint_elecpot_grad
#endif

#ifdef HAVE_LIBINT_2CENTER
   subroutine libint_2center(amA,contrdepthA,A,alphaA,cA, &
                             amC,contrdepthC,C,alphaC,cC, &
                             rcut,eriAC) bind(C)
     use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE
     integer(C_INT),value         :: amA,contrdepthA
     real(C_DOUBLE),intent(in)    :: A(*)
     real(C_DOUBLE),intent(in)    :: alphaA(*)
     real(C_DOUBLE),intent(in)    :: cA(*)
     integer(C_INT),value         :: amC,contrdepthC
     real(C_DOUBLE),intent(in)    :: C(*)
     real(C_DOUBLE),intent(in)    :: alphaC(*)
     real(C_DOUBLE),intent(in)    :: cC(*)
     real(C_DOUBLE),intent(in),value :: rcut
     real(C_DOUBLE),intent(inout)    :: eriAC(*)

   end subroutine libint_2center
#endif

 end interface



 interface transform_libint_to_molgw
   module procedure transform_libint_to_molgw_2d
!   module procedure transform_libint_to_molgw_3d
!   module procedure transform_libint_to_molgw_4d
 end interface



contains


!=========================================================================
subroutine transform_libint_to_molgw_2d(gaussian_type,am1,am2,array_in,matrix_out)
 implicit none
 character(len=4),intent(in)      :: gaussian_type
 integer,intent(in)               :: am1,am2
 real(C_DOUBLE),intent(in)        :: array_in(:)
 real(dp),allocatable,intent(out) :: matrix_out(:,:)
!=====
 integer :: n1,n2,n1c,n2c
 integer :: i1c,i2c,i12c
 real(dp),allocatable :: matrix_tmp(:,:)
!=====

 n1c = number_basis_function_am('CART',am1)
 n2c = number_basis_function_am('CART',am2)
 n1  = number_basis_function_am(gaussian_type,am1)
 n2  = number_basis_function_am(gaussian_type,am2)

 if( .NOT. ALLOCATED(matrix_out) ) allocate(matrix_out(n1,n2))
 allocate(matrix_tmp(n1,n2c))

 matrix_tmp(:,:) = TRANSPOSE( MATMUL( RESHAPE( array_in(:) , (/ n2c , n1c /) ) , cart_to_pure_norm(am1)%matrix(1:n1c,1:n1) ) )

 matrix_out(:,:) = MATMUL( matrix_tmp(:,:) , cart_to_pure_norm(am2)%matrix(:,:) )

 deallocate(matrix_tmp)

end subroutine transform_libint_to_molgw_2d


!=========================================================================
subroutine transform_libint_to_molgw_3d(gaussian_type_left,am1,gaussian_type_right,am2,am3,array_in,matrix_out)
 implicit none
 character(len=4),intent(in)      :: gaussian_type_left,gaussian_type_right
 integer,intent(in)               :: am1,am2,am3
 real(C_DOUBLE),intent(in)        :: array_in(:)
 real(dp),allocatable,intent(out) :: matrix_out(:,:,:)
!=====
 integer :: n1,n2,n3,n1c,n2c,n3c
 integer :: i1c,i2c,i12c
 real(dp),allocatable :: matrix_tmp(:,:,:)
!=====

 n1c = number_basis_function_am('CART',am1)
 n2c = number_basis_function_am('CART',am2)
 n3c = number_basis_function_am('CART',am3)
 n1  = number_basis_function_am(gaussian_type_left,am1)
 n2  = number_basis_function_am(gaussian_type_right,am2)
 n3  = number_basis_function_am(gaussian_type_right,am3)

 if( am2 + am3 >= am1 ) then


 else

 endif



! if( .NOT. ALLOCATED(matrix_out) ) allocate(matrix_out(n1,n2,n3))
! allocate(matrix_tmp(n1,n2c))
!
! matrix_tmp(:,:) = TRANSPOSE( MATMUL( RESHAPE( array_in(:) , (/ n2c , n1c /) ) , cart_to_pure_norm(am1)%matrix(1:n1c,1:n1) ) )
!
! matrix_out(:,:) = MATMUL( matrix_tmp(:,:) , cart_to_pure_norm(am2)%matrix(:,:) )
!
! deallocate(matrix_tmp)

end subroutine transform_libint_to_molgw_3d




end module m_libint_tools


!=========================================================================
