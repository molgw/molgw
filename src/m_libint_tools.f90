!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to interface with the C++ library LIBINT
! * explicit interfaces with C binding
! * subroutines to pass from row-major (C convention) to column-major (Fortran convention)
!
!=========================================================================
#include "molgw.h"
#if !defined(NO_LIBINT)
#include<libint2/libint2_params.h>
#endif

module m_libint_tools
  use m_definitions
  use m_cart_to_pure
  use m_basis_set


  !=========================================================================
  ! interfaces to libint wrappers in C/C++ with iso_c_binding
  !=========================================================================
  interface

#if (LIBINT2_SUPPORT_ONEBODY)
    subroutine libint_overlap(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              overlapAB) bind(C)
      import :: C_INT,C_DOUBLE
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

    subroutine libint_kinetic(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              kineticAB) bind(C)
      import :: C_INT,C_DOUBLE
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

    subroutine libint_elecpot(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              C,elecpotAB) bind(C)
      import :: C_INT,C_DOUBLE
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
#endif

#if (LIBINT2_DERIV_ONEBODY_ORDER > 0)
    subroutine libint_overlap_grad(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              overlapABx,overlapABy,overlapABz) bind(C)
      import :: C_INT,C_DOUBLE
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

    subroutine libint_kinetic_grad(amA,contrdepthA,A,alphaA,cA, &
                                   amB,contrdepthB,B,alphaB,cB, &
                                   kineticABx,kineticABy,kineticABz) bind(C)
      import :: C_INT,C_DOUBLE
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

    subroutine libint_elecpot_grad(amA,contrdepthA,A,alphaA,cA, &
                                   amB,contrdepthB,B,alphaB,cB, &
                                   C,elecpotAx,elecpotAy,elecpotAz, &
                                     elecpotBx,elecpotBy,elecpotBz) bind(C)
      import :: C_INT,C_DOUBLE
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
#if (LIBINT2_DERIV_ERI_ORDER >0)
    subroutine libint_4center_grad(amA,contrdepthA,A,alphaA,cA, &
                                   amB,contrdepthB,B,alphaB,cB, &
                                   amC,contrdepthC,C,alphaC,cC, &
                                   amD,contrdepthD,D,alphaD,cD, &
                                   rcut, &
                                   eriAx,eriAy,eriAz, &
                                   eriBx,eriBy,eriBz, &
                                   eriCx,eriCy,eriCz, &
                                   eriDx,eriDy,eriDz) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value         :: amA,contrdepthA
      real(C_DOUBLE),intent(in)    :: A(*)
      real(C_DOUBLE),intent(in)    :: alphaA(*)
      real(C_DOUBLE),intent(in)    :: cA(*)
      integer(C_INT),value         :: amB,contrdepthB
      real(C_DOUBLE),intent(in)    :: B(*)
      real(C_DOUBLE),intent(in)    :: alphaB(*)
      real(C_DOUBLE),intent(in)    :: cB(*)
      integer(C_INT),value         :: amC,contrdepthC
      real(C_DOUBLE),intent(in)    :: C(*)
      real(C_DOUBLE),intent(in)    :: alphaC(*)
      real(C_DOUBLE),intent(in)    :: cC(*)
      integer(C_INT),value         :: amD,contrdepthD
      real(C_DOUBLE),intent(in)    :: D(*)
      real(C_DOUBLE),intent(in)    :: alphaD(*)
      real(C_DOUBLE),intent(in)    :: cD(*)
      real(C_DOUBLE),intent(in),value :: rcut
      real(C_DOUBLE),intent(inout) :: eriAx(*)
      real(C_DOUBLE),intent(inout) :: eriAy(*)
      real(C_DOUBLE),intent(inout) :: eriAz(*)
      real(C_DOUBLE),intent(inout) :: eriBx(*)
      real(C_DOUBLE),intent(inout) :: eriBy(*)
      real(C_DOUBLE),intent(inout) :: eriBz(*)
      real(C_DOUBLE),intent(inout) :: eriCx(*)
      real(C_DOUBLE),intent(inout) :: eriCy(*)
      real(C_DOUBLE),intent(inout) :: eriCz(*)
      real(C_DOUBLE),intent(inout) :: eriDx(*)
      real(C_DOUBLE),intent(inout) :: eriDy(*)
      real(C_DOUBLE),intent(inout) :: eriDz(*)

    end subroutine libint_4center_grad
#endif

    subroutine libint_init(ammax,has_onebody,has_gradient) bind(C)
      import :: C_INT,C_BOOL
      integer(C_INT),intent(out)  :: ammax
      logical(C_BOOL),intent(out) :: has_onebody,has_gradient
    end subroutine libint_init

    subroutine libint_2center(amA,contrdepthA,A,alphaA,cA, &
                              amC,contrdepthC,C,alphaC,cC, &
                              rcut,eriAC) bind(C)
      import ::  C_INT,C_DOUBLE
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

    subroutine libint_3center(amA,contrdepthA,A,alphaA,cA, &
                              amC,contrdepthC,C,alphaC,cC, &
                              amD,contrdepthD,D,alphaD,cD, &
                              rcut,eriACD) bind(C)
      import ::  C_INT,C_DOUBLE
      integer(C_INT),value         :: amA,contrdepthA
      real(C_DOUBLE),intent(in)    :: A(*)
      real(C_DOUBLE),intent(in)    :: alphaA(*)
      real(C_DOUBLE),intent(in)    :: cA(*)
      integer(C_INT),value         :: amC,contrdepthC
      real(C_DOUBLE),intent(in)    :: C(*)
      real(C_DOUBLE),intent(in)    :: alphaC(*)
      real(C_DOUBLE),intent(in)    :: cC(*)
      integer(C_INT),value         :: amD,contrdepthD
      real(C_DOUBLE),intent(in)    :: D(*)
      real(C_DOUBLE),intent(in)    :: alphaD(*)
      real(C_DOUBLE),intent(in)    :: cD(*)
      real(C_DOUBLE),intent(in),value :: rcut
      real(C_DOUBLE),intent(inout)    :: eriACD(*)

    end subroutine libint_3center


    subroutine libint_4center(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              amC,contrdepthC,C,alphaC,cC, &
                              amD,contrdepthD,D,alphaD,cD, &
                              rcut,eriABCD) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value         :: amA,contrdepthA
      real(C_DOUBLE),intent(in)    :: A(*)
      real(C_DOUBLE),intent(in)    :: alphaA(*)
      real(C_DOUBLE),intent(in)    :: cA(*)
      integer(C_INT),value         :: amB,contrdepthB
      real(C_DOUBLE),intent(in)    :: B(*)
      real(C_DOUBLE),intent(in)    :: alphaB(*)
      real(C_DOUBLE),intent(in)    :: cB(*)
      integer(C_INT),value         :: amC,contrdepthC
      real(C_DOUBLE),intent(in)    :: C(*)
      real(C_DOUBLE),intent(in)    :: alphaC(*)
      real(C_DOUBLE),intent(in)    :: cC(*)
      integer(C_INT),value         :: amD,contrdepthD
      real(C_DOUBLE),intent(in)    :: D(*)
      real(C_DOUBLE),intent(in)    :: alphaD(*)
      real(C_DOUBLE),intent(in)    :: cD(*)
      real(C_DOUBLE),intent(in),value :: rcut
      real(C_DOUBLE),intent(inout)    :: eriABCD(*)

    end subroutine libint_4center


  end interface

  interface transform_molgw_to_molgw
    module procedure transform_molgw_to_molgw_2d
  end interface


  interface transform_libint_to_molgw
    module procedure transform_libint_to_molgw_2d
    module procedure transform_libint_to_molgw_3d
    module procedure transform_libint_to_molgw_4d
  end interface



contains


!=========================================================================
subroutine transform_molgw_to_molgw_2d(gaussian_type,am1,am2,array_in,matrix_out)
  implicit none
  character(len=4),intent(in)      :: gaussian_type
  integer,intent(in)               :: am1,am2
  real(C_DOUBLE),intent(in)        :: array_in(:)
  real(dp),allocatable,intent(out) :: matrix_out(:,:)
  !=====
  integer :: n1,n2,n1c,n2c
  integer :: gt_tag
  real(dp),allocatable :: matrix_tmp(:,:)
  !=====

  gt_tag = get_gaussian_type_tag(gaussian_type)
  n1c = number_basis_function_am('CART',am1)
  n2c = number_basis_function_am('CART',am2)
  n1  = number_basis_function_am(gaussian_type,am1)
  n2  = number_basis_function_am(gaussian_type,am2)

  if( .NOT. ALLOCATED(matrix_out) ) allocate(matrix_out(n1,n2))
  allocate(matrix_tmp(n1,n2c))

  ! Transform the 1st index
  matrix_tmp(:,:) = TRANSPOSE( MATMUL( RESHAPE( array_in(:) , (/ n2c , n1c /) ) , cart_to_pure(am1,gt_tag)%matrix(1:n1c,1:n1) ) )

  ! Transform the 2nd index
  matrix_out(:,:) = MATMUL( matrix_tmp(:,:) , cart_to_pure(am2,gt_tag)%matrix(:,:) )

  deallocate(matrix_tmp)

end subroutine transform_molgw_to_molgw_2d


!=========================================================================
subroutine transform_libint_to_molgw_2d(gaussian_type,am1,am2,array_in,matrix_out)
  implicit none
  character(len=4),intent(in)      :: gaussian_type
  integer,intent(in)               :: am1,am2
  real(C_DOUBLE),intent(in)        :: array_in(:)
  real(dp),allocatable,intent(out) :: matrix_out(:,:)
  !=====
  integer :: n1,n2,n1c,n2c
  integer :: gt_tag
  real(dp),allocatable :: matrix_tmp(:,:)
  !=====

  gt_tag = get_gaussian_type_tag(gaussian_type)
  n1c = number_basis_function_am('CART',am1)
  n2c = number_basis_function_am('CART',am2)
  n1  = number_basis_function_am(gaussian_type,am1)
  n2  = number_basis_function_am(gaussian_type,am2)

  if( .NOT. ALLOCATED(matrix_out) ) allocate(matrix_out(n1,n2))
  allocate(matrix_tmp(n1,n2c))

  ! Transform the 1st index
  matrix_tmp(:,:) = TRANSPOSE( MATMUL( RESHAPE( array_in(:) , (/ n2c , n1c /) ) , &
                    cart_to_pure_norm(am1,gt_tag)%matrix(1:n1c,1:n1) ) )

  ! Transform the 2nd index
  matrix_out(:,:) = MATMUL( matrix_tmp(:,:) , cart_to_pure_norm(am2,gt_tag)%matrix(:,:) )

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
  integer :: i1,i2,i3,ii
  integer :: gt_tagl,gt_tagr
  real(dp),allocatable :: matrix_tmp1(:,:)
  real(dp),allocatable :: matrix_tmp2(:,:)
  !=====

  gt_tagl = get_gaussian_type_tag(gaussian_type_left)
  gt_tagr = get_gaussian_type_tag(gaussian_type_right)
  n1c = number_basis_function_am('CART',am1)
  n2c = number_basis_function_am('CART',am2)
  n3c = number_basis_function_am('CART',am3)
  n1  = number_basis_function_am(gaussian_type_left,am1)
  n2  = number_basis_function_am(gaussian_type_right,am2)
  n3  = number_basis_function_am(gaussian_type_right,am3)

  if( .NOT. ALLOCATED(matrix_out) ) allocate(matrix_out(n1,n2,n3))

  ! When CART -> CART, just normalize (no unitary transform needed)
  if( gt_tagl == CARTG .AND. gt_tagr == CARTG ) then
    ii = 0
    do i1=1,n1
      do i2=1,n2
        do i3=1,n3
          ii = ii + 1
          matrix_out(i1,i2,i3) = array_in(ii) * cart_to_pure_norm(am1,CARTG)%matrix(i1,i1) &
                                              * cart_to_pure_norm(am2,CARTG)%matrix(i2,i2) &
                                              * cart_to_pure_norm(am3,CARTG)%matrix(i3,i3) 
        enddo
      enddo
    enddo
  else

    allocate(matrix_tmp1(n1,n2c*n3c))
    allocate(matrix_tmp2(n2,n3c))

    ! Transform the 1st index
    matrix_tmp1(:,:) = TRANSPOSE( MATMUL( RESHAPE( array_in(:) , (/ n2c * n3c , n1c /) ) , &
                        cart_to_pure_norm(am1,gt_tagl)%matrix(1:n1c,1:n1) ) )

    do i1=1,n1
      ! Transform the 2nd index
      matrix_tmp2(1:n2,:) = TRANSPOSE( MATMUL( RESHAPE( matrix_tmp1(i1,:) , (/ n3c , n2c /) ) ,  &
                                                cart_to_pure_norm(am2,gt_tagr)%matrix(1:n2c,1:n2) ) )

      ! Transform the 3rd index
      matrix_out(i1,:,:) = MATMUL( matrix_tmp2(:,:) , cart_to_pure_norm(am3,gt_tagr)%matrix(:,:) )
    enddo

    deallocate(matrix_tmp1,matrix_tmp2)

  endif

end subroutine transform_libint_to_molgw_3d


!=========================================================================
subroutine transform_libint_to_molgw_4d(gaussian_type,am1,am2,am3,am4,array_in,matrix_out)
  implicit none
  character(len=4),intent(in)      :: gaussian_type
  integer,intent(in)               :: am1,am2,am3,am4
  real(C_DOUBLE),intent(in)        :: array_in(:)
  real(dp),allocatable,intent(out) :: matrix_out(:,:,:,:)
  !=====
  integer :: n1,n2,n3,n4,n1c,n2c,n3c,n4c
  integer :: i1,i2
  integer :: gt_tag
  real(dp),allocatable :: matrix_tmp0(:,:)
  real(dp),allocatable :: matrix_tmp1(:,:)
  real(dp),allocatable :: matrix_tmp2(:,:)
  real(dp),allocatable :: matrix_tmp3(:,:)
  !=====

  gt_tag = get_gaussian_type_tag(gaussian_type)
  n1c = number_basis_function_am('CART',am1)
  n2c = number_basis_function_am('CART',am2)
  n3c = number_basis_function_am('CART',am3)
  n4c = number_basis_function_am('CART',am4)
  n1  = number_basis_function_am(gaussian_type,am1)
  n2  = number_basis_function_am(gaussian_type,am2)
  n3  = number_basis_function_am(gaussian_type,am3)
  n4  = number_basis_function_am(gaussian_type,am4)

  if( .NOT. ALLOCATED(matrix_out) ) allocate(matrix_out(n1,n2,n3,n4))

  allocate(matrix_tmp1(n1,n2c*n3c*n4c))
  allocate(matrix_tmp2(n2,n3c*n4c))
  allocate(matrix_tmp3(n3,n4c))


  ! Transform the 1st index
  allocate(matrix_tmp0(n2c*n3c*n4c,n1c))
  matrix_tmp0(:,:) = RESHAPE( array_in(:) , (/ n2c * n3c * n4c , n1c /) )

  matrix_tmp1(1:n1,:) = TRANSPOSE( MATMUL( matrix_tmp0(:,:) , cart_to_pure_norm(am1,gt_tag)%matrix(1:n1c,1:n1) ) )
  deallocate(matrix_tmp0)

  do i1=1,n1
    ! Transform the 2nd index
    matrix_tmp2(1:n2,:) = TRANSPOSE( MATMUL( RESHAPE( matrix_tmp1(i1,:) , (/ n3c * n4c , n2c /) ) ,  &
                                              cart_to_pure_norm(am2,gt_tag)%matrix(1:n2c,1:n2) ) )
    do i2=1,n2
      ! Transform the 3rd index
      matrix_tmp3(1:n3,:) = TRANSPOSE( MATMUL( RESHAPE( matrix_tmp2(i2,:) , (/ n4c , n3c /) ) ,  &
                                              cart_to_pure_norm(am3,gt_tag)%matrix(1:n3c,1:n3) ) )

      ! Transform the 4th index
      matrix_out(i1,i2,:,:) = MATMUL( matrix_tmp3(:,:) , cart_to_pure_norm(am4,gt_tag)%matrix(:,:) )
    enddo
  enddo

  deallocate(matrix_tmp1,matrix_tmp2,matrix_tmp3)

end subroutine transform_libint_to_molgw_4d


!=========================================================================
! Here index 1 stands for the projector which is always 'PURE' and already normalized
subroutine transform_libint_to_molgw_gth_projector(gaussian_type,am1,am2,array_in,matrix_out)
  implicit none
  character(len=4),intent(in)      :: gaussian_type
  integer,intent(in)               :: am1,am2
  real(C_DOUBLE),intent(in)        :: array_in(:)
  real(dp),allocatable,intent(out) :: matrix_out(:,:)
  !=====
  integer :: n1,n2,n1c,n2c
  integer :: gt_tag
  real(dp),allocatable :: matrix_tmp(:,:)
  !=====

  gt_tag = get_gaussian_type_tag(gaussian_type)
  n1c = number_basis_function_am('CART',am1)
  n2c = number_basis_function_am('CART',am2)
  n1  = number_basis_function_am('PURE',am1)
  n2  = number_basis_function_am(gaussian_type,am2)

  if( .NOT. ALLOCATED(matrix_out) ) allocate(matrix_out(n1,n2))
  allocate(matrix_tmp(n1,n2c))

  ! Transform the 1st index
  matrix_tmp(:,:) = TRANSPOSE( MATMUL( RESHAPE( array_in(:) , (/ n2c , n1c /) ) , &
                    cart_to_pure(am1,gt_tag)%matrix(1:n1c,1:n1) ) )

  ! Transform the 2nd index
  matrix_out(:,:) = MATMUL( matrix_tmp(:,:) , cart_to_pure_norm(am2,gt_tag)%matrix(:,:) )

  deallocate(matrix_tmp)

end subroutine transform_libint_to_molgw_gth_projector


!=========================================================================
subroutine set_libint_shell(shell,amA,contrdepthA,A,alphaA,cA)
  implicit none

  type(shell_type),intent(in)            :: shell
  integer(C_INT),intent(out)             :: amA,contrdepthA
  real(C_DOUBLE),intent(out)             :: A(3)
  real(C_DOUBLE),allocatable,intent(out) :: alphaA(:),cA(:)
  !=====
  !=====

  allocate(alphaA(shell%ng))
  allocate(cA(shell%ng))

  contrdepthA = shell%ng
  amA         = shell%am
  A(:)        = shell%x0(:)
  alphaA(:)   = shell%alpha(:)
  cA(:)       = shell%coeff(:)

end subroutine set_libint_shell


!=========================================================================
end module m_libint_tools
!=========================================================================
