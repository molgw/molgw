!=========================================================================
! This file is part of MOLGW.
! Authors: Fabien Bruneval, Ivan Maliyov
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian
! with no distribution of the memory
!
!=========================================================================
#include<libint2/libint2_params.h>

module m_hamiltonian_onebody
  use m_definitions
  use m_tddft_variables
  use m_timing
  use m_mpi
  use m_scalapack
  use m_warning
  use m_memory
  use m_cart_to_pure
  use m_inputparam,only: nspin,spin_fact,scalapack_block_min
  use m_basis_set
  use m_libint_tools
  use m_libcint_tools
  use m_io




contains


!=========================================================================
subroutine setup_overlap(basis,s_matrix)
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: s_matrix(basis%nbf,basis%nbf)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable :: array_cart(:)
 integer(C_INT)             :: amA,contrdepthA
 real(C_DOUBLE)             :: A(3)
 real(C_DOUBLE),allocatable :: alphaA(:)
 real(C_DOUBLE),allocatable :: cA(:)
 integer(C_INT)             :: amB,contrdepthB
 real(C_DOUBLE)             :: B(3)
 real(C_DOUBLE),allocatable :: alphaB(:)
 real(C_DOUBLE),allocatable :: cB(:)
!=====
 integer :: i_cart,j_cart,ij
 integer :: ibf_cart,jbf_cart
#if defined(HAVE_LIBCINT)
 integer :: info
 integer :: shls(2)
#endif
!=====

 call start_clock(timing_overlap)
#if defined(HAVE_LIBCINT)
 write(stdout,'(/,a)') ' Setup overlap matrix S (LIBCINT)'
#elif defined(LIBINT2_SUPPORT_ONEBODY)
 write(stdout,'(/,a)') ' Setup overlap matrix S (LIBINT)'
#else
 write(stdout,'(/,a)') ' Setup overlap matrix S (internal)'
#endif

 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=jshell,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart(ni_cart*nj_cart))

#if defined(HAVE_LIBCINT)
     shls(1) = ishell-1  ! C convention starts with 0
     shls(2) = jshell-1  ! C convention starts with 0
     info = cint1e_ovlp_cart(array_cart, shls, atm, LIBCINT_natm, bas, LIBCINT_nbas, env)

     call transform_libcint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)

#elif defined(LIBINT2_SUPPORT_ONEBODY)
     call libint_overlap(amA,contrdepthA,A,alphaA,cA, &
                         amB,contrdepthB,B,alphaB,cB, &
                         array_cart)

     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#else
     ij = 0
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         ij = ij + 1
         ibf_cart = basis%shell(ishell)%istart_cart + i_cart - 1
         jbf_cart = basis%shell(jshell)%istart_cart + j_cart - 1
         call overlap_basis_function(basis%bfc(ibf_cart),basis%bfc(jbf_cart),array_cart(ij))
       enddo
     enddo
     call transform_molgw_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#endif

     deallocate(alphaA,cA)


     s_matrix(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)
     s_matrix(jbf1:jbf2,ibf1:ibf2) = TRANSPOSE(matrix(:,:))

     deallocate(array_cart,matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 title='=== Overlap matrix S ==='
 call dump_out_matrix(.FALSE.,title,s_matrix)


 call stop_clock(timing_overlap)


end subroutine setup_overlap


!=========================================================================
subroutine setup_overlap_mixedbasis(basis1,basis2,s_matrix)
 implicit none
 type(basis_set),intent(in) :: basis1,basis2
 real(dp),intent(out)       :: s_matrix(basis1%nbf,basis2%nbf)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable :: array_cart(:)
 integer(C_INT)             :: amA,contrdepthA
 real(C_DOUBLE)             :: A(3)
 real(C_DOUBLE),allocatable :: alphaA(:)
 real(C_DOUBLE),allocatable :: cA(:)
 integer(C_INT)             :: amB,contrdepthB
 real(C_DOUBLE)             :: B(3)
 real(C_DOUBLE),allocatable :: alphaB(:)
 real(C_DOUBLE),allocatable :: cB(:)
!=====
 integer :: i_cart,j_cart,ij
 integer :: ibf_cart,jbf_cart
!=====

 call start_clock(timing_overlap)
#if defined(LIBINT2_SUPPORT_ONEBODY)
 write(stdout,'(/,a)') ' Setup mixed overlap matrix S (LIBINT)'
#else
 write(stdout,'(/,a)') ' Setup mixed overlap matrix S (internal)'
#endif

 if( basis1%gaussian_type /= basis2%gaussian_type ) call die('setup_overlap_mixedbasis_libint: case not implemented')

 do jshell=1,basis2%nshell
   lj      = basis2%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis2%gaussian_type,lj)
   jbf1    = basis2%shell(jshell)%istart
   jbf2    = basis2%shell(jshell)%iend

   call set_libint_shell(basis2%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=1,basis1%nshell
     li      = basis1%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis1%gaussian_type,li)
     ibf1    = basis1%shell(ishell)%istart
     ibf2    = basis1%shell(ishell)%iend

     call set_libint_shell(basis1%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart(ni_cart*nj_cart))

#if defined(LIBINT2_SUPPORT_ONEBODY)
     call libint_overlap(amA,contrdepthA,A,alphaA,cA, &
                         amB,contrdepthB,B,alphaB,cB, &
                         array_cart)

     call transform_libint_to_molgw(basis1%gaussian_type,li,lj,array_cart,matrix)
#else
     ij = 0
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         ij = ij + 1
         ibf_cart = basis1%shell(ishell)%istart_cart + i_cart - 1
         jbf_cart = basis2%shell(jshell)%istart_cart + j_cart - 1
         call overlap_basis_function(basis1%bfc(ibf_cart),basis2%bfc(jbf_cart),array_cart(ij))
       enddo
     enddo
     call transform_molgw_to_molgw(basis1%gaussian_type,li,lj,array_cart,matrix)
#endif

     deallocate(alphaA,cA)


     s_matrix(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)

     deallocate(array_cart,matrix)

   enddo
   deallocate(alphaB,cB)
 enddo


 call stop_clock(timing_overlap)


end subroutine setup_overlap_mixedbasis


!=========================================================================
subroutine setup_overlap_grad(basis,s_matrix_grad)
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: s_matrix_grad(basis%nbf,basis%nbf,3)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable        :: array_cart_gradx(:)
 real(C_DOUBLE),allocatable        :: array_cart_grady(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradz(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
!=====

 call start_clock(timing_overlap)
 write(stdout,'(/,a)') ' Setup gradient of the overlap matrix S (LIBINT)'

 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=1,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)

     allocate(array_cart_gradx(ni_cart*nj_cart))
     allocate(array_cart_grady(ni_cart*nj_cart))
     allocate(array_cart_gradz(ni_cart*nj_cart))

#if (LIBINT2_DERIV_ONEBODY_ORDER > 0)
     call libint_overlap_grad(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              array_cart_gradx,array_cart_grady,array_cart_gradz)
#else
     call die('overlap gradient not implemented without LIBINT one-body gradient terms')
#endif

     deallocate(alphaA,cA)

     ! X
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradx,matrix)
     s_matrix_grad(ibf1:ibf2,jbf1:jbf2,1) = matrix(:,:)

     ! Y
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_grady,matrix)
     s_matrix_grad(ibf1:ibf2,jbf1:jbf2,2) = matrix(:,:)

     ! Z
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradz,matrix)
     s_matrix_grad(ibf1:ibf2,jbf1:jbf2,3) = matrix(:,:)


     deallocate(array_cart_gradx)
     deallocate(array_cart_grady)
     deallocate(array_cart_gradz)
     deallocate(matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 title='=== Overlap matrix S (LIBINT) X ==='
 call dump_out_matrix(.FALSE.,title,s_matrix_grad(:,:,1))
 title='=== Overlap matrix S (LIBINT) Y ==='
 call dump_out_matrix(.FALSE.,title,s_matrix_grad(:,:,2))
 title='=== Overlap matrix S (LIBINT) Z ==='
 call dump_out_matrix(.FALSE.,title,s_matrix_grad(:,:,3))

 call stop_clock(timing_overlap)


end subroutine setup_overlap_grad


!=========================================================================
subroutine setup_kinetic(basis,hamiltonian_kinetic)
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic(basis%nbf,basis%nbf)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable :: array_cart(:)
 integer(C_INT)             :: amA,contrdepthA
 real(C_DOUBLE)             :: A(3)
 real(C_DOUBLE),allocatable :: alphaA(:)
 real(C_DOUBLE),allocatable :: cA(:)
 integer(C_INT)             :: amB,contrdepthB
 real(C_DOUBLE)             :: B(3)
 real(C_DOUBLE),allocatable :: alphaB(:)
 real(C_DOUBLE),allocatable :: cB(:)
!=====
 integer :: i_cart,j_cart,ij
 integer :: ibf_cart,jbf_cart
!=====

 call start_clock(timing_hamiltonian_kin)
#if defined(LIBINT2_SUPPORT_ONEBODY)
 write(stdout,'(/,a)') ' Setup kinetic part of the Hamiltonian (LIBINT)'
#else
 write(stdout,'(/,a)') ' Setup kinetic part of the Hamiltonian (internal)'
#endif


 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=jshell,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart(ni_cart*nj_cart))


#if defined(LIBINT2_SUPPORT_ONEBODY)
     call libint_kinetic(amA,contrdepthA,A,alphaA,cA, &
                         amB,contrdepthB,B,alphaB,cB, &
                         array_cart)
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#else
     ij = 0
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         ij = ij + 1
         ibf_cart = basis%shell(ishell)%istart_cart + i_cart - 1
         jbf_cart = basis%shell(jshell)%istart_cart + j_cart - 1
         call kinetic_basis_function(basis%bfc(ibf_cart),basis%bfc(jbf_cart),array_cart(ij))
       enddo
     enddo
     call transform_molgw_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#endif
     deallocate(alphaA,cA)



     hamiltonian_kinetic(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)
     hamiltonian_kinetic(jbf1:jbf2,ibf1:ibf2) = TRANSPOSE(matrix(:,:))


     deallocate(array_cart,matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 title='===  Kinetic energy contribution ==='
 call dump_out_matrix(.FALSE.,title,hamiltonian_kinetic)

 call stop_clock(timing_hamiltonian_kin)

end subroutine setup_kinetic


!=========================================================================
subroutine setup_kinetic_grad(basis,hamiltonian_kinetic_grad)
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic_grad(basis%nbf,basis%nbf,3)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable        :: array_cart_gradx(:)
 real(C_DOUBLE),allocatable        :: array_cart_grady(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradz(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
!=====

 call start_clock(timing_hamiltonian_kin)
 write(stdout,'(/,a)') ' Setup gradient of the kinetic part of the Hamiltonian (LIBINT)'

 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=1,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart_gradx(ni_cart*nj_cart))
     allocate(array_cart_grady(ni_cart*nj_cart))
     allocate(array_cart_gradz(ni_cart*nj_cart))

#if LIBINT2_DERIV_ONEBODY_ORDER > 0
     call libint_kinetic_grad(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              array_cart_gradx,array_cart_grady,array_cart_gradz)
#else
     call die('kinetic operator gradient not implemented without LIBINT one-body terms')
#endif
     deallocate(alphaA,cA)

     ! X
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradx,matrix)
     hamiltonian_kinetic_grad(ibf1:ibf2,jbf1:jbf2,1) = matrix(:,:)

     ! Y
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_grady,matrix)
     hamiltonian_kinetic_grad(ibf1:ibf2,jbf1:jbf2,2) = matrix(:,:)

     ! Z
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradz,matrix)
     hamiltonian_kinetic_grad(ibf1:ibf2,jbf1:jbf2,3) = matrix(:,:)


     deallocate(array_cart_gradx)
     deallocate(array_cart_grady)
     deallocate(array_cart_gradz)
     deallocate(matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 title='===  Kinetic energy contribution (LIBINT) X ==='
 call dump_out_matrix(.FALSE.,title,hamiltonian_kinetic_grad(:,:,1))
 title='===  Kinetic energy contribution (LIBINT) Y ==='
 call dump_out_matrix(.FALSE.,title,hamiltonian_kinetic_grad(:,:,2))
 title='===  Kinetic energy contribution (LIBINT) Z ==='
 call dump_out_matrix(.FALSE.,title,hamiltonian_kinetic_grad(:,:,3))

 call stop_clock(timing_hamiltonian_kin)

end subroutine setup_kinetic_grad


!=========================================================================
subroutine setup_nucleus(basis,hamiltonian_nucleus,atom_list)
 use m_atoms
 implicit none
 type(basis_set),intent(in)  :: basis
 real(dp),intent(out)        :: hamiltonian_nucleus(basis%nbf,basis%nbf)
 integer,intent(in),optional :: atom_list(:)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: icenter
 real(dp),allocatable :: matrix(:,:)
 real(C_DOUBLE),allocatable        :: array_cart(:)
 real(C_DOUBLE),allocatable        :: array_cart_C(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
 real(C_DOUBLE)                    :: C(3)
 integer  :: i_cart,j_cart,ij
 integer  :: ibf_cart,jbf_cart
 real(dp) :: nucleus
!=====

 if( in_tddft_loop ) then
   call start_clock(timing_tddft_hamiltonian_nuc)
 else
   call start_clock(timing_hamiltonian_nuc)
 end if

#if defined(LIBINT2_SUPPORT_ONEBODY)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian (LIBINT)'
#else
   write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian (internal)'
#endif

 if( PRESENT(atom_list) ) then
   write(stdout,'(1x,a,i5,a)') 'Only calculate the contribution from ',SIZE(atom_list),' nucleus/nuclei'
 endif

 hamiltonian_nucleus(:,:) = 0.0_dp

 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   if( MODULO(jshell-1,world%nproc) /= world%rank ) cycle

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   !$OMP PARALLEL PRIVATE(li,ni_cart,ni,ibf1,ibf2,amA,contrdepthA,A,alphaA,cA,array_cart,array_cart_C,C,matrix, &
   !$OMP&                 ij,ibf_cart,jbf_cart,nucleus)
   !$OMP DO
   do ishell=jshell,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart(ni_cart*nj_cart))
     allocate(array_cart_C(ni_cart*nj_cart))
     array_cart(:) = 0.0_dp

     do icenter=1,ncenter_nuclei
       ! Skip the contribution if icenter is not contained in the list
       if( PRESENT(atom_list) ) then
         if( ALL(atom_list(:) /= icenter ) ) cycle
       endif

       C(:) = xatom(:,icenter)
#if defined(LIBINT2_SUPPORT_ONEBODY)
       call libint_elecpot(amA,contrdepthA,A,alphaA,cA, &
                           amB,contrdepthB,B,alphaB,cB, &
                           C,array_cart_C)
       array_cart(:) = array_cart(:) - zvalence(icenter) * array_cart_C(:)
#else
       ij = 0
       do i_cart=1,ni_cart
         do j_cart=1,nj_cart
           ij = ij + 1
           ibf_cart = basis%shell(ishell)%istart_cart + i_cart - 1
           jbf_cart = basis%shell(jshell)%istart_cart + j_cart - 1
           call nucleus_basis_function(basis%bfc(ibf_cart),basis%bfc(jbf_cart),zvalence(icenter),xatom(:,icenter),nucleus)
           array_cart(ij) = array_cart(ij) + nucleus
         enddo
       enddo
#endif

     enddo
     deallocate(alphaA,cA)

#if defined(LIBINT2_SUPPORT_ONEBODY)
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#else
     call transform_molgw_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#endif

     hamiltonian_nucleus(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)
     hamiltonian_nucleus(jbf1:jbf2,ibf1:ibf2) = TRANSPOSE(matrix(:,:))


     deallocate(array_cart,array_cart_C,matrix)

   enddo
   !$OMP END DO
   !$OMP END PARALLEL
   deallocate(alphaB,cB)
 enddo

 !
 ! Reduce operation
 call world%sum(hamiltonian_nucleus)

 call dump_out_matrix(.FALSE.,'===  Nucleus potential contribution ===',hamiltonian_nucleus)

 if( in_tddft_loop ) then
   call stop_clock(timing_tddft_hamiltonian_nuc)
 else
   call stop_clock(timing_hamiltonian_nuc)
 endif

end subroutine setup_nucleus


!=========================================================================
subroutine setup_nucleus_grad(basis,hamiltonian_nucleus_grad)
 use m_atoms
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus_grad(basis%nbf,basis%nbf,ncenter_nuclei+1,3)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: natom_local
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: icenter
 character(len=100)   :: title
 real(dp),allocatable :: matrixA(:,:)
 real(dp),allocatable :: matrixB(:,:)

 real(C_DOUBLE),allocatable        :: array_cart_gradAx(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradAy(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradAz(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBx(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBy(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBz(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
 real(C_DOUBLE)                    :: C(3)
!=====

 if( ncenter_basis /= ncenter_nuclei ) call die('setup_nucleus_grad: not implemented with ghosts or projectiles')

 call start_clock(timing_hamiltonian_nuc)

 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian gradient (LIBINT)'
 if( world%nproc > 1 ) then
   natom_local=0
   do icenter=1,ncenter_nuclei
     if( world%rank /= MODULO(icenter-1,world%nproc) ) cycle
     natom_local = natom_local + 1
   enddo
   write(stdout,'(a)')         '   Parallelizing over atoms'
   write(stdout,'(a,i5,a,i5)') '   this proc treats ',natom_local,' over ',ncenter_nuclei
 endif

 hamiltonian_nucleus_grad(:,:,:,:) = 0.0_dp

 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=1,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)

     allocate(array_cart_gradAx(ni_cart*nj_cart))
     allocate(array_cart_gradAy(ni_cart*nj_cart))
     allocate(array_cart_gradAz(ni_cart*nj_cart))
     allocate(array_cart_gradBx(ni_cart*nj_cart))
     allocate(array_cart_gradBy(ni_cart*nj_cart))
     allocate(array_cart_gradBz(ni_cart*nj_cart))


     do icenter=1,ncenter_nuclei
       if( world%rank /= MODULO(icenter-1,world%nproc) ) cycle


       C(:) = xatom(:,icenter)

#if LIBINT2_DERIV_ONEBODY_ORDER > 0
       call libint_elecpot_grad(amA,contrdepthA,A,alphaA,cA, &
                                amB,contrdepthB,B,alphaB,cB, &
                                C,                           &
                                array_cart_gradAx,array_cart_gradAy,array_cart_gradAz, &
                                array_cart_gradBx,array_cart_gradBy,array_cart_gradBz)
#else
       call die('nuclear potential gradient not implemented without LIBINT one-body and gradient terms')
#endif
       array_cart_gradAx(:) = array_cart_gradAx(:) * (-zvalence(icenter))
       array_cart_gradAy(:) = array_cart_gradAy(:) * (-zvalence(icenter))
       array_cart_gradAz(:) = array_cart_gradAz(:) * (-zvalence(icenter))
       array_cart_gradBx(:) = array_cart_gradBx(:) * (-zvalence(icenter))
       array_cart_gradBy(:) = array_cart_gradBy(:) * (-zvalence(icenter))
       array_cart_gradBz(:) = array_cart_gradBz(:) * (-zvalence(icenter))

       ! X
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAx,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBx,matrixB)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,icenter  ,1) = -matrixA(:,:) - matrixB(:,:)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,ncenter_nuclei+1,1) = &
                  hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,ncenter_nuclei+1,1) + matrixA(:,:)

       ! Y
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAy,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBy,matrixB)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,icenter  ,2) = -matrixA(:,:) - matrixB(:,:)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,ncenter_nuclei+1,2) = &
                  hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,ncenter_nuclei+1,2) + matrixA(:,:)

       ! Z
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAz,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBz,matrixB)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,icenter  ,3) = -matrixA(:,:) - matrixB(:,:)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,ncenter_nuclei+1,3) = &
                  hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,ncenter_nuclei+1,3) + matrixA(:,:)

     enddo
     deallocate(alphaA,cA)


     deallocate(array_cart_gradAx)
     deallocate(array_cart_gradAy)
     deallocate(array_cart_gradAz)
     deallocate(array_cart_gradBx)
     deallocate(array_cart_gradBy)
     deallocate(array_cart_gradBz)
     deallocate(matrixA)
     deallocate(matrixB)

   enddo
   deallocate(alphaB,cB)

 enddo

 !
 ! Reduce operation
 call world%sum(hamiltonian_nucleus_grad)

 title='===  Nucleus potential contribution (LIBINT) C1X==='
 call dump_out_matrix(.FALSE.,title,hamiltonian_nucleus_grad(:,:,1,1))
 title='===  Nucleus potential contribution (LIBINT) C1Y==='
 call dump_out_matrix(.FALSE.,title,hamiltonian_nucleus_grad(:,:,1,2))
 title='===  Nucleus potential contribution (LIBINT) C1Z==='
 call dump_out_matrix(.FALSE.,title,hamiltonian_nucleus_grad(:,:,1,3))

 call stop_clock(timing_hamiltonian_nuc)

end subroutine setup_nucleus_grad


!=========================================================================
subroutine calculate_dipole_ao(basis,dipole_ao)
 implicit none
 type(basis_set),intent(in)         :: basis
 real(dp),allocatable,intent(out)   :: dipole_ao(:,:,:)
!=====
 integer              :: gt
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
 integer              :: li,lj,ni_cart,nj_cart,i_cart,j_cart
 integer              :: idir
 real(dp),allocatable :: dipole_cart(:,:,:)
#if defined(HAVE_LIBCINT)
 integer :: info
 integer :: shls(2)
 real(dp),allocatable :: array_cart(:,:)
 real(dp),allocatable :: matrix(:,:)
#endif
!=====

#if defined(HAVE_LIBCINT)
 write(stdout,'(/,a)') ' Setup dipole matrix (LIBCINT)'
#else
 write(stdout,'(/,a)') ' Setup dipole matrix (internal)'
#endif
 gt = get_gaussian_type_tag(basis%gaussian_type)

 allocate(dipole_ao(basis%nbf,basis%nbf,3))


 do jshell=1,basis%nshell
   lj        = basis%shell(jshell)%am
   nj_cart   = number_basis_function_am('CART',lj)
   jbf1      = basis%shell(jshell)%istart
   jbf1_cart = basis%shell(jshell)%istart_cart
   jbf2      = basis%shell(jshell)%iend

   do ishell=1,basis%nshell
     li        = basis%shell(ishell)%am
     ni_cart   = number_basis_function_am('CART',li)
     ibf1      = basis%shell(ishell)%istart
     ibf1_cart = basis%shell(ishell)%istart_cart
     ibf2      = basis%shell(ishell)%iend


#if defined(HAVE_LIBCINT)
     allocate(array_cart(ni_cart*nj_cart,3))
     shls(1) = ishell-1  ! C convention starts with 0
     shls(2) = jshell-1  ! C convention starts with 0
     info = cint1e_r_cart(array_cart, shls, atm, LIBCINT_natm, bas, LIBCINT_nbas, env)

     do idir=1,3
       call transform_libcint_to_molgw(basis%gaussian_type,li,lj,array_cart(:,idir),matrix)
       dipole_ao(ibf1:ibf2,jbf1:jbf2,idir) = matrix(:,:)
       dipole_ao(jbf1:jbf2,ibf1:ibf2,idir) = TRANSPOSE(matrix(:,:))
     enddo
     deallocate(matrix)

     deallocate(array_cart)
#else
     allocate(dipole_cart(3,ni_cart,nj_cart))

     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call basis_function_dipole(basis%bfc(ibf1_cart+i_cart-1),basis%bfc(jbf1_cart+j_cart-1),dipole_cart(:,i_cart,j_cart))
       enddo
     enddo

     do idir=1,3
       dipole_ao(ibf1:ibf2,jbf1:jbf2,idir) = MATMUL( TRANSPOSE( cart_to_pure(li,gt)%matrix(:,:) ) , &
             MATMUL(  dipole_cart(idir,:,:) , cart_to_pure(lj,gt)%matrix(:,:) ) )
     enddo

     deallocate(dipole_cart)
#endif

   enddo
 enddo

 call dump_out_matrix(.FALSE.,'===  Dipole AO X ===',dipole_ao(:,:,1))
 call dump_out_matrix(.FALSE.,'===  Dipole AO Y ===',dipole_ao(:,:,2))
 call dump_out_matrix(.FALSE.,'===  Dipole AO Z ===',dipole_ao(:,:,3))


end subroutine calculate_dipole_ao


!=========================================================================
subroutine calculate_quadrupole_ao(basis,quadrupole_ao)
 implicit none
 type(basis_set),intent(in)         :: basis
 real(dp),allocatable,intent(out)   :: quadrupole_ao(:,:,:,:)
!=====
 integer              :: gt
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
 integer              :: li,lj,ni_cart,nj_cart,i_cart,j_cart
 integer              :: idir,jdir,ijdir
 real(dp),allocatable :: quadrupole_cart(:,:,:,:)
#if defined(HAVE_LIBCINT)
 integer :: info
 integer :: shls(2)
 real(dp),allocatable :: array_cart(:,:)
 real(dp),allocatable :: matrix(:,:)
#endif
!=====

#if defined(HAVE_LIBCINT)
 write(stdout,'(/,a)') ' Setup quadrupole matrix (LIBCINT)'
#else
 write(stdout,'(/,a)') ' Setup quadrupole matrix (internal)'
#endif
 gt = get_gaussian_type_tag(basis%gaussian_type)

 allocate(quadrupole_ao(basis%nbf,basis%nbf,3,3))


 do jshell=1,basis%nshell
   lj        = basis%shell(jshell)%am
   nj_cart   = number_basis_function_am('CART',lj)
   jbf1      = basis%shell(jshell)%istart
   jbf1_cart = basis%shell(jshell)%istart_cart
   jbf2      = basis%shell(jshell)%iend

   do ishell=1,basis%nshell
     li        = basis%shell(ishell)%am
     ni_cart   = number_basis_function_am('CART',li)
     ibf1      = basis%shell(ishell)%istart
     ibf1_cart = basis%shell(ishell)%istart_cart
     ibf2      = basis%shell(ishell)%iend


#if defined(HAVE_LIBCINT)
     allocate(array_cart(ni_cart*nj_cart,9))
     shls(1) = ishell-1  ! C convention starts with 0
     shls(2) = jshell-1  ! C convention starts with 0
     info = cint1e_rr_cart(array_cart, shls, atm, LIBCINT_natm, bas, LIBCINT_nbas, env)

     ijdir=0
     do jdir=1,3
       do idir=1,3
         ijdir=ijdir+1
         call transform_libcint_to_molgw(basis%gaussian_type,li,lj,array_cart(:,ijdir),matrix)
         quadrupole_ao(ibf1:ibf2,jbf1:jbf2,idir,jdir) = matrix(:,:)
         quadrupole_ao(jbf1:jbf2,ibf1:ibf2,idir,jdir) = TRANSPOSE(matrix(:,:))
       enddo
     enddo
     deallocate(matrix)

     deallocate(array_cart)
#else
     allocate(quadrupole_cart(ni_cart,nj_cart,3,3))

     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call basis_function_quadrupole(basis%bfc(ibf1_cart+i_cart-1),basis%bfc(jbf1_cart+j_cart-1), &
                                        quadrupole_cart(i_cart,j_cart,:,:))
       enddo
     enddo

     do jdir=1,3
       do idir=1,3
         quadrupole_ao(ibf1:ibf2,jbf1:jbf2,idir,jdir) = MATMUL( TRANSPOSE( cart_to_pure(li,gt)%matrix(:,:) ) , &
               MATMUL(  quadrupole_cart(:,:,idir,jdir) , cart_to_pure(lj,gt)%matrix(:,:) ) )
       enddo
     enddo
     deallocate(quadrupole_cart)
#endif


   enddo
 enddo


end subroutine calculate_quadrupole_ao


!=========================================================================
subroutine calculate_gos_ao(basis,qvec,gos_ao)
 implicit none
 type(basis_set),intent(in)          :: basis
 real(dp),intent(in)                 :: qvec(3)
 complex(dp),allocatable,intent(out) :: gos_ao(:,:)
!=====
 integer              :: gt
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
 integer              :: li,lj,ni_cart,nj_cart,i_cart,j_cart
 complex(dp),allocatable :: gos_cart(:,:)
!=====

 gt = get_gaussian_type_tag(basis%gaussian_type)

 allocate(gos_ao(basis%nbf,basis%nbf))

 !$OMP PARALLEL PRIVATE(li,lj,ni_cart,nj_cart,ibf1,ibf1_cart,ibf2,jbf1,jbf1_cart,jbf2,gos_cart)
 !$OMP DO
 do jshell=1,basis%nshell
   lj        = basis%shell(jshell)%am
   nj_cart   = number_basis_function_am('CART',lj)
   jbf1      = basis%shell(jshell)%istart
   jbf1_cart = basis%shell(jshell)%istart_cart
   jbf2      = basis%shell(jshell)%iend

   do ishell=1,basis%nshell
     li        = basis%shell(ishell)%am
     ni_cart   = number_basis_function_am('CART',li)
     ibf1      = basis%shell(ishell)%istart
     ibf1_cart = basis%shell(ishell)%istart_cart
     ibf2      = basis%shell(ishell)%iend


     allocate(gos_cart(ni_cart,nj_cart))

     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call basis_function_gos(basis%bfc(ibf1_cart+i_cart-1),basis%bfc(jbf1_cart+j_cart-1),qvec,gos_cart(i_cart,j_cart))
       enddo
     enddo

     gos_ao(ibf1:ibf2,jbf1:jbf2) = MATMUL( TRANSPOSE( cart_to_pure(li,gt)%matrix(:,:) ) , &
             MATMUL(  gos_cart(:,:) , cart_to_pure(lj,gt)%matrix(:,:) ) )

     deallocate(gos_cart)

   enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL


end subroutine calculate_gos_ao


!=========================================================================
subroutine setup_nucleus_ecp(basis,hamiltonian_nucleus)
 use m_atoms
 use m_dft_grid
 use m_ecp
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(inout)     :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!=====
 integer              :: ibf,jbf
 integer              :: icenter

 integer              :: iecp
 integer              :: iproj,nproj
 integer              :: mm
 real(dp)             :: rr(3)
 real(dp)             :: basis_function_r(basis%nbf)
 integer              :: iradial
 integer              :: i1,n1
 real(dp)             :: xtmp,phi,cos_theta
 real(dp)             :: wxa(nradial_ecp),xa(nradial_ecp)
 real(dp)             :: w1(nangular_ecp),x1(nangular_ecp),y1(nangular_ecp),z1(nangular_ecp)
 real(dp),allocatable :: int_fixed_r(:,:)
 real(dp),external    :: real_spherical_harmonics
 integer              :: necp,ie
 character(len=100)   :: title
 logical              :: element_has_ecp
 real(dp)             :: r1,r2
 real(dp),allocatable :: vr(:),ur(:)
 real(dp),allocatable :: kb(:,:)
 integer              :: ir
!=====

 ! Check if there are some ECP
 if( nelement_ecp == 0 ) return


 call start_clock(timing_ecp)

 !
 ! Since there will be an allreduce operation in the end,
 ! anticipate by dividing the input value of Hnucl by the number of procs
 if( world%nproc > 1 ) then
   hamiltonian_nucleus(:,:) = hamiltonian_nucleus(:,:) / world%nproc
 endif

 n1 = nangular_ecp
 select case(nangular_ecp)
 case(6)
   call ld0006(x1,y1,z1,w1,n1)
 case(14)
   call ld0014(x1,y1,z1,w1,n1)
 case(26)
   call ld0026(x1,y1,z1,w1,n1)
 case(38)
   call ld0038(x1,y1,z1,w1,n1)
 case(50)
   call ld0050(x1,y1,z1,w1,n1)
 case(74)
   call ld0074(x1,y1,z1,w1,n1)
 case(86)
   call ld0086(x1,y1,z1,w1,n1)
 case(110)
   call ld0110(x1,y1,z1,w1,n1)
 case(146)
   call ld0146(x1,y1,z1,w1,n1)
 case(170)
   call ld0170(x1,y1,z1,w1,n1)
 case(230)
   call ld0230(x1,y1,z1,w1,n1)
 case(302)
   call ld0302(x1,y1,z1,w1,n1)
 case(434)
   call ld0434(x1,y1,z1,w1,n1)
 case default
   write(stdout,*) 'grid points: ',nangular_ecp
   call die('setup_nucleus_ecp: Lebedev grid is not available')
 end select






 do icenter=1,ncenter_nuclei
   element_has_ecp = .FALSE.
   do ie=1,nelement_ecp
     if( ABS( element_ecp(ie) - zatom(icenter) ) < 1.0e-5_dp ) then
       element_has_ecp = .TRUE.
       exit
     endif
   enddo

   if( .NOT. element_has_ecp ) cycle

   necp = ecp(ie)%necp

   select case(ecp(ie)%ecp_format)
   case(ECP_PSP6,ECP_PSP8)
     do iradial=1,nradial_ecp
       xa(iradial)  = ecp(ie)%rad(iradial)
       if(iradial == 1) then
         wxa(iradial) = ecp(ie)%rad(iradial+1)-ecp(ie)%rad(iradial)
       else
         wxa(iradial) = 0.5_dp * ( ecp(ie)%rad(iradial+1)-ecp(ie)%rad(iradial-1) )
       endif
     enddo
   case default
     do iradial=1,nradial_ecp
       xtmp = ( iradial - 0.5_dp ) / REAL(nradial_ecp,dp)
       xa(iradial)   = -5.0_dp * log( 1.0_dp - xtmp**3)
       wxa(iradial)  = 3.0_dp * 5.0_dp * xtmp**2 / ( 1.0_dp - xtmp**3 ) / REAL(nradial_ecp,dp)
     enddo
   end select

   nproj = 0
   do iecp=1,necp
     if( ecp(ie)%lk(iecp) /= -1 ) then   ! -1 encodes a local component
       nproj = nproj + number_basis_function_am('PURE',ecp(ie)%lk(iecp))
     endif
   enddo
   select case(ecp(ie)%ecp_format)
   case(ECP_PSP6)
     allocate(vr(necp),ur(necp))
     allocate(kb(basis%nbf,nproj))
     kb(:,:) = 0.0_dp
   case(ECP_PSP8)
     allocate(vr(necp))
     allocate(kb(basis%nbf,nproj))
     kb(:,:) = 0.0_dp
   case(ECP_NWCHEM)
     allocate(int_fixed_r(basis%nbf,nproj))
   end select


   do iradial=1,nradial_ecp
     if( MODULO(iradial-1,world%nproc) /= world%rank ) cycle

     ! avoid x=0, divergence in the local component and 0 weight anyway
     if( ABS(xa(iradial)) < 1.0e-10_dp ) cycle

     select case(ecp(ie)%ecp_format)
     case(ECP_PSP6,ECP_PSP8)
        vr(:) = ecp(ie)%vpspll(iradial,:)
        if( ALLOCATED(ur) ) then
          ur(:) = ecp(ie)%wfll(iradial,:)
        endif
       ! remove the Coulomb part for the local part: it is treated in the regular routine setup_nucleus
       do iecp=1,necp
         if( ecp(ie)%lk(iecp) == -1 ) then   ! -1 encodes a local component
           vr(iecp) = vr(iecp) + zvalence(icenter) /  xa(iradial)
         endif
       enddo
     end select

     if( ALLOCATED(int_fixed_r) ) int_fixed_r(:,:) = 0.0_dp
     do i1=1,nangular_ecp
       rr(1) = xa(iradial) * x1(i1) + xatom(1,icenter)
       rr(2) = xa(iradial) * y1(i1) + xatom(2,icenter)
       rr(3) = xa(iradial) * z1(i1) + xatom(3,icenter)
       call calculate_basis_functions_r(basis,rr,basis_function_r)

       cos_theta = z1(i1)
       phi       = ATAN2(y1(i1),x1(i1))

       iproj = 0
       do iecp=1,necp

         !
         ! Treat the local potential:
         !
         if( ecp(ie)%lk(iecp) == -1 ) then
           select case(ecp(ie)%ecp_format)
           case(ECP_NWCHEM)
             do jbf=1,basis%nbf
               do ibf=1,basis%nbf
                 hamiltonian_nucleus(ibf,jbf) = hamiltonian_nucleus(ibf,jbf)  &
                     + basis_function_r(ibf) * basis_function_r(jbf) * w1(i1) * 4.0_dp * pi  &
                        * wxa(iradial) * xa(iradial)**2  &
                        * ecp(ie)%dk(iecp) * EXP( -ecp(ie)%zetak(iecp) * xa(iradial)**2 ) * xa(iradial)**(ecp(ie)%nk(iecp)-2)
               enddo
             enddo
           case(ECP_PSP6,ECP_PSP8)
             do jbf=1,basis%nbf
               do ibf=1,basis%nbf
                 hamiltonian_nucleus(ibf,jbf) = hamiltonian_nucleus(ibf,jbf)  &
                     + basis_function_r(ibf) * basis_function_r(jbf) * w1(i1) * 4.0_dp * pi  &
                        * wxa(iradial) * xa(iradial)**2  &
                        * vr(iecp)
               enddo
             enddo
           end select

         else
           !
           ! non local case
           !
           select case(ecp(ie)%ecp_format)
           case(ECP_NWCHEM)
             do mm=-ecp(ie)%lk(iecp),ecp(ie)%lk(iecp)
               iproj = iproj + 1
               int_fixed_r(:,iproj) = int_fixed_r(:,iproj) + basis_function_r(:) &
                                         * real_spherical_harmonics(ecp(ie)%lk(iecp),mm,cos_theta,phi) &
                                            * w1(i1) * 4.0_dp * pi
             enddo
           case(ECP_PSP6)
             do mm=-ecp(ie)%lk(iecp),ecp(ie)%lk(iecp)
               iproj = iproj + 1
               kb(:,iproj) = kb(:,iproj) + basis_function_r(:) &
                                 * real_spherical_harmonics(ecp(ie)%lk(iecp),mm,cos_theta,phi) &
                                    * w1(i1) * 4.0_dp * pi *  wxa(iradial) * xa(iradial) &
                                        * vr(iecp) * ur(iecp)
             enddo
           case(ECP_PSP8)
             do mm=-ecp(ie)%lk(iecp),ecp(ie)%lk(iecp)
               iproj = iproj + 1
               kb(:,iproj) = kb(:,iproj) + basis_function_r(:) &
                                 * real_spherical_harmonics(ecp(ie)%lk(iecp),mm,cos_theta,phi) &
                                    * w1(i1) * 4.0_dp * pi *  wxa(iradial) * xa(iradial) &
                                        * vr(iecp)
             enddo
           end select
         endif

       enddo ! iecp
     enddo ! (theta, phi) points

     ! non-local ECP contribution
     if( ecp(ie)%ecp_format == ECP_NWCHEM ) then
       iproj = 0
       do iecp=1,necp
         if( ecp(ie)%lk(iecp) /= -1 ) then
           do mm=-ecp(ie)%lk(iecp),ecp(ie)%lk(iecp)
             iproj = iproj + 1
             do jbf=1,basis%nbf
               do ibf=1,basis%nbf
                 hamiltonian_nucleus(ibf,jbf) = hamiltonian_nucleus(ibf,jbf)  &
                     + int_fixed_r(ibf,iproj) * int_fixed_r(jbf,iproj) * wxa(iradial) * xa(iradial)**2  &
                        * ecp(ie)%dk(iecp) * EXP( -ecp(ie)%zetak(iecp) * xa(iradial)**2 ) * xa(iradial)**(ecp(ie)%nk(iecp)-2)
               enddo
             enddo
           enddo
         endif
       enddo
     endif

   enddo  ! iradial


   ! non-local numerical grid contribution
   select case(ecp(ie)%ecp_format)
   case(ECP_PSP6,ECP_PSP8)
     iproj = 0
     do iecp=1,necp
       if( ecp(ie)%lk(iecp) /= -1 ) then
         do mm=-ecp(ie)%lk(iecp),ecp(ie)%lk(iecp)
           iproj = iproj + 1
           do jbf=1,basis%nbf
             do ibf=1,basis%nbf
               hamiltonian_nucleus(ibf,jbf) =  hamiltonian_nucleus(ibf,jbf)  &
                   + kb(ibf,iproj) * kb(jbf,iproj) * ecp(ie)%ekb(iecp)
             enddo
           enddo
         enddo
       endif
     enddo
   end select


   if( ALLOCATED(kb) ) deallocate(kb)
   if( ALLOCATED(int_fixed_r) ) deallocate(int_fixed_r)
   if( ALLOCATED(ur) ) deallocate(ur)
   if( ALLOCATED(vr) ) deallocate(vr)

 enddo ! ie

 call world%sum(hamiltonian_nucleus)

 title='=== ECP Nucleus potential contribution ==='
 call dump_out_matrix(.FALSE.,title,hamiltonian_nucleus)

 call stop_clock(timing_ecp)

end subroutine setup_nucleus_ecp


end module m_hamiltonian_onebody
!=========================================================================
